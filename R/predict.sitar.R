#' Predict SITAR model
#'
#' Predict method for \code{sitar} objects, based on \code{predict.lme}.
#'
#' When \code{deriv = 1} the returned velocity is in units of \code{yfun(y)}
#' per \code{xfun(x)} where \code{xfun} and \code{yfun} default to
#' \code{identity}. If \code{x} and/or \code{y} include transformations,
#' e.g. \code{x = log(age)}, velocity
#' in the original units is obtained by specifying \code{xfun}
#' and/or \code{yfun} to back-transform them appropriately. By default this
#' is done automatically by \code{ifun}, so if for example \code{x = log(age)}
#' then \code{xfun} is set to \code{exp}. In this way velocity is in the
#' untransformed units by default.
#'
#' @param object an object inheriting from class \code{sitar}.
#' @param newdata an optional data frame to be used for obtaining the
#' predictions, defaulting to the data used to fit \code{object}.
#' It requires named columns for \code{x}, and also for \code{id} if
#' \code{level = 1}, matching the names in \code{object}. Variables with the
#' reserved names \code{x = .x} or \code{id = .id} take precedence over the model
#' \code{x} and \code{id} variables. Any covariates in
#' \code{a.formula}, \code{b.formula}, \code{c.formula} or \code{d.formula} can also be included.
#' By default their values are set to the mean, so when \code{level = 0} the
#' prediction represents the mean curve.
#' @param level an optional integer vector giving the level(s) of grouping to be used
#' in obtaining the predictions, level 0 corresponding to the population
#' predictions. Defaults to level 1, and \code{level = 0:1} fits both levels.
#' @param \dots other optional arguments: \code{asList}, \code{na.action} and
#' \code{naPattern}.
#' @param deriv an optional integer specifying predictions corresponding to
#' the fitted curve and/or its derivative. \code{deriv = 0} (default)
#' specifies the distance curve, \code{deriv = 1} the velocity curve and
#' \code{deriv = 0:1} fits both.
#' @param abc an optional named vector containing values of a subset of
#' \code{a}, \code{b}, \code{c} and \code{d}, default \code{NULL}. Ignored if
#' \code{level = 0}. It gives predictions for a single subject with the
#' specified values of \code{a}, \code{b}, \code{c} and \code{d}, where missing values
#' are set to 0. Alternatively \code{abc} can contain the value for a single id.
#' @param xfun an optional function to apply to \code{x} to convert it back to
#' the original scale, e.g. if x = log(age) then xfun = exp. Only relevant if
#' \code{deriv == 1} - see Details.
#' @param yfun an optional function to apply to \code{y} to convert it back to
#' the original scale, e.g. if y = sqrt(height) then yfun = function(z) z^2.
#' @return A vector of the predictions, or a list of vectors if \code{asList =
#' TRUE} and \code{level == 1}, or a data frame if \code{length(level) * length(deriv) > 1}.
#' The data frame column names are (a subset of): \code{(id, predict.fixed,
#' predict.id, vel.predict.fixed, vel.predict.id)}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{ifun}} for a way to generate the functions \code{xfun}
#' and \code{yfun} automatically from the \code{sitar} model call.
#' @examples
#'
#' data(heights)
#' ## fit model
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)
#'
#' ## predictions at level 0
#' predict(m1, newdata=data.frame(age=9:16), level=0)
#'
#' ## predictions at level 1 for subject 5
#' predict(m1, newdata=data.frame(age=9:16, id=5), level=1)
#'
#' ## velocity predictions for subjects with early and late puberty
#' vel1 <- predict(m1, deriv=1, abc=c(b=-1))
#' mplot(age, vel1, id, heights, col=id)
#' vel1 <- predict(m1, deriv=1, abc=c(b=1))
#' mplot(age, vel1, id, heights, col=id, add=TRUE)
#'
#' @importFrom rlang .data
#' @importFrom dplyr arrange mutate n nest_by pull
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr unnest
#' @export
  predict.sitar <- function(object, newdata = getData(object), level = 1L, ...,
                            deriv = 0L, abc = NULL,
                            xfun = identity, yfun = identity) {
# obtain distance and velocity predictions
    predictions <- function(object, newdata, level, dx = 1e-5, ...) {
      xfun <- ifun(object$call.sitar$x)
      yfun <- ifun(object$call.sitar$y)
      # offset for mean curve
      xy.id <- xyadj(object, x = x, id = id, abc = re.mean)
      # convert object from sitar to nlme
      class(object) <- setdiff(class(object), 'sitar')
      # calculate predictions by level
      map_dfr(setNames(level, c('fixed', 'id')[level + 1L]), \(ilevel){
        newdata %>%
          rownames_to_column('rowname') %>%
          as_tibble %>%
          mutate(x = if (ilevel == 0L) xy.id$x - xoffset else x,
                 xc = xfun(x + xoffset),
                 y = predict(object, ., level = ilevel),
                 ylo = predict(object, . |> mutate(x = x - dx), level = ilevel),
                 yhi = predict(object, . |> mutate(x = x + dx), level = ilevel),
                 y = if (ilevel == 0L) .data$y - xy.id$y else .data$y,
                 predict = yfun(.data$y),
                 # calculate velocity as dy/dx
                 vel.predict = (.data$yhi - .data$ylo) / dx / 2 /
                   Dxy(object, .data$predict, 'y') * Dxy(object, .data$xc, 'x')) %>%
          select(c(id, .data$rowname, predict, .data$vel.predict))
      }, .id = 'level') %>%
        pivot_wider(names_from = 'level', values_from = c(predict, .data$vel.predict),
                    names_sep = '.') %>%
        mutate(across(ends_with('.fixed'), unname),
               across(ends_with('.id'), ~setNames(., id))) %>%
        select(-.data$rowname)
    }
# apply differential of x or y to variable
    Dxy <- function(object, var, which = c('x', 'y')) {
      which <- match.arg(which)
      call <- object$call.sitar[[which]]
      name <- all.vars(call)
      eval(D(call, name), setNames(as.data.frame(var), name))
    }
    mc <- match.call()
# match call for sitar
    oc <- object$call.sitar
# random effects
    re <- ranef(object)
# ensure level is subset of 0:1
    level <- intersect(as.integer(level), 0L:1L)
# ensure deriv is subset of 0:1
    deriv <- intersect(as.integer(deriv), 0L:1L)
# check if old-style object lacking fitnlme
    if (!'fitnlme' %in% names(object)) {
      warning('fitnlme missing - best to refit model')
      object <- update(object, control=nlmeControl(maxIter=0, pnlsMaxIter=0, msMaxIter=0L))
    }
# attach object for fitnlme
    on.exit(detach(object))
    eval(parse(text='attach(object)'))
# identify sitar formula covariates
    covnames <- all.vars(asOneFormula(oc$a.formula, oc$b.formula, oc$c.formula, oc$d.formula))
# covariates also in newdata
    yesnames <- intersect(covnames, names(newdata))
# covariates not in newdata
    notnames <- setdiff(covnames, yesnames)
    newdata[, notnames] <- NA
# if groups add linear contrasts to newdata
    if (length(yesnames) > 0L) {
      groupnames <- yesnames[sapply(yesnames, function(x) !is.numeric(getData(object)[[x]]))]
      if (length(groupnames) > 0L) {
# convert newdata to factors
        for (i in groupnames) {
          storage.mode(newdata[[i]]) <- storage.mode(getData(object)[[i]])
          attributes(newdata[[i]]) <- attributes(getData(object)[[i]])
        }
        extra <- formula(paste(c("~1", groupnames), collapse = "+"))
        extra <- as_tibble(model.matrix(extra, newdata))[-1]
# ensure valid names to match sitar model names
        names(extra) <- make.names(names(extra), unique = TRUE)
        newdata <- bind_cols(newdata, extra)
      }
# centre covariates in newdata (using means from sitar)
      gd <- update(object, returndata=TRUE)
      covmeans <- attr(gd, 'scaled:center')
# covariates and contrasts in model
      allnames <- setdiff(formalArgs(object$fitnlme), names(coef(object)))[-1]
      for (i in allnames)
        if (is.numeric(newdata[, i]))
          newdata[, i] <- newdata[, i] - covmeans[i]
    }
# set covariates not in newdata to 0
    if (length(notnames) > 0L) {
# check if subset from plot
      subset <- attr(newdata, 'subset')
# centre covariates not in newdata to mean gd
      if (is.null(subset)) {
        newdata[, notnames] <- 0
      } else {
        if (!exists('gd'))
          gd <- update(object, returndata=TRUE)
        for (i in notnames)
          newdata[, i] <- mean(gd[subset, i])
      }
    }
# centre random effects for level 0
    re <- scale(re, scale = FALSE)
    re.mean <- data.frame(t(attr(re, 'scaled:center')))
    re.mean <- re.mean[rep(1, nrow(newdata)), , drop = FALSE]
# create x in newdata
    x <- if ('.x' %in% names(newdata))
      newdata$.x
    else
      eval(oc$x, newdata)
    if (is.null(xoffset <- object$xoffset)) {
      xoffset <- mean(getCovariate(object))
      warning('xoffset set to mean(x) - best to refit model')
    }
    newdata$x <- x - xoffset
# check abc and set level accordingly
    if (!is.null(abc)) {
      if (is.null(names(abc)) && length(abc) == 1 &&
          as.character(abc) %in% rownames(re))
      # abc is id
        level <- 1L
      else if (!is.null(names(abc)) && all(names(abc) %in% names(re.mean))) {
        # abc is offset from re.mean
        level <- 0L
        for (i in names(abc))
          re.mean[, i] <- re.mean[, i] + abc[i]
        abc <- NULL
      } else
        stop('abc unrecognised')
    }
# create id in newdata
    newdata$id <- if ('.id' %in% names(newdata))
      newdata$.id
    else
      if (any(level == 1L)) {
        if (!is.null(abc))
          factor(1, labels = abc)
        else
          eval(oc$id, newdata)
      } else
        factor(1, labels = rownames(re)[1])
    id <- newdata$id
# predictions
    output <- predictions(object, newdata, level = level)
# return columns of data frame
    # single column
    if (length(level) * length(deriv) == 1L) {
      output <- output %>%
        pull(deriv + 2L)
    # split by id?
      asList <- ifelse(is.null(asList <- eval(mc$asList)), FALSE, asList)
      if (asList && level == 1L)
        output <- split(output, id)
    }
    # multiple columns
    else if (length(level) == 2L && length(deriv) == 1L) {
      output <- output %>%
      select(1L, 2L:3L + deriv * 2L)
    }
    return(output)
  }
