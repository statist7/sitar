#' Predict SITAR model
#'
#' Predict method for \code{sitar} objects, based on \code{predict.lme}.
#'
#' When \code{deriv = 1} the returned velocity is in units of \code{yfun(y)}
#' per \code{xfun(x)}. So if \code{x} and/or \code{y} are transformed, velocity
#' in units of \code{y} per \code{x} can be obtained by specifying \code{xfun}
#' and/or \code{yfun} to back-transform them appropriately.
#'
#' @param object an object inheriting from class \code{sitar}.
#' @param newdata an optional data frame to be used for obtaining the
#' predictions, defaulting to the data used to fit \code{object}.
#' It requires named columns for \code{x}, and for \code{id} if
#' \code{level = 1}, matching the names in \code{object}. Variables with the
#' reserved names \code{x=.x} or \code{id=.id} take precedence over the model
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
#' either the fitted curve or its derivative. \code{deriv = 0} (default)
#' specifies the distance curve, \code{deriv = 1} the velocity curve and
#' \code{deriv = 2} the acceleration curve.
#' @param abc an optional named vector containing values of a subset of
#' \code{a}, \code{b}, \code{c} and \code{d}, default \code{NULL}. Ignored if
#' \code{level = 0}. It gives predictions for a single subject with the
#' specified values of \code{a}, \code{b}, \code{c} and \code{d}, where missing values
#' are set to 0. Alternatively \code{abc} can contain the value for a single id.
#' @param xfun an optional function to apply to \code{x} to convert it back to
#' the original scale, e.g. if x = log(age) then xfun = exp. Only relevant if
#' \code{deriv > 0} - see Details.
#' @param yfun an optional function to apply to \code{y} to convert it back to
#' the original scale, e.g. if y = sqrt(height) then yfun = function(z) z^2.
#' @return A vector of the predictions, or a list of vectors if \code{asList =
#' TRUE} and \code{level == 1}, or a data frame if \code{length(level) > 1}.
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
#' @importFrom tidyr unnest
#' @export
  predict.sitar <- function(object, newdata=getData(object), level=1L, ...,
                            deriv=0L, abc=NULL,
                            xfun=identity, yfun=identity) {

# derive x-y velocity curve
    get_vel <- function(x, y, deriv = 1) {
      smooth.spline(x, y) |>
        predict(x, deriv = deriv)
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
# ensure level is integral
    level <- as.integer(level)
# ensure deriv is integral and unique
    deriv <- as.integer(max(deriv))
# check if old-style object lacking fitnlme
    if (!'fitnlme' %in% names(object)) {
      warning('fitnlme missing - best to refit model')
      object <- update(object, control=nlmeControl(maxIter=0, pnlsMaxIter=0, msMaxIter=0))
    }
# attach object for fitnlme
    on.exit(detach(object))
    eval(parse(text='attach(object)'))
# identify sitar formula covariates in newdata
    covnames <- all.vars(asOneFormula(oc$a.formula, oc$b.formula, oc$c.formula, oc$d.formula))
    covnames <- covnames[covnames %in% names(newdata)]
# if non-numeric add linear contrasts to newdata
    factornames <- covnames[unlist(lapply(covnames, function(x) !is.numeric(with(newdata, get(x)))))]
    if (length(factornames) > 0L) {
      extra <- eval(parse(text = paste(c("~0", factornames), collapse = "+"))[[1]])
      extra <- as_tibble(model.matrix(extra, newdata))
      newdata <- bind_cols(newdata, extra)
      covnames <- c(covnames, names(extra))
    }
# identify covariates in model (not x or coef)
    argnames <- names(formals(fitnlme))
    argnames <- argnames[!argnames %in% names(coef(object))][-1]
    if (length(argnames) > 0) {
# drop any factors in covnames
      covnames <- names(newdata)
      covnames <- covnames[covnames %in% argnames]
# set to 0 covariates not in newdata
      notnames <- argnames[!argnames %in% covnames]
      newdata[, notnames] <- 0
# centre covariates in newdata (using means from sitar)
      if (length(covnames) > 0) {
        gd <- update(object, returndata=TRUE)
        covmeans <- attr(gd, 'scaled:center')
        for (i in covnames)
          newdata[, i] <- newdata[, i] - covmeans[i]
      }
    }
# check if subset from plot
    subset <- attr(newdata, 'subset')
# centre covariates not in newdata to mean gd
    if (!is.null(subset)) {
      if (exists('notnames') && length(notnames) > 0) {
        if (!exists('gd'))
          gd <- update(object, returndata=TRUE)
        for (i in notnames)
          newdata[, i] <- mean(gd[subset, i])
      }
    }
# centre random effects
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
# create id in newdata
    id <- rownames(re)
    newdata$id <- if ('.id' %in% names(newdata))
      newdata$.id
    else {
      if (any(level == 1L) && is.null(abc))
        eval(oc$id, newdata)
      else
        factor(1, labels=id[1])
    }
# check if abc is id
    if (all(level == 0L))
      abc <- NULL
    else if (!is.null(abc) && is.null(names(abc)) &&
             length(abc) == 1 && as.character(abc) %in% rownames(re)) {
      newdata$id <- abc
      abc <- NULL
    }
    id <- newdata$id
# adjust abc for mean ranef
    if (!is.null(abc)) {
      abc.t <- re.mean
      for (i in names(abc.t))
        if (!is.na(abc[i]))
          abc.t[, i] <- abc.t[, i] + abc[i]
      abc <- abc.t
    }
# create nlme object
    object.nlme <- structure(object, class = c('nlme', 'lme'))
# DISTANCE
# level 1 prediction
    if (is.null(abc))
      pred <- predict(object.nlme, newdata)
    else {
      xy.id <- xyadj(object, x = x, id = id, abc = abc)
      newdata$x <- xy.id$x - xoffset
      pred <- predict(object.nlme, newdata, level = 0L) - xy.id$y
    }
    pred <- yfun(pred)
# level 0 prediction
    xy.id <- xyadj(object, x = x, id = id, abc = re.mean)
    newdata$x <- xy.id$x - xoffset
    pred0.raw <- predict(object.nlme, newdata, level=0L) - xy.id$y
    pred0 <- yfun(pred0.raw)
    if (deriv > 0L) {
# VELOCITY
# level 0 prediction
      if (!'.groups' %in% names(newdata)) {
        pred0 <- get_vel(xfun(x), pred0, deriv)$y
# or if grouped opt dv
      } else {
        pred0 <- newdata %>%
          mutate(x = !!x,
                 pred0 = !!pred0,
                 n = 1:n()) %>%
          nest_by(.data$.groups) %>%
          mutate(vel = list(with(.data$data, get_vel(xfun(x), pred0, deriv)$y))) %>%
          unnest(cols = c(.data$data, .data$vel)) %>%
          arrange(n) %>%
          pull(.data$vel)
      }
      if (any(level == 1L)) {
# level 1 prediction
# velocity curve on back-transformed axes
    # shift x to mean curve equivalents
        x.id <- xyadj(object, x = x, id = id, abc = abc)$x
    # unique mean spline curve on transformed x-y scales
        pred <- get_vel(x, pred0.raw, deriv) %>% as_tibble() %>% unique() %>%
    # fit mean velocity curve
          spline(., method = 'natural', xout = x.id) %>% as_tibble() %>% pull(.data$y) %>%
    # shift and scale to individual velocities
          xyadj(object, x = 0, v = ., id = id, abc = abc, tomean = FALSE) %>% as_tibble() %>%
    # adjust for y and x transformations
          mutate(v = .data$v / Dxy(object, pred, 'y') * Dxy(object, xfun(!!x), 'x')) %>% pull(.data$v)
      }
    }
# return data frame if level 0:1
    if (length(level) > 1L)
      return(data.frame(id=factor(id), predict.fixed=pred0, predict.id=pred))
# add names or split by id if level 1
    if (level == 0L)
      pred <- pred0
    asList <- ifelse(is.null(asList <- eval(mc$asList)), FALSE, asList)
    if (asList)
      pred <- split(pred, id)
    else
      names(pred) <- id
    attr(pred, 'label') <- 'Predicted values'
    return(pred)
  }
