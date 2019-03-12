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
#' \code{a.formula}, \code{b.formula} or \code{c.formula} can also be included.
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
#' \code{a}, \code{b} and \code{c}, default \code{NULL}. Ignored if
#' \code{level = 0}. It gives predictions for a single subject with the
#' specified values of \code{a}, \code{b} and \code{c}, where missing values
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
#' @export
  predict.sitar <- function(object, newdata=getData(object), level=1L, ...,
                            deriv=0L, abc=NULL,
                            xfun=function(x) x, yfun=function(y) y) {
# create x in newdata
    oc <- object$call.sitar
    x <- if ('.x' %in% names(newdata))
      newdata$.x
    else
      eval(oc$x, newdata)
    if (is.null(xoffset <- object$xoffset)) {
      xoffset <- mean(getCovariate(object))
      warning('xoffset set to mean(x) - best to refit model')
    }
    newdata$x <- x - xoffset
# check abc as id
    id <- rownames(ranef(object))
    if (all(level == 0L))
      abc <- NULL
    else
    if (!is.null(abc)) {
      if (is.null(names(abc))) {
        if (is.vector(abc) && length(abc) == 1 && as.character(abc) %in% id) {
          newdata$id <- abc
          newdata$.id <- abc <- NULL
        }
        else
          stop('abc unrecognised as id')
      } else {
# check abc as length 3 vector or data frame
        if (is.vector(abc))
          abc <- data.frame(t(abc))
        if (!all(names(abc) %in% letters[1:3]))
          stop('abc unrecognised as data frame')
        abc[, letters[1:3][!letters[1:3] %in% names(abc)]] <- 0 # fill with zeros
      }
    }
# create id in newdata
    newdata$id <- if ('.id' %in% names(newdata))
      newdata$.id
    else {
      if (any(level == 1L) && (is.null(abc) || nrow(abc) > 1))
        eval(oc$id, newdata)
      else
        factor(1, labels=id[1])
    }
    id <- newdata$id
    if (!is.null(abc))
      abc <- abc[id, ] # expand with id
# check if old-style object lacking fitnlme
    if (!'fitnlme' %in% names(object)) {
      warning('fitnlme missing - best to refit model')
      object <- update(object, control=nlmeControl(maxIter=0, pnlsMaxIter=0, msMaxIter=0))
    }
# attach object for fitnlme
    on.exit(detach(object))
    eval(parse(text='attach(object)'))
# identify covariates in model (not x or coef)
    argnames <- names(formals(fitnlme))
    argnames <- argnames[!argnames %in% names(coef(object))][-1]
    if (length(argnames) > 0) {
# identify model covariates in newdata
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
# check if newdata subsetted (from plot)
    subset <- attr(newdata, 'subset')
    if (!is.null(subset)) {
      if (!is.null(abc))
        stop('use subset or abc but not both')
# create abc for subset
      re <- ranef(object)
      abc <- apply(re[rownames(re) %in% getGroups(object)[subset], , drop=FALSE], 2, mean)
      abc <- data.frame(t(abc))
      level <- 1
# set to mean gd, covariates not in newdata
      if (exists('notnames') && length(notnames) > 0) {
        if (!exists('gd'))
          gd <- update(object, returndata=TRUE)
        for (i in notnames)
          newdata[, i] <- mean(gd[subset, i])
      }
    }
# set class to nlme
    class(object) <- class(object)[-1]
# ensure deriv integral
    deriv <- as.integer(deriv)
# simple prediction
    if (all(deriv == 0L) && is.null(abc))
      return(yfun(predict(object, newdata, level=level, ...)))
# DISTANCE
# level 0 prediction
    pred0 <- yfun(predict(object, newdata, level=0L))
    xy.id <- xyadj(object, x=x, y=0, id=id, abc=abc)
# level 1 prediction
    newdata$x <- xy.id$x - xoffset
    pred <- yfun(predict(object, newdata, level=0L) - xy.id$y)
    if (any(deriv > 0L)) {
# VELOCITY
# level 0 prediction
      ss0 <- smooth.spline(xfun(x), pred0)
      vel0 <- predict(ss0, xfun(x), deriv=max(deriv))
      pred0 <- vel0$y
# velocity curve on back-transformed axes
      if (any(level == 1L)) {
# level 1 prediction
        if (body(xfun) == as.name('x') &&
            body(yfun) == as.name('x')) {
# x and y untransformed
            vel <- spline(vel0, method='natural', xout=xfun(xy.id$x))$y
            if (is.null(abc))
              abc <- ranef(object)[id, , drop=FALSE]
            if (!is.null(abc$c))
              vel <- vel * exp(abc$c)
          } else {
# x or y transformed
          newdata$pred <- pred
          newdata$xorig <- xfun(x)
          vel <- by(newdata, id, function(z) {
            with(z, {
              if (length(xorig) >= 4) {
                ss <- smooth.spline(xorig, pred, df=min(20, length(xorig)))
                predict(ss, xorig, deriv=max(deriv))$y
              } else
                predict(ss0, xorig, deriv=max(deriv))$y
            })
          })
          vel <- do.call('c', vel)
        }
        pred <- vel
      }
    }
# return data frame if level 0:1
    if (length(level) > 1)
      return(data.frame(id=factor(id), predict.fixed=pred0, predict.id=pred))
# add names or split by id if level 1
    if (level == 0L)
      pred <- pred0
    asList <- ifelse(is.null(asList <- list(...)$asList), FALSE, asList)
    if (asList)
      pred <- split(pred, id)
    else
      names(pred) <- id
    attr(pred, 'label') <- 'Predicted values'
    return(pred)
  }
