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
#' predictions. It requires named columns for \code{x}, and for \code{id} if
#' \code{level = 1}, matching the names in \code{object}. Any covariates in
#' \code{a.formula}, \code{b.formula} or \code{c.formula} can also be included.
#' By default their values are set to the mean, so when \code{level = 0} the
#' prediction represents the mean curve.
#' @param level an optional integer giving the level(s) of grouping to be used
#' in obtaining the predictions, level 0 corresponding to the population
#' predictions. Defaults to level 1.
#' @param \dots other optional arguments: \code{asList}, \code{na.action} and
#' \code{naPattern}.
#' @param deriv an optional integer specifying predictions corresponding to
#' either the fitted curve or its derivative. \code{deriv = 0} (default)
#' specifies the distance curve, \code{deriv = 1} the velocity curve and
#' \code{deriv = 2} the acceleration curve.
#' @param abc an optional named vector containing values of a subset of
#' \code{a}, \code{b} and \code{c}, default \code{NULL}. If \code{abc} is set,
#' \code{level} is set to 0. It gives predictions for a single subject with the
#' specified values of \code{a}, \code{b} and \code{c}, where missing values
#' are set to 0.
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
  predict.sitar <- function(object, newdata=getData(object), level=1, ...,
                            deriv=0, abc=NULL,
                            xfun=function(x) x, yfun=function(y) y) {
# create x and id variables in newdata
    oc <- object$call.sitar
    if (is.null(newdata$x)) newdata$x <- eval(oc$x, newdata)
    x <- newdata$x
    if (is.null(xoffset <- object$xoffset)) {
      xoffset <- mean(getCovariate(object))
      warning('xoffset set to mean(x) - best to refit model')
    }
    newdata$x <- newdata$x - xoffset
# create id in newdata
    if (is.null(newdata$id)) {
      if (any(level == 1)) newdata$id <- eval(oc$id, newdata)
      else newdata$id <- rep.int(getGroups(object)[1], nrow(newdata))
    }
    id <- newdata$id <- factor(newdata$id)
# check abc as length-3 vector or 1-row data frame
    if (is.null(abc)) abc <- ranef(object)
    else if (is.vector(abc)) {
#	abc is named vector
      if (!is.null(names(abc)) && all(unique(names(abc)) %in% letters[1:3])) {
        abc <- data.frame(t(abc))
        abc[, letters[1:3][!letters[1:3] %in% names(abc)]] <- 0 # fill with zeros
      } else
#	abc is id level
      if (length(abc) == 1) {
        . <- rownames(ranef(object)) %in% abc
        if (!any(.)) stop(paste('id', abc, 'not found'))
        abc <- ranef(object)[., , drop=FALSE]
      }
    }
    if (is.null(nrow(abc)))
      stop('abc should be either a single id level or up to three named random effect values')
    else if (abcset <- nrow(abc) == 1) {
      level <- 0
    } else
      abc <- abc[id, , drop=FALSE]
# check if old-style object lacking fitnlme
    if(!'fitnlme' %in% names(object)) {
      warning('fitnlme missing - best to refit model')
      object <- update(object, control=nlmeControl(maxIter=0, pnlsMaxIter=0, msMaxIter=0))
    }
# attach object for fitnlme
    on.exit(detach(object))
    eval(parse(text='attach(object)'))
# identify covariates needed in newdata, omitting x, fixed effects and random effects
    argnames <- names(formals(fitnlme))
    argnames <- argnames[!argnames %in% names(fixef(object))][-1]
    argnames <- argnames[!argnames %in% names(ranef(object))]
    if (length(argnames) > 0) {
# check if newdata subsetted (from plot)
      if (is.null(subset <- attr(newdata, 'subset'))) {
# identify covariates in newdata other than x and id
    		covnames <- names(newdata)
    		covnames <- covnames[!covnames %in% c('x', 'id')]
# set to 0 covariates not in newdata
    		newdata[, argnames[!argnames %in% covnames]] <- 0
# centre each needed covariate in newdata
    		covnames <- covnames[covnames %in% argnames]
    		if (length(covnames) > 0) {
    		  gd <- getData(object)
    		  for (i in covnames) {
# continuous variable
    		    if (i %in% argnames) newdata[[i]] <- newdata[[i]] - mean(gd[[i]])
    	      else {
# factor as instrumental variable(s)
    	        lev <- levels(gd[[i]])
    	        for (j in 2:length(lev)) {
    	          k <- paste0(i, lev[j])
    	          newdata[[k]] <- as.numeric(newdata[[i]] == lev[j]) - mean(gd[[i]] == lev[j])
    	        }
    	      }
    		  }
    		}
      }
# newdata subsetted (in plot)
      else {
        gd <- update(object, returndata=TRUE)[subset, argnames, drop=FALSE]
        argnames <- unlist(lapply(gd, mean))
        newdata <- data.frame(newdata, t(argnames))
      }
    }
# set class to nlme
    class(object) <- class(object)[-1]
# ensure deriv integral
    deriv <- as.integer(deriv)
# simple prediction
    if (deriv == 0 && !abcset) {
      pred <- yfun(predict(object=object, newdata=newdata, level=level, ...))
      return(pred)
    }
# complex prediction
    else { # deriv != 0 || abcset
# mean distance curve
      pred <- predict(object=object, newdata=newdata, level=0, ...)
# x changed to reflect individual b and c
      xout <- xyadj(object=object, x=x, id=id, abc=abc)$x
# DISTANCE
      if (deriv == 0) { # abcset
# level 1 prediction
        pred <- spline(list(x=x, y=pred), method='natural', xout=xout)$y
# add individual a to prediction (inexact when yfun != y)
        if (!is.null(abc$a)) pred <- yfun(pred + abc$a)
      }
# VELOCITY
      else { # deriv != 0
# mean velocity curve on back-transformed axes
        vel0 <- predict(smooth.spline(xfun(x), yfun(pred)), xfun(x), deriv=deriv)
        if (any(level == 0) && !abcset) pred0 <- pred <- vel0$y
        if (any(level == 1) || abcset) {
# level 1 prediction
          pred <- spline(vel0, method='natural', xout=xfun(xout))$y
# multiply velocity by individual c (inexact when xfun != x)
          if (!is.null(abc$c)) pred <- pred * exp(abc$c)
        }
      }
# return data frame if level 0:1
      if (length(level) > 1) return(data.frame(id=id, predict.fixed=pred0, predict.id=pred))
# add names or split by id if level 1
      if (level == 1) {
        asList <- ifelse(is.null(asList <- list(...)$asList), FALSE, asList)
        if (asList) pred <- split(pred, id) else names(pred) <- id
      }
      attr(pred, 'label') <- 'Predicted values'
      return(pred)
    }
  }
