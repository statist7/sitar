  predict.sitar <- function(object, newdata, level=1, ..., deriv=0,
                            abc=NULL, xfun=I, yfun=I) {
# function expand
    expand <- function(z) if (is.null(z)) 0 else z[id]
# create x and id variables in newdata
    if (missing(newdata)) newdata <- getData(object)
    oc <- object$call.sitar
    if (is.null(newdata$x)) newdata$x <- eval(oc$x, newdata)
    x <- newdata$x
    xoffset <- object$xoffset
    if (is.null(xoffset)) {
      xoffset <- mean(getCovariate(object))
      warning('xoffset set to mean(x) - best to refit model')
    }
    newdata$x <- newdata$x - xoffset
# check abc
    if (!is.null(abc)) {
      abcset <- TRUE
      if (!'data.frame' %in% class(abc)) abc <- data.frame(t(abc))
      if (nrow(abc) == 1) level <- 0
        else if (nrow(abc) == nrow(ranef(object))) abcset <- FALSE
                                   else stop('abc not a variate of length up to 3')
    }
    else {
      abcset <- FALSE
      abc <- ranef(object)
    }
# check for id in newdata
    if (id.ok <- deparse(oc$id) %in% names(newdata)) newdata$id <- factor(eval(oc$id, newdata))
    else newdata$id <- factor(rep.int(getGroups(object)[1], nrow(newdata)))
    id <- factor(newdata$id)
# check if old-style object lacking fitnlme
    if(!'fitnlme' %in% names(object)) {
      cat('need to update object to obtain fitnlme\n')
      object <- update(object, control=nlmeControl(maxIter=0, pnlsMaxIter=0, msMaxIter=0))
    }
# attach object for fitnlme
    on.exit(detach(object))
    eval(parse(text='attach(object)'))
# identify extra covariates needed in newdata, omitting x and fixed effects
		argnames <- names(formals(fitnlme))
		argnames <- argnames[!match(argnames, names(fixef(object)), 0)][-1]
# zero or centre each covariate
		if (length(argnames) > 0) {
		  if (any(match(argnames, names(newdata), 0))) gd <- getData(object)
		  for (i in argnames) {
		    if (match(i, names(newdata), 0)) newdata[[i]] <- newdata[[i]] - mean(gd[[i]])
		    else newdata[[i]] <- 0
		  }
		}
# set class to nlme
    class(object) <- class(object)[-1]
# simple prediction
    if (deriv == 0 && !abcset) {
      pred <- yfun(predict(object=object, newdata=newdata, level=level, ...))
    }
# complex prediction
    else {
# mean distance curve
      pred <- predict(object=object, newdata=newdata, level=0, ...)
# DISTANCE
      if (deriv == 0) {
# mean distance curve back-transformed
        if (level == 0 && !abcset) pred <- yfun(pred)
# level 1 prediction based on x changed to reflect individual b and c
        else {
          pred <- spline(list(x=x, y=pred), method='natural',
                         xout=xyadj(x=x, id=id, object=object, abc=abc)$x)$y
# add individual a to prediction
          if (!is.null(ranef(object)$a)) pred <- yfun(pred + expand(abc$a))
        }
      }
# VELOCITY
      else {
# mean velocity curve on back-transformed axes
        vel0 <- predict(makess(x, pred, xfun=xfun, yfun=yfun), xfun(x), deriv=1)
        if (level == 0 && !abcset) pred <- vel0$y
# level 1 prediction based on x changed to reflect individual b and c
        else {
          pred <- spline(vel0, method='natural',
                         xout=xfun(xyadj(x=x, id=id, object=object, abc=abc)$x))$y
# multiply velocity by individual c (inexact when xfun != I)
          if (!is.null(ranef(object)$c)) pred <- pred * exp(expand(abc$c))
        }
      }
    }
    if (id.ok && is.null(dim(pred))) names(pred) <- id
    return(pred)
  }

  getData.sitar <- function(object) {
    object$call <- object$call.sitar
    class(object) <- 'lme'
    getData(object)
  }

  getVarCov.sitar <- function(obj, ...) {
    class(obj) <- 'lme'
    getVarCov(obj)
  }

  getCovariate.sitar <- function (object, ...)
  {
    eval(object$call.sitar$x, getData(object))
  }
