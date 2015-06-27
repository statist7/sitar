  predict.sitar <- function(object, newdata, level=1, ..., deriv=0,
                            abc=ranef(object), xfun=I, yfun=I) {
# function expand
    expand <- function(z) if (is.null(z)) 0 else if (nrow(abc) == 1) z else z[id]
# create x and id variables in newdata
    if (missing(newdata)) newdata <- getData(object)
    oc <- object$call.sitar
    if (!is.null(newdata$x)) x <- newdata$x else
      newdata$x <- x <- eval(oc$x, newdata)
    if (!is.null(newdata$id)) id <- factor(newdata$id) else
      newdata$id <- id <- if (any(level == 1)) factor(eval(oc$id, newdata)) else
        rep.int(getGroups(object)[1], nrow(newdata))
# attach object for fitnlme
    on.exit(detach(object))
    eval(parse(text='attach(object)'))
    if (!exists('fitnlme'))
      stop('could not find function "fitnlme": please update model')
# omit fitnlme args already in newdata
		argnames <- names(formals(fitnlme))
		args <- setNames(vector('integer', length=length(argnames)), argnames)
    args <- args[!match(argnames, names(newdata), 0)]
# add other fitnlme args to newdata
    newdata <- data.frame(newdata, t(args))
# set class to nlme
    class(object) <- class(object)[-1]
# simple prediction
    if (deriv == 0 && missing(abc)) {
      pred <- yfun(predict(object=object, newdata=newdata, level=level, ...))
    }
# complex prediction
    else {
# mean distance curve
      pred <- predict(object=object, newdata=newdata, level=0, ...)
# DISTANCE
      if (deriv == 0) {
# mean distance curve back-transformed
        if (level == 0) pred <- yfun(pred)
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
        if (level == 0) pred <- vel0$y
# level 1 prediction based on x changed to reflect individual b and c
        else {
          pred <- spline(vel0, method='natural',
                         xout=xfun(xyadj(x=x, id=id, object=object, abc=abc)$x))$y
# multiply velocity by individual c (inexact when xfun != I)
          if (!is.null(ranef(object)$c)) pred <- pred * exp(expand(abc$c))
        }
      }
    }
    attributes(pred) <- NULL
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
