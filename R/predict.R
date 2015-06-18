  predict.sitar <- function(object, newdata, level=1, ..., deriv=0, xfun=I, yfun=I) {
# create x and id variables in newdata
    if (missing(newdata)) newdata <- getData(object)
    oc <- object$call.sitar
    if (!is.null(newdata$x)) x <- newdata$x else
      newdata$x <- x <- eval(oc$x, newdata)
    if (!is.null(newdata$id)) id <- newdata$id else
      newdata$id <- id <- if (any(level == 1)) eval(oc$id, newdata) else
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
    newdata <- data.frame(newdata, t(args))
    class(object) <- class(object)[-1]
    pred <- predict(object=object, newdata=newdata, level=level, ...)
    if (deriv == 0) {
      pred <- yfun(pred)
      attributes(pred) <- NULL
      return(pred)
    }
# deriv = 1
    if (level == 1) pred <- predict(object=object, newdata=newdata, level=0, ...)
# level 0 velocity predictions
    vel0 <- predict(makess(x, pred, xfun=xfun, yfun=yfun), xfun(x), deriv=1)
    if (level == 0) return(vel0$y)
# level 1 velocity predictions
    vel1 <- spline(vel0, method='natural',
                   xout=xfun(xyadj(x=x, id=id, object=object)$x))
    if (!is.null(c <- ranef(object)$c)) vel1$y <- vel1$y * exp(c[id])
    return(vel1$y)
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
