  predict.sitar <- function(object, newdata, level=1, ...) {
    mcall <- match.call()[-1]
    if (!missing(newdata)) {
      if (!'x' %in% names(newdata)) stop('newdata lacks x variable')
      if (!'id' %in% names(newdata)) {
        if (level != 0) stop('newdata lacks id variable')
        newdata <- data.frame(newdata, id=object$groups$id[1])
      }
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
    }
    class(object) <- class(object)[class(object) != 'sitar']
    pred <- predict(object=object, newdata=newdata, level=level, ...)
    attributes(pred) <- NULL
    pred
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
