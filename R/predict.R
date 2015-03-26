  predict.sitar <- function(object, newdata, level, ...) {
    mcall <- match.call()[-1]
    if (!missing(newdata)) {
      on.exit(if (exists('.fitnlme')) rm(.fitnlme, envir=globalenv()))
      scall <- object$call.sitar
      if (!'x' %in% names(newdata)) stop('newdata lacks x variable')
      if (!'id' %in% names(newdata)) {
        if (level != 0) stop('newdata lacks id variable')
        newdata <- data.frame(newdata, id=object$groups$id[1])
      }
      fe <- fixef(object)
      fe <- fe[!match(names(fe), names(newdata), 0)]
      newdata <- data.frame(newdata, t(fe))
      mcall$newdata <- quote(newdata)
      if (exists('object$.fitnlme')) .fitnlme <<- object$.fitnlme else
        .fitnlme <<- update(object, returnsub=TRUE)
    }
    do.call(nlme:::predict.nlme, list(object, newdata, level, ...))
  }

  getData.sitar <- function(object) {
    object$call <- object$call.sitar
    nlme:::getData.nlme(object)
  }

  getVarCov.sitar <- function(obj, ...) {
    class(obj) <- 'lme'
    getVarCov(obj)
  }
