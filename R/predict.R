  predict.sitar <- function(object, newdata, level, ...) {
    mcall <- match.call()[-1]
    if (!missing(newdata)) {
      on.exit(rm(.fitnlme, envir=globalenv()))
      scall <- object$call.sitar
      xname <- deparse(scall$x)
      if (!xname %in% names(newdata)) stop('newdata lacks x variable')
      names(newdata)[match(xname, names(newdata), 0)] <- 'x'
      idname <- deparse(scall$id)
      if (!idname %in% names(newdata)) {
        if (level != 0) stop('newdata lacks id variable')
        newdata <- data.frame(newdata, id=object$groups$id[1])
        names(newdata)[match('id', names(newdata), 0)] <- idname
      }
      fe <- fixef(object)
      fe <- fe[!match(names(fe), names(newdata), 0)]
      newdata <- data.frame(newdata, t(fe))
      mcall$newdata <- quote(newdata)
      .fitnlme <<- object$.fitnlme
    }
    do.call(nlme:::predict.nlme, as.list(mcall))
  }

  getData.sitar <- function(object) {
    object$call <- object$call.sitar
    nlme:::getData.nlme(object)
  }

  getVarCov.sitar <- function(obj, ...) {
    class(obj) <- 'lme'
    getVarCov(obj)
  }
