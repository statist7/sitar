  predict.sitar <- function(object, newdata, level=1, ...) {
    mcall <- match.call()[-1]
    if (!missing(newdata)) {
      scall <- object$call.sitar
      if (!'x' %in% names(newdata)) stop('newdata lacks x variable')
      if (!'id' %in% names(newdata)) {
        if (level != 0) stop('newdata lacks id variable')
        newdata <- data.frame(newdata, id=object$groups$id[1])
      }
      fe <- fixef(object)
# omit fixed effects already in newdata
      fe <- fe[!match(names(fe), names(newdata), 0)]
      newdata <- data.frame(newdata, t(fe))
      mcall$newdata <- quote(newdata)
# attach object for fitnlme
      attach(object)
      on.exit(detach(object))
    }
    pred <- tryCatch(do.call(nlme:::predict.nlme, as.list(mcall)), error=function(e) {
        print(e)
        if ('could not find function "fitnlme"' %in% e)
          return(paste('need to refit model', deparse(mcall[[1]])))
    })
    attributes(pred) <- NULL
    pred
  }

  getData.sitar <- function(object) {
    object$call <- object$call.sitar
    nlme:::getData.nlme(object)
  }

  getVarCov.sitar <- function(obj, ...) {
    class(obj) <- 'lme'
    getVarCov(obj)
  }
