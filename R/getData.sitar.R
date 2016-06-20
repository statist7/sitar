#' Extract elements of fitted SITAR models
#'
#' getData, getCovariate and getVarCov methods for \code{sitar} objects,
#' based on \code{lme}.
#'
#'
#' @param object,obj an object inheriting from class \code{sitar}.
#' @param \dots other optional arguments.
#' @return Respectively the data frame and \code{x} variable
#' used in the fit, and the returned variance-covariance matrix.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @importFrom nlme getData
#' @method getData sitar
#' @export
  getData.sitar <- function(object) {
    object$call <- object$call.sitar
    class(object) <- 'lme'
    getData(object)
  }


#' @rdname getData.sitar
#' @importFrom nlme getCovariate
#' @method getCovariate sitar
#' @export
  getCovariate.sitar <- function(object, ...)
  {
    eval(object$call.sitar$x, getData(object))
  }


#' @rdname getData.sitar
#' @importFrom nlme getVarCov
#' @method getVarCov sitar
#' @export
  getVarCov.sitar <- function(obj, ...) {
    class(obj) <- 'lme'
    getVarCov(obj)
  }
