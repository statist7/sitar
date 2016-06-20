#' Compare Likelihoods of Fitted SITAR Objects
#'
#' anova method for \code{sitar} objects, based on \code{anova.lme}.
#'
#'
#' @param object an object inheriting from class \code{sitar}.
#' @param \dots other optional fitted model objects.
#' @param test an optional logical value controlling whether likelihood ratio
#' tests should be used.
#' @param type an optional character string specifying the type of sum of
#' squares to be used.
#' @param adjustSigma see \code{\link{anova.lme}}.
#' @param Terms see \code{\link{anova.lme}}.
#' @param L see \code{\link{anova.lme}}.
#' @param verbose an optional logical value.
#' @return a data frame inheriting from class "anova.lme".
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @export
anova.sitar <- function (object, ..., test = TRUE, type = c("sequential", "marginal"),
                         adjustSigma = TRUE, Terms, L, verbose = FALSE) {
  name <- deparse(substitute(object))
  class(object) <- class(object)[class(object) != 'sitar']
  assign(name, object)
  dots <- match.call(expand.dots=FALSE)$...
  if (length(dots) > 0) for (obj in dots) {
    name <- deparse(obj)
    obj <- eval(obj)
    class(obj) <- class(obj)[class(obj) != 'sitar']
    assign(name, obj)
  }
  do.call('anova', as.list(match.call()[-1]))
}
