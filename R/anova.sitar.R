anova.sitar <- function (object, ..., test = TRUE, type = c("sequential", "marginal"),
                         adjustSigma = TRUE, Terms, L, verbose = FALSE) {
  mcall <- match.call()
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
  mcall[[1]] <- as.name("anova")
  eval(mcall)
}
