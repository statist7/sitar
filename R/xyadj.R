xyadj <- function(object, data) {

#	returns x and y adjusted for random effects a, b and c
  mcall <- object$call.sitar
  if (missing(data)) {
    data <- eval(mcall$data)
    subset <- eval(mcall$subset, data)
    if (!is.null(subset)) data <- data[subset,]
  }
  xoffset <- object$xoffset
  if (is.null(xoffset)) xoffset <- 0
  if (!is.na(fixef(object)['b'])) xoffset <- xoffset + fixef(object)['b']
  x.adj <- eval(mcall$x, data) - xoffset
  id <- factor(eval(mcall$id, data))
  re <- ranef(object)
  if (!is.null(re$b)) x.adj <- x.adj - re$b[id]
  if (!is.null(re$c)) x.adj <- x.adj * exp(re$c[id])
  x.adj <- x.adj + xoffset
  y.adj <- try(eval(mcall$y, data), silent=TRUE)
  if (class(y.adj) == 'try-error') y.adj <- NULL else
    if (!is.null(re$a)) y.adj <- y.adj - re$a[id]
  return(list(x=x.adj, y=y.adj))
}
