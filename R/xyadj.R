xyadj <- function(object) {
  #	returns x and y adjusted for random effects a, b and c
  random <- as.character(object$call$random)[[2]]
  mcall <- object$call.sitar
  data <- eval(mcall$data)
  subset <- eval(mcall$subset, data)
  if (!is.null(subset)) data <- data[subset,]
  x <- eval(mcall$x, data)
  y <- eval(mcall$y, data)
  id <- eval(mcall$id, data)
  nf <- length(fitted(object))
  if (nf != length(y)) stop(paste('model (length=', nf, ') incompatible with data (rows=', length(y), ')', sep=''))
  xoffset <- object$xoffset
  if (is.null(xoffset)) xoffset <- 0
  if (!is.na(fixef(object)['b'])) xoffset <- xoffset + fixef(object)['b']
  x.adj <- x - xoffset
  if (grepl('b', random)) x.adj <- x.adj - ranef(object)[factor(id),'b']
  if (grepl('c', random)) x.adj <- x.adj * exp(ranef(object)[factor(id),'c'])
  x.adj <- x.adj + xoffset
  y.adj <- y
  if (grepl('a', random)) y.adj <- y.adj - ranef(object)[factor(id),'a']
  invisible(data.frame(x=x.adj, y=y.adj))
}
