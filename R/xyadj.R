xyadj <- function(x, y=NULL, id, object, abc=NULL, tomean=TRUE) {
#	returns x and y adjusted for random effects a, b and c
  re <- ranef(object)
  if (is.null(abc)) abc <- re
  abc[, letters[1:3][!letters[1:3] %in% names(abc)]] <- 0 # fill with zeros
  abc[, letters[1:3][!letters[1:3] %in% names(re)]] <- 0 # omit not in model
  abc <- abc[id, ]
  if (is.null(xoffset <- object$xoffset)) xoffset <- mean(getCovariate(object))
  if (!is.na(b0 <- fixef(object)['b'])) xoffset <- xoffset + b0
  if (tomean) {
    x.adj <- (x - xoffset - abc$b) * exp(abc$c) + xoffset
    y.adj <- y - abc$a
  } else {
    x.adj <- (x - xoffset) / exp(abc$c) + xoffset + abc$b
    y.adj <- y + abc$a
  }
  return(list(x=x.adj, y=y.adj))
}
