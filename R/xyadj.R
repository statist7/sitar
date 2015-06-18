xyadj <- function(x, y=NULL, id, object, abc=ranef(object), tomean=TRUE) {
#	returns x and y adjusted for random effects a, b and c
  if (is.null(xoffset <- object$xoffset)) xoffset <- 0
  if (!is.na(b0 <- fixef(object)['b'])) xoffset <- xoffset + b0
  if (nrow(abc) > 1) abc <- abc[id, ]
  abc[, letters[1:3][!letters[1:3] %in% names(abc)]] <- 0
  if (tomean) {
    x.adj <- (x - xoffset - abc$b) * exp(abc$c) + xoffset
    y.adj <- y - abc$a
  } else {
    x.adj <- (x - xoffset) / exp(abc$c) + xoffset + abc$b
    y.adj <- y + abc$a
  }
  return(list(x=x.adj, y=y.adj))
}
