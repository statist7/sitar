#' Adjust x and y variables for SITAR random effects
#'
#' \code{xyadj} Adjusts \code{x} and \code{y} and optionally \code{v} values for subject-specific
#' random effects from a SITAR model.
#'
#' When \code{tomean = TRUE} the x and y and v values are adjusted to
#' \deqn{(x - xoffset - b<fixed> - b<random>) * exp(c<random>) + xoffset + b<fixed>}
#' \deqn{y - a<random> - d<random> * x}
#' \deqn{(v - d<random>) / exp(c<random>)}
#' When \code{tomean = FALSE} they are adjusted to
#' \deqn{(x - xoffset - b<fixed>) / exp(c<random>) + xoffset + b<fixed> + b<random>}
#' \deqn{y + a<random> + d<random> * x}
#' \deqn{v * exp(c<random>) + d<random>}
#' In each case missing values of the fixed or random effects are
#' set to zero.
#'
#' @param object a SITAR model.
#' @param x a vector of x coordinates. If missing, \code{x} and
#' \code{y} and \code{id} are obtained from \code{object}.
#' @param y a vector of y coordinates (default 0).
#' @param v a vector of velocity coordinates (default 0).
#' @param id a factor denoting the subject levels corresponding to \code{x} and
#' \code{y} and \code{v}.
#' @param abc a data frame containing random effects for a, b, c and d (default
#' \code{ranef(object)[id, ]}).
#' @param tomean a logical defining the direction of adjustment. TRUE (default)
#' indicates that individual curves are translated and rotated to match the
#' mean curve, while FALSE indicates the reverse, the mean curve being
#' translated and rotated to match individual curves.
#' @return The list of adjusted values: \item{x}{numeric vector.}
#' \item{y}{numeric vector the same length as x, or NULL.}
#' \item{v}{numeric vector the same length as x, or NULL.}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#'
#' data(heights)
#' ## fit sitar model for height
#' m1 <- sitar(x = age, y = height, id = id, data = heights, df = 5)
#'
#' ## plot unadjusted data as growth curves
#' plot(m1, opt='u')
#'
#' ## overplot with adjusted data as points
#' with(heights, points(xyadj(m1), col='red', pch = 19))
#'
#' @export
xyadj <- function(object, x, y = 0, v = 0, id, abc = NULL, tomean = TRUE) {
#	returns x and y adjusted for random effects a, b and c
  if (missing(x)) {
    x <- getCovariate(object)
    y <- getResponse(object)
    id <- getGroups(object)
  }
  # add missing columns
  if (is.null(abc)) {
    re <- ranef(object)
    abc <- re[match(id, rownames(re)), , drop = FALSE]
  }
  abc <- as.data.frame(abc)
  for (i in letters[1:4])
    if (!i %in% names(abc))
      abc[, i] <- 0
  xoffset <- object$xoffset
  if (!is.na(b0 <- fixef(object)['b']))
    xoffset <- xoffset + b0
  x <- x - xoffset
  if (tomean) {
    x.adj <- (x - abc$b) * exp(abc$c) + xoffset
    y.adj <- y - abc$a - abc$d * x
    v.adj <- (v - abc$d) / exp(abc$c)
  } else {
    x.adj <- x / exp(abc$c) + abc$b + xoffset
    y.adj <- y + abc$a + abc$d * x
    v.adj <- v * exp(abc$c) + abc$d
  }
  return(list(x = x.adj, y = y.adj, v = v.adj))
}
