#' Adjust x and y variables for SITAR random effects
#'
#' \code{xyadj} Adjusts \code{x} and \code{y} values for subject-specific
#' random effects from a SITAR model.
#'
#' When tomean = TRUE the x and y values are adjusted to \deqn{(x - xoffset -
#' b<fixed> - b<random>) * exp(c<random>) + xoffset + b<fixed>} \deqn{y -
#' a<random>} When tomean = FALSE they are adjusted to \deqn{(x - xoffset -
#' b<fixed>) / exp(c<random>) + xoffset + b<fixed> + b<random>} \deqn{y +
#' a<random>} In each case missing values of the fixed or random effects are
#' set to zero.
#'
#' @param object a SITAR model.
#' @param x a vector of x coordinates. If missing, \code{x} and
#' \code{y} and \code{id} are obtained from \code{object}.
#' @param y a vector of y coordinates (default NULL).
#' @param id a factor denoting the subject levels corresponding to \code{x} and
#' \code{y}.
#' @param abc a data frame containing random effects for a, b and c (default
#' \code{ranef(object)[id, ]}).
#' @param tomean a logical defining the direction of adjustment. TRUE (default)
#' indicates that individual curves are translated and rotated to match the
#' mean curve, while FALSE indicates the reverse, the mean curve being
#' translated and rotated to match individual curves.
#' @return The list of adjusted values: \item{x}{numeric vector.}
#' \item{y}{numeric vector the same length as x, or NULL.}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#'
#' data(heights)
#' ## fit sitar model for height
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)
#'
#' ## plot unadjusted data as growth curves
#' plot(m1, opt='u')
#'
#' ## overplot with adjusted data as points
#' with(heights, points(xyadj(m1), col='red', pch=19))
#'
#' @export xyadj
xyadj <- function(object, x, y=NULL, id, abc=ranef(object)[id, , drop=FALSE], tomean=TRUE) {
#	returns x and y adjusted for random effects a, b and c
  if (missing(x)) {
    x <- getCovariate(object)
    y <- getResponse(object)
    id <- getGroups(object)
  }
  abc[, letters[1:3][!letters[1:3] %in% names(ranef(object))]] <- 0 # omit not in model
  xoffset <- object$xoffset
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
