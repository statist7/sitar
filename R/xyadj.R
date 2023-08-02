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
#' @importFrom rlang .data %||%
#' @export
xyadj <- function(object, x, y = 0, v = 0, id, abc = NULL, tomean = TRUE) {
#	returns x, y and v adjusted for random effects a, b, c and d
  if (missing(x)) {
    x <- getCovariate(object)
    y <- getResponse(object)
    id <- getGroups(object)
  }
  xoffset <- object$xoffset
  if (!is.na(b0 <- fixef(object)['b']))
    xoffset <- xoffset + b0
  # set up abc
  if (is.null(abc)) {
    re <- ranef(object)
    abc <- re[match(id, rownames(re)), , drop = FALSE]
  }
  abc <- as.data.frame(abc)
  # add missing a:d columns
  abc[, setdiff(letters[1:4], names(abc))] <- 0 # chatGPT suggestion
  abc <- abc %>%
    mutate(x = x - xoffset,
           d.adjusted = attr(object, 'd.adjusted') %||% FALSE)
  if (tomean) {
    abc <- abc %>%
      mutate(x.adj = (x - .data$b) * exp(.data$c) + xoffset,
             y.adj = y - .data$a - .data$d * if_else(.data$d.adjusted,
                                                     .data$x.adj - xoffset,
                                                     x),
             v.adj = if_else(.data$d.adjusted,
                             v / exp(.data$c) - .data$d,
                             (v - .data$d) / exp(.data$c)))
  } else {
    abc <- abc %>%
      mutate(x.adj = x / exp(.data$c) + .data$b + xoffset,
             y.adj = y + .data$a + .data$d * if_else(.data$d.adjusted,
                                                     .data$x.adj - xoffset,
                                                     x),
             v.adj = if_else(.data$d.adjusted,
                             (v + .data$d) * exp(.data$c),
                             v * exp(.data$c) + .data$d))
  }
  return(list(x = abc$x.adj, y = abc$y.adj, v = abc$v.adj))
}
