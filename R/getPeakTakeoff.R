#' Identify peak or takeoff point of velocity curve
#'
#' Given vectors \code{x} and \code{y}, returns their values at the (highest)
#' peak or (lowest) takeoff point of the smooth (e.g. cubic spline) velocity
#' curve \code{y ~ x}.
#'
#' The takeoff point is defined as the lowest velocity before the peak, so if
#' the curve has no peak there is no takeoff.
#'
#' @param x vector.
#' @param y vector.
#' @param peak logical determining whether peak or takeoff is returned.
#' @return A length-2 vector containing the values of \code{x} and \code{y} at
#' the peak or takeoff. If no peak/takeoff is identified NA's are returned.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' ## create mean height velocity curve
#' data(heights)
#' m1 <- sitar(age, height, id, heights, 4)
#' ## plot velocity curve
#' plot(m1, 'v')
#' ## mark peak and takeoff
#' xy <- plot_v(m1)
#' points(t(getPeak(xy)), pch=17)
#' points(t(getTakeoff(xy)), pch=25, col=2, bg=2)
#' @export getPeakTakeoff
getPeakTakeoff <- function(x, y = NULL, peak = TRUE) {
  xy <- xy.coords(x, y)
  # sort data and remove duplicates
  xy <- unique(as.data.frame(xy[1:2])[order(xy$x),])
  x <- xy$x
  y <- xy$y
  # find peak
  ddy <- diff(diff(y) > 0)
  # turning point(s)
  tp <- which(ddy == -1) + 1
  # return if no peak
  if (length(tp) == 0)
    return(c(x = NA, y = NA))
  # find highest peak
  tp <- tp[which.max(y[tp])]
  # find takeoff(s) earlier than peak
  if (!peak) {
    tp <- which(ddy[seq_len(tp)] == 1) + 1
    # find lowest takeoff
    tp <- tp[which.min(y[tp])]
  }
  # return if no takeoff
  if (length(tp) == 0)
    return(c(x = NA, y = NA))
# quadratic in x
  n <- 0
  repeat {
    n <- n + 1
    if (tp == n || tp + n > nrow(xy))
      break
    curve <- with(xy[(tp - n):(tp + n),],
                  lm(y ~ poly(x, 2, raw = TRUE)))
    if (curve$rank == 3)
      break
  }
  # x and y at tp
  x <- -curve$coef[[2]] / curve$coef[[3]] / 2
  y <- unname(predict(curve, data.frame(x = x)))
  return(c(x = x, y = y))
}
#' @rdname getPeakTakeoff
#' @export
getPeak <- function(x, y = NULL, peak = TRUE) {

}
#' @rdname getPeakTakeoff
#' @export
getTakeoff <- function(x, y = NULL, peak = FALSE) {

}
body(getTakeoff) <- body(getPeak) <- body(getPeakTakeoff)
