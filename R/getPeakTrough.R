#' Identify peak or trough on curve
#'
#' Given vectors \code{x} and \code{y}, returns their values at the peak or
#' trough of the smooth (e.g. cubic spline) curve \code{y ~ x}.
#'
#' Optionally the trough can be specified as takeoff, which is defined
#' for a growth velocity curve as the lowest velocity before the pubertal peak,
#' and if there is no peak then there is by definition no takeoff.
#'
#' @param x vector.
#' @param y vector.
#' @param peak logical determining whether peak or trough is returned.
#' @param takeoff logical determining whether, if \code{peak} is FALSE, the
#' trough is takeoff.
#' @return A length-2 vector containing the values of \code{x} and \code{y} at
#' the peak or trough. If none are identified NA's are returned.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' ## create mean height velocity curve
#' data(heights)
#' m1 <- sitar(age, height, id, heights, 4)
#' ## plot velocity curve
#' plot(m1, 'v')
#' ## mark peak, trough and takeoff
#' xy <- plot_v(m1)
#' points(t(getPeak(xy)), pch=17)
#' points(t(getTrough(xy)), pch=25, col=2, bg=2)
#' points(t(getTakeoff(xy)), pch=25, col=3, bg=3)
#' @export
getPeakTrough <- function(x, y = NULL, peak = TRUE, takeoff = FALSE) {
  xy <- xy.coords(x, y)
  # sort data and remove duplicates
  xy <- unique(as.data.frame(xy[1:2])[order(xy$x),])
  x <- xy$x
  y <- xy$y
  # find peak
  ddy <- diff(diff(y) > 0)
  # turning point(s)
  tp <- which(ddy == -1) + 1
  # find highest peak
  tp <- tp[which.max(y[tp])]
  # return if no peaks and peak or takeoff
  if (length(tp) == 0 && (peak || takeoff))
    return(c(x = NA, y = NA))
  # find trough(s)
  if (!peak) {
    # if takeoff, search earlier than peak
    if (!takeoff)
      tp <- length(ddy)
    tp <- which(ddy[seq_len(tp)] == 1) + 1
    # find lowest trough
    tp <- tp[which.min(y[tp])]
    # return if no trough
    if (length(tp) == 0)
      return(c(x = NA, y = NA))
  }
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
#' @rdname getPeakTrough
#' @export
getPeak <- function(x, y = NULL, peak = TRUE, takeoff = FALSE) {

}
#' @rdname getPeakTrough
#' @export
getTrough <- function(x, y = NULL, peak = FALSE, takeoff = FALSE) {

}
#' @rdname getPeakTrough
#' @export
getTakeoff <- function(x, y = NULL, peak = FALSE, takeoff = TRUE) {

}
body(getTakeoff) <- body(getTrough) <- body(getPeak) <- body(getPeakTrough)
