#' Identify peak or trough of curve
#'
#' Given vectors \code{x} and  \code{y}, returns their values at the (highest)
#' peak or (lowest) trough of the curve \code{y ~ x}, where \code{Dy} is zero and
#' \code{D2y} is negative or positive respectively.
#'
#' Note that turning points where the curvature is close to zero (and the curve
#' close to linear) are ignored.
#'
#' @param x vector.
#' @param y vector.
#' @param peak logical determining whether peak or trough is returned.
#' @return A length-2 vector containing the values of \code{x} and  \code{y}
#' at the peak or trough. If no peak/trough is identified \code{x} and  \code{y}
#' are set to NA.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' ## create mean height velocity curve
#' data(heights)
#' m1 <- sitar(age, height, id, heights, 4)
#' xy <- plot_v(m1)
#' ## and plot it
#' plot(xy, type='l', xlab='age', ylab='height velocity')
#' points(t(getPeak(xy)), pch=17)
#' points(t(getTrough(xy)), pch=25, bg=1)
#' @export getPeakTrough
	getPeakTrough <- function(x, y=NULL, peak=TRUE) {
#	returns values of x and y at peak/trough, i.e. where dy/dx=0
	xy <- xy.coords(x, y)
	ox <- order(xy$x)
	x <- xy$x[ox]
	y <- xy$y[ox]
  sign_change <- ifelse(peak, -1, +1)
	tp <- which(c(FALSE, diff(diff(y) / diff(x) > 0) == sign_change, FALSE)) # turning point(s)
	if (length(tp) == 0)
	  return(c(x=NA, y=NA))
# turning point(s) found
  d2y <- c(0, diff(diff(y) / diff(x)), 0) # curvature
  d2y <- abs(d2y) / max(abs(d2y)) # scale curvature
  tp <- tp[d2y[tp] > 0.1] # exclude noisy linear segments
  if (length(tp) == 0)
    return(c(x=NA, y=NA))
  if (length(tp) > 1) {
    best <- y[tp]
    best <- ifelse(peak, which.max(best), which.min(best)) # point of tp
    tp <- tp[best]
  }
	region <- max(1, tp - 2):min(tp + 2, length(x)) # region of tp
	x <- x[region]
	y <- y[region]
	curve <- lm(y ~ poly(x, 2, raw=TRUE)) # quadratic in x
	x <- - curve$coef[2] / curve$coef[3] / 2 # x at tp
	y <- predict(curve, data.frame(x=x)) # y at tp
	setNames(c(x, y), c('x', 'y'))
}
#' @rdname getPeakTrough
#' @export
	getPeak <- function(x, y=NULL, peak=TRUE) {}

#' @rdname getPeakTrough
#' @export
	 getTrough <- function(x, y=NULL, peak=FALSE) {}

	 body(getTrough) <- body(getPeak) <- body(getPeakTrough)

