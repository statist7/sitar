#' Identify peak or takeoff point of velocity curve
#'
#' Given vectors \code{x} and  \code{y}, returns their values at the (highest)
#' peak or (lowest) takeoff point of the smooth (e.g. cubic spline) velocity
#' curve \code{y ~ x}, where \code{Dy} is zero and \code{D2y} is
#' negative or positive respectively. The takeoff point is required to occur
#' before the peak.
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
#' ## plot it
#' plot(m1, 'v')
#' ## mark peak and takeoff
#' xy <- plot_v(m1)
#' points(t(getPeak(xy)), las=1, pch=17)
#' points(t(getTakeoff(xy)), pch=25, col=2, bg=2)
#' @export getPeaktakeoff
	getPeakTakeoff <- function(x, y=NULL, peak=TRUE) {
#	returns values of x and y at peak/takeoff, i.e. where dy/dx=0
	xy <- xy.coords(x, y)
# sort data and remove duplicates
	xy <- unique(as.data.frame(xy[1:2])[order(xy$x), ])
	x <- xy$x
	y <- xy$y
  sign_change <- ifelse(peak, -1, +1)
  tp <- which(c(FALSE, diff(diff(y) / diff(x) > 0) == sign_change, FALSE)) # turning point(s)
  if (length(tp) == 0)
	  return(c(x=NA, y=NA))
# turning point(s) found
  if (length(tp) > 1) {
# multiple turning points - find the best
    best <- y[tp]
    best <- ifelse(peak, which.max(best), which.min(best)) # point of tp
    tp <- tp[best]
  }
	region <- max(1, tp - 2):min(tp + 2, length(x)) # region of tp
	x <- x[region]
	y <- y[region]
	curve <- lm(y ~ poly(x, 2, raw=TRUE)) # quadratic in x
	x <- - curve$coef[[2]] / curve$coef[[3]] / 2 # x at tp
	y <- unname(predict(curve, data.frame(x=x))) # y at tp
	return(c(x=x, y=y))
	}
#' @rdname getPeakTakeoff
#' @export
	getPeak <- function(x, y=NULL, peak=TRUE, use="value") {}

#' @rdname getPeakTakeoff
#' @export
	 getTakeoff <- function(x, y=NULL, peak=FALSE, use="value") {}

	 body(getTakeoff) <- body(getPeak) <- body(getPeakTakeoff)

