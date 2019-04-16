#' Identify peak or trough of curve
#'
#' Given vectors \code{x} and  \code{y}, returns their values at the (highest)
#' peak or (lowest) trough of the smooth curve \code{y ~ x} (typically a fitted
#' cubic spline curve), where \code{Dy} is zero and \code{D2y} is negative or
#' positive respectively.
#'
#' @param x vector.
#' @param y vector.
#' @param peak logical determining whether peak or trough is returned.
#' @param use character string defining the criterion to use - "value"
#' (default) or "curvature".
#' @return A length-2 vector containing the values of \code{x} and \code{y} at
#' the peak or trough, defined as the maximum/minimum value or minimum/maximum
#' curvature, depending on \code{use}. If no peak/trough is identified NA's
#' are returned.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' ## create mean height velocity curve
#' data(heights)
#' m1 <- sitar(age, height, id, heights, 4)
#' ## plot it
#' plot(m1, 'v')
#' ## mark peak and trough
#' xy <- plot_v(m1)
#' points(t(getPeak(xy)), las=1, pch=17)
#' points(t(getTrough(xy)), pch=25, col=2, bg=2)
#' ## wrong trough - use curvature criterion instead
#' points(t(getTrough(xy, use='curvature')), pch=25, bg=1)
#' @export getPeakTrough
	getPeakTrough <- function(x, y=NULL, peak=TRUE, use="value") {
#	returns values of x and y at peak/trough, i.e. where dy/dx=0
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
    if (use == 'curvature') {
    	yt <- -c(NA, diff(diff(y) / diff(x)) / diff(x, 2), NA) # minus curvature
    } else if (use == 'value') {
      yt <- y # value
    } else
      return(c(x=NA, y=NA))
    best <- yt[tp]
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
#' @rdname getPeakTrough
#' @export
	getPeak <- function(x, y=NULL, peak=TRUE, use="value") {}

#' @rdname getPeakTrough
#' @export
	 getTrough <- function(x, y=NULL, peak=FALSE, use="value") {}

	 body(getTrough) <- body(getPeak) <- body(getPeakTrough)

