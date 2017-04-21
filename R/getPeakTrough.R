#' Identify peaks and troughs of curve
#'
#' Given vectors \code{x} and  \code{y}, returns their values at the peak or
#' trough of the curve, where dy/dx = 0.
#'
#' @param x vector.
#' @param y vector.
#' @param peak logical determining whether peak or trough is returned.
#' @return A length-2 vector containing the values of \code{x} and  \code{y}
#' at the peak or trough. If no peak/trough is identified NULL is returned invisibly.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' ## create mean height velocity curve
#' data(heights)
#' m1 <- sitar(age, height, id, heights, 4)
#' x <- getCovariate(m1)
#' y <- fitted(m1, level=0)
#' y <- predict(smooth.spline(x, y), x, deriv=1)$y
#'
#' ## and plot it
#' . <- order(x)
#' plot(y[.] ~ x[.], type='l', xlab='age', ylab='height')
#' points(t(getPeakTrough(x, y)), pch=17)
#' points(t(getPeakTrough(x, y, peak=FALSE)), pch=25)
#' @export getPeakTrough
	getPeakTrough <- function(x, y, peak=TRUE) {
#	returns values of x and y at peak/trough, i.e. where dy/dx=0
	. <- unique(data.frame(x, y)[order(x), ])
	x <- .$x
	y <- .$y
  . <- ifelse(peak, -1, 1)
	. <- c(FALSE, diff(diff(y) > 0) == ., FALSE) # find turning point(s)
	if (any(., na.rm=TRUE)) { # tp(s) found
	  . <- y * .
	  .[!.] <- NA
	  . <- ifelse(peak, which.max(.), which.min(.)) # point of tp
		. <- max(1, . - 2):min(. + 2, length(x)) # region of tp
		x <- x[.]
		y <- y[.]
		. <- lm(y ~ poly(x, 2, raw=TRUE)) # quadratic in x
		x <- - .$coef[2] / .$coef[3] / 2 # x at tp
		y <- predict(., data.frame(x=x)) # y at tp
		setNames(c(x, y), c('x', 'y'))
	} else
	  invisible(NULL)
}
