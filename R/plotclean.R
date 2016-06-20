#' Plot multiple growth curves to identify outliers
#' 
#' A version of \code{mplot} to plot growth curves and identify outliers. When
#' outliers are clicked on, and if id is specified, the corresponding growth
#' curve is highlighted.  If id is not specified the selected point is
#' highlighted. Use right-click to exit.
#' 
#' 
#' @param x vector of x coordinates.
#' @param y vector of y coordinates.
#' @param id factor of subject levels indexing each growth curve.
#' @param data optional dataframe containing \code{x}, \code{y} and \code{id}.
#' @param n maximum number of points to be identified.
#' @param par.out list of optional graphical parameters to control appearance
#' of selected outlying points and lines.
#' @param \dots Further graphical parameters (see \code{\link{par}}) may also
#' be supplied as arguments for lines and points, particularly line type, lty,
#' line width, lwd and color, col.
#' @return \code{plotclean} returns either a vector \code{rows} (if data is not
#' specified) or a list: \item{rows}{a vector of row numbers corresponding to
#' the selected points.} \item{data}{a subset of \code{data} consisting of rows
#' \code{rows}, and columns \code{id}, \code{x} and \code{y}.}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' 
#' if (interactive()) plotclean(age, height, id, heights)
#' 
#' @export plotclean
	plotclean <- function(x, y, id=NULL, data=parent.frame(), n=length(x), par.out=list(pch=20), ...) {
#	plot growth curves to identify outlying points and curves
#	plot y ~ x by id with data
#	identify up to n points
#	par.out passes graphical par args for outliers
	mcall <- match.call()
	dots <- as.list(match.call(expand.dots = FALSE)$...)
	data.null <- identical(data, parent.frame())
	data <- eval(mcall$data)
	xt <- substitute(x)
	xlab <- deparse(xt)
	if (!"xlab" %in% names(dots)) dots <- c(dots, list(xlab=xlab))
	x <- eval(mcall$x, data)
	yt <- substitute(y)
	ylab <- deparse(yt)
	if (!"ylab" %in% names(dots)) dots <- c(dots, list(ylab=ylab))
	y <- eval(mcall$y, data)
	xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
	id <- eval(mcall$id, data)
	if (!"pch" %in% names(dots)) dots <- c(dots, list(pch=46))
	if (is.null(id)) {
		idlab <- NULL
		do.call('plot', c(list(x, y), dots))
	}
	else {
		idlab <- deparse(mcall$id)
		id <- factor(eval(mcall$id, data))
		do.call('plot', c(list(x, y, type='n'), dots))
		if (!"col" %in% names(dots)) dots <- c(dots, list(col='gray'))
		if (!"type" %in% names(dots)) dots <- c(dots, list(type='o'))
		tt <- by(data, id, function (z) do.call('lines', c(list(eval(yt, z) ~ eval(xt, z)), dots)))
	}
	sel <- rep(FALSE, length(x))
	res <- integer(0)
	title('click on outliers in plot - then right-click to escape')
	while (sum(sel) < n) {
		ans <- identify(x[!sel], y[!sel], n=1, plot=FALSE)
		if (!length(ans)) break
		ans <- which(!sel)[ans]
		do.call('points', c(list(x[ans], y[ans]), par.out))
		if (!is.null(id)) {
			text(x[ans], y[ans], label=id[ans], adj=c(0.5, 0.1))
			idt <- id==id[ans]
			ox <- order(x[idt])
			do.call('lines', c(list(x[idt][ox], y[idt][ox]), par.out))
		}	
		sel[ans] <- TRUE
		res <- c(res, ans)
	}
	res <- res[order(res)]
	if (data.null) res 
	else list(rows=res, data=data[res, c(idlab, xlab, ylab)])
}
