#' Plot multiple growth curves
#' 
#' Function to plot multiple growth curves indexed by subject id.
#' 
#' The arguments \code{x}, \code{y} and \code{id} can be given as character
#' strings. The \code{\link{par}} parameters can be functions of vector
#' variables in \code{data}, e.g. to colour curves separately by \code{id} use:
#' \code{col = id}.
#' 
#' @param x vector of x coordinates.
#' @param y vector of y coordinates.
#' @param id factor denoting subject levels.
#' @param data optional dataframe containing \code{x}, \code{y} and \code{id}.
#' @param subset optional logical defining a subset of rows in \code{data}.
#' @param add optional logical defining whether the plot is pre-existing (TRUE)
#' or new (FALSE).
#' @param \dots Further graphical parameters (see \code{\link{par}}) may also
#' be supplied as arguments, particularly background colour \code{bg},
#' character expansion \code{cex}, colour \code{col}, line type \code{lty},
#' line width \code{lwd} and character \code{pch}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{y2plot}}
#' @examples
#' 
#' mplot(age, height, id, heights, col=id)
#' 
#' @export mplot
	mplot <- function(x, y, id, data=parent.frame(), subset=NULL, add=FALSE, ...) {
#	plots y ~ x by id with data
#	x and y can be name or character
#	subset defines a subset of rows
#	add TRUE suppresses plot axes
#	... parameters where bg, cex, col, lty, lwd, pch can depend on id

#	save x y id
	mcall <- match.call()[-1]
	df <- as.data.frame(lapply(as.list(mcall[1:3]), function(z) {
		if (is.character(z)) with(data,	get(z, inherits=TRUE))
		else eval(z, envir = data, enclos = parent.frame())
	}))
	if (length(deparse(mcall)) == 1) {
    names(df) <- lapply(as.list(mcall[1:3]), function(z) {
  		if (is.character(z)) z else deparse(z)
  	})
	}

#	extract and save vector par args: bg cex col lty lwd pch
	if (length(dots <- match.call(expand.dots=FALSE)$...) > 0) {
		ARG <- lapply(as.list(dots), eval, envir = data, enclos = parent.frame())
		cnames <- names(ARG)[lapply(ARG, length) == nrow(df)]
		df[, cnames] <- ARG[cnames]
		ARG[cnames] <- NULL
	} else {
		ARG <- list()
		cnames <- NULL
	}

#	subset data
	subset <- eval(substitute(subset), data, parent.frame())
	if (!is.null(subset)) {
		if (length(subset) != nrow(df)) stop('subset wrong length for data')
		subset <- ifelse(is.na(df[, 1]) | is.na(df[, 2]), FALSE, subset)
		df <- df[subset, ]
	}
	if (nrow(df) == 0) stop("no data to plot")

#	plot axes if new graph
	if (!add) {
		if (!"xlab" %in% names(ARG)) ARG <- c(ARG, list(xlab=quote(names(df)[1])))
		if (!"ylab" %in% names(ARG)) ARG <- c(ARG, list(ylab=quote(names(df)[2])))
		type <- match(names(ARG), "type", 0)
		do.call("plot", c(list(x=df[, 1], y=df[, 2], type='n'), ARG[!type]))
	}

#	draw growth curves
	tt <- by(df, df[, 3], function(z) {
#	sort by x
		ox <- order(z[, 1])
#	restore vector ... args
		if (length(cnames) > 0) ARG[cnames] <- as.list(as.data.frame(z[ox, cnames], stringsAsFactors=FALSE))
#	lines(x, y, ...)
		do.call("lines", c(list(x=z[ox, 1], y=z[ox, 2]), ARG))
	})
}
