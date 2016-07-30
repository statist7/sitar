#' Plot SITAR model
#'
#' \code{plot} and \code{lines} methods for objects of class \code{sitar},
#' providing various flavours of plot of the fitted growth curves.
#'
#' For option 'dv' (the default) the velocity curve plot (with right axis) can
#' be annotated with \code{par} parameters given as a named list called
#' \code{y2par}. To suppress the legend that comes with it set \code{xlegend =
#' NULL}.
#'
#' @aliases plot.sitar lines.sitar
#' @param x object of class \code{sitar}.
#' @param opt character string containing a subset of letters corresponding to
#' the options: 'd' for fitted Distance curve, 'v' for fitted Velocity curve,
#' 'e' for fitted fixed Effects distance curve, 'D' for individual fitted
#' Distance curves, 'V' for individual fitted Velocity curves, 'u' for
#' Unadjusted individual growth curves, and 'a' for Adjusted individual growth
#' curves. Options 'dveDV' give spline curves, while 'ua' give data curves made
#' up as line segments. If both distance and velocity curves are specified, the
#' axis for the velocity curve appears on the right side of the plot (y2), and
#' a legend identifying the distance and velocity curves is provided.
#' @param labels optional character vector containing plot labels for \code{x},
#' \code{y} and \code{y} velocity from the original SITAR model. The first two
#' elements can alternatively be provided via \code{\link{par}} parameters
#' \code{xlab} and \code{ylab}, and the third element via \code{y2par} (see
#' Details). Default labels are the names of \code{x} and \code{y}, and
#' "\code{y} velocity", suitably adjusted to reflect any back-transformation
#' via \code{xfun} and \code{yfun}.
#' @param apv optional logical specifying whether or not to calculate the age
#' at peak velocity from the velocity curve. If TRUE, age at peak velocity is
#' calculated as the age when the second derivative of the fitted curve changes
#' sign (after applying \code{xfun} and/or \code{yfun}). Age at peak velocity
#' is marked in the plot with a vertical dotted line, and its value, along with
#' peak velocity, is printed and returned.
#' @param xfun optional function to be applied to the x variable prior to
#' plotting. Defaults to NULL, which translates to \code{ifun(x$call.sitar$x)}
#' and inverts any transformation applied to x in the original SITAR model
#' call. To plot on the transformed scale set \code{xfun} to \code{I}.
#' @param yfun optional function to be applied to the y variable prior to
#' plotting. Defaults to NULL, which translates to \code{ifun(x$call.sitar$y)}
#' and inverts any transformation applied to y in the original SITAR model
#' call. To plot on the transformed scale set \code{yfun} to \code{I}.
#' @param subset optional logical vector of length \code{x} defining a subset
#' of \code{data} rows to be plotted.
#' @param abc vector of named values of random effects a, b and c used to
#' define an individual growth curve, e.g. abc=c(a=1, c=-0.1). Alternatively a
#' single character string defining an \code{id} level whose random effect
#' values are used. If \code{abc} is set, \code{level} is ignored. If
#' \code{abc} is NULL (default), or if a, b or c values are missing, values of
#' zero are assumed.
#' @param add optional logical defining if the plot is pre-existing (TRUE) or
#' new (FALSE). TRUE is equivalent to using \code{lines}.
#' @param nlme optional logical which set TRUE plots the model as an
#' \code{nlme} object, using \code{plot.nlme} arguments.
#' @param \dots Further graphical parameters (see \code{par}) may also be
#' supplied as arguments, e.g. axis labels \code{xlab} and \code{ylab}, line
#' type \code{lty}, line width \code{lwd}, and color \code{col}. For the
#' velocity (y2) plot \code{y2par} can be used (see Details).
#' @return Returns invisibly a list of three objects:
#' \item{ss}{\code{smooth.spline} object corresponding to the fitted distance
#' curve.} \item{usr}{value of \code{par('usr')} for the main plot.}
#' \item{usr2}{the value of \code{par('usr')} for the velocity (y2) plot.}
#' \item{apv}{if argument \code{apv} is TRUE a named list giving the age at
#' peak velocity (apv) and peak velocity (pv) from the fitted velocity curve.}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{mplot}},
#' \code{\link{plotclean}}, \code{\link{y2plot}}, \code{\link{ifun}}
#' @keywords aplot
#' @examples
#'
#' ## fit sitar model
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=7)
#'
#' ## draw fitted distance and velocity curves
#' ## with velocity curve in blue
#' ## adding age at peak velocity
#' plot(m1, y2par=list(col='blue'), apv=TRUE)
#'
#' ## draw individually coloured growth curves adjusted for random effects
#' ## using same x-axis limits as for previous plot
#' plot(m1, opt='a', col=id, xlim=xaxsd())
#'
#' ## add mean curve in red
#' lines(m1, opt='d', col='red', lwd=2)
#'
#' ## add mean curve for a, b, c = -1 SD
#' lines(m1, opt='d', lwd=2, abc=-sqrt(diag(getVarCov(m1))))
#'
#' @importFrom grDevices xy.coords
#' @importFrom graphics axis identify legend lines locator par text title
#' @export
plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun=NULL, yfun=NULL, subset=NULL,
	                       abc=NULL, add=FALSE, nlme=FALSE, ...) {
	if (nlme) {
	  do.call('plot.lme', as.list(match.call()[-1]))
	}
	else {
    options <- c('d', 'e', 'u', 'a', 'D', 'v', 'V')
    optaxis <- c(1, 1, 1, 1, 1, 2, 2) # default y axis
    optmult <- c(F, F, T, T, T, F, T) # multiple curves?
    axismin <- 3; axismax <- 0
    for (i in 1:nchar(opt)) {
      no <- match(substr(opt, i, i), options, NA)
      if (is.na(no)) next
      if (optaxis[no] > axismax) axismax <- optaxis[no]
      if (optaxis[no] < axismin) axismin <- optaxis[no]
    }
    if (axismin == 2) {
      optaxis <- optaxis - 1
      axismax <- axismin <- 1
    }
    model <- x
		data <- getData(model)
		mcall <- model$call.sitar
		x <- getCovariate(model)
		y <- getResponse(model)
		id <- getGroups(model)
		nf <- length(fitted(model))
		if (nf != length(y)) stop(paste0('model (length=', nf, ') incompatible with data (rows=', length(y), ')'))

#	extract list(...)
		ccall <- match.call()[-1]
#	subset to plot model
		subset <- if (is.null(subset))
		  rep_len(TRUE, nf)
		else
		  eval(ccall$subset, data, parent.frame())
# ... args
		dots <- match.call(expand.dots=FALSE)$...
		ARG <- if(!is.null(dots))
		  lapply(as.list(dots), eval, data, parent.frame())
		else
		  NULL
#	if xlab not specified replace with label or x name (depending on xfun)
		if (is.null(ARG$xlab)) {
		  ARG$xlab <- if(!missing(labels))
		    labels[1]
		  else {
		    if(!is.null(xfun))
		      paste0('(', deparse(substitute(xfun)), ')(', deparse(mcall$x), ")")
		    else
		      ifun(mcall$x)$varname
		  }
		}
#	if ylab not specified replace with label or y name (depending on yfun)
		if (is.null(ARG$ylab)) {
		  ARG$ylab <- if(!missing(labels))
		    labels[2]
		  else {
		    if(!is.null(yfun))
		      paste0('(', deparse(substitute(yfun)), ')(', deparse(mcall$y), ")")
		    else
		      ifun(mcall$y)$varname
		  }
		}
# if labels not specified create it
		if (missing(labels))
		  labels <- c(ARG$xlab, ARG$ylab, paste(ARG$ylab, 'velocity'))

#	create output list
		xy <- list()

# derive xfun and yfun
		if (is.null(xfun))
		  xfun <- ifun(mcall$x)$fn
		if (is.null(yfun))
		  yfun <- ifun(mcall$y)$fn

#	plot y vs t by subject
		if (grepl("u", opt)) {
		  xt <- x
		  yt <- y
		  do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG))
		  add <- TRUE
		}

		xseq <- function(x, n=101) {
		  # n is the number of points across the x range
		  rx <- range(x, na.rm=TRUE)
		  seq(rx[1], rx[2], length.out=n)
		}

		stackage <- function(x, id, n=101) {
		  # generate x and id values across the x range to plot spline curves
		  npt <- n / diff(range(x))
		  xid <- by(data.frame(x=x, id=id), id, function(z) {
		    nt <- floor(npt * diff(range(z$x))) + 1
		    data.frame(x=seq(min(z$x), to=max(z$x), length.out=nt), id=rep.int(z$id[[1]], nt))
		  })
		  df <- xid[[1]][FALSE, ]
		  for (dft in xid) df <- rbind(df, dft)
		  df
		}

#	plot fitted curves by subject
		if (grepl("D", opt)) {
		  newdata=stackage(x[subset], id[subset])
		  newdata <- cbind(newdata, y=predict(model, newdata=newdata, xfun=xfun, yfun=yfun))
		  do.call("mplot", c(list(x=xfun(newdata[, 1]), y=newdata[, 3], id=newdata[, 2],
		                          data=newdata, add=add), ARG))
		  add <- TRUE
		}

#	plot fitted velocity curves by subject
		if (grepl("V", opt)) {
		  newdata=stackage(x[subset], id[subset])
		  newdata <- cbind(newdata, y=predict(model, newdata=newdata, deriv=1, xfun=xfun, yfun=yfun))
		  ARG$ylab <- labels[3]
		  do.call("mplot", c(list(x=xfun(newdata[, 1]), y=newdata[, 3], id=newdata[, 2],
                              data=newdata, add=add), ARG))
      add <- TRUE
		}

#	plot fitted distance and velocity curves
		if (grepl("d", opt) || grepl("v", opt) || apv) {
			xt <- xseq(x[subset])
  		newdata <- data.frame(x=xt)
# if subset, flag for predict
      if (!identical(subset, rep(TRUE, nf))) attr(newdata, 'subset') <- subset

#	adjust for abc
			if (!is.null(abc)) {
#	abc is id level
				if (length(abc) == 1) {
					idabc <- rownames(ranef(model)) %in% abc
					if (sum(idabc) == 0) stop(paste('id', abc, 'not found'))
					abc <- ranef(model)[idabc, ]
				}
#	abc is named vector
			  else if (length(abc) > 3 || is.null(names(abc)))
			    stop('abc should be either single id level or up to three named random effect values')
			}
  		else abc <- ranef(model)

			yt <- yfun(predict(object=model, newdata=newdata, level=0, abc=abc))
			vt <- predict(object=model, newdata=newdata, level=0, deriv=1, abc=abc, xfun=xfun, yfun=yfun)
#	derive cubic smoothing spline curve
			xt <- xfun(xt)
			xy$ss <- ss <- makess(xt, yt)

			if (grepl("d", opt) && grepl("v", opt)) {
				xy <- do.call("y2plot", c(list(x=xt, y1=yt, y2=vt, labels=labels, add=add, xy=xy), ARG))
				add <- TRUE
			} else
			if (grepl("d", opt)) {
				xy <- do.call("y2plot", c(list(x=xt, y1=yt, add=add, xy=xy), ARG))
				add <- TRUE
			} else
			if (grepl("v", opt)) {
				ARG$ylab <- labels[3]
				xy <- do.call("y2plot", c(list(x=xt, y1=vt, labels=labels[c(1,3)], add=add, xy=xy), ARG))
				add <- TRUE
			}
		}

#	plot fixed effects distance curve
		if (grepl("e", opt)) {
  		xt <- xseq(x[subset])
      yt <- predict(model$ns, newdata=data.frame(x=xt - model$xoffset))
			ox <- order(xt)
			xy <- do.call("y2plot", c(list(x=xfun(xt[ox]), y1=yfun(yt[ox]), add=add, xy=xy), ARG))
			add <- TRUE
		}

#	plot adjusted y vs adjusted t by subject
		if (grepl("a", opt)) {
			yt <- xyadj(x, y, id, model)
			xt <- yt$x
			yt <- yt$y
    	do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG))
			add <- TRUE
		}

#	plot vertical line at age of peak velocity
		if (apv) {
			xy$apv <- ss$apv
			if (!is.na(opt)) print(signif(xy$apv, 4))
			if (add) {
				if (is.null(ARG$y2par$lty)) ARG$y2par$lty <- 3
				do.call('abline', c(list(v=xy$apv[1]), ARG$y2par))
			}
		}
		invisible(xy)
	}
}
#############################
#
#	lines.sitar
#
#############################

#' @rdname plot.sitar
#' @export
lines.sitar <- function (x, ...)
{
  mcall <- match.call()
  mcall[[1]] <- as.name("plot")
  mcall[['add']] <- TRUE
  eval(mcall, parent.frame())
}

