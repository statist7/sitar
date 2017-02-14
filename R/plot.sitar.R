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
#' Details). The latter take precedence. Default labels are the names of \code{x} and \code{y}, and
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
#' @param ns scalar defining the number of points for spline curves
#' (default 101).
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
#' @importFrom graphics axis identify legend lines locator par text title mtext abline
#' @export
plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun=NULL, yfun=NULL, subset=NULL,
	                     ns=101, abc=NULL, add=FALSE, nlme=FALSE, ...) {
	if (nlme) {
	  do.call('plot.lme', as.list(match.call()[-1]))
	}
	else {
    options <- c('d', 'e', 'u', 'a', 'D', 'v', 'V')
    optaxis <- c( 1,   1,   1,   1,   1,   2,   2 ) # default y1=1, y2=2
    optmult <- c( F,   F,   T,   T,   T,   F,   T ) # multiple curves
    opts <- na.omit(match(unlist(strsplit(opt, '')), options))
    if (length(opts) == 0) stop('option(s) not recognised')
    y2 <- diff(range(optaxis[opts])) == 1 # need both y1 and y2
    mult <- any(optmult[opts]) # multiple curves
# identify axis ranges then draw axes
    # axismin <- 3; axismax <- 0
    # for (i in 1:nchar(opt)) {
    #   no <- match(substr(opt, i, i), options, NA)
    #   if (is.na(no)) next
    #   if (optaxis[no] > axismax) axismax <- optaxis[no]
    #   if (optaxis[no] < axismin) axismin <- optaxis[no]
    # }
    # if (axismin == 2) {
    #   optaxis <- optaxis - 1
    #   axismax <- axismin <- 1
    # }
    model <- x
		data <- getData(model)
		mcall <- model$call.sitar
		x <- getCovariate(model)
		y <- getResponse(model)
		id <- getGroups(model)
		nf <- length(fitted(model))
		if (nf != length(y))
		  stop(paste0('model (length=', nf, ') incompatible with data (nrows=', length(y), ')'))
#	extract list(...)
		ccall <- match.call(expand.dots=FALSE)
#	subset to plot model
		subset <- eval(ccall$subset, data, parent.frame())
		if (is.null(subset)) subset <- rep_len(TRUE, nf)
# ... args
		dots <- ccall$...
		ARG <- if(!is.null(dots))
		  lapply(as.list(dots), eval, data, parent.frame())
		else
		  NULL
# create missing labels
		if (missing(labels))
		  labels <- vector('character', 3)
		else if (length(labels) < 3)
		  labels <- c(labels, '', '')
# test for velocity plot
		optv <- all(unique(strsplit(tolower(opt), '')[[1]]) == 'v')
# test for velocity plot label via y2par
		if (!is.null(ARG$y2par$ylab)) {
		  labels[3] <- ARG$y2par$ylab
		  if (optv)
		    ARG$ylab <- labels[2] <- ARG$y2par$ylab
		}
#	if xlab not specified replace with label or x name (depending on xfun)
		if (is.null(ARG$xlab)) {
		  ARG$xlab <- if(labels[1] != '')
		    labels[1]
		  else {
		    if (!is.null(xfun))
		      paste0('(', deparse(substitute(xfun)), ')(', deparse(mcall$x), ")")
		    else
		      attr(ifun(mcall$x), 'varname')
		  }
		} else
		  labels[1] <- ARG$xlab
#	if ylab not specified replace with label or y name (depending on yfun)
		if (is.null(ARG$ylab)) {
		  if(labels[2] != '')
		    ARG$ylab <- labels[2]
		  else {
		    ARG$ylab <- if (!is.null(yfun))
		      paste0('(', deparse(substitute(yfun)), ')(', deparse(mcall$y), ")")
		    else
		      attr(ifun(mcall$y), 'varname')
	      labels[2] <- ARG$ylab
		  }
	    if (labels[3] == '')
	      labels[3] <- paste(ARG$ylab, 'velocity')
		} else {
   	  labels[2] <- ARG$ylab
   	  if (optv)
   	    labels[3] <- ARG$ylab
  		else
  		  if (labels[3] == '')
  		    labels[3] <- paste(ARG$ylab, 'velocity')
   	}

#	create output list
		xy <- list()

# derive xfun and yfun
		if (is.null(xfun))
		  xfun <- ifun(mcall$x)
		if (is.null(yfun))
		  yfun <- ifun(mcall$y)

#	plot y vs t by subject
		if (grepl("u", opt)) {
		  xt <- x
		  yt <- y
		  do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG))
		  add <- TRUE
		}

		xseq <- function(x, n=ns) {
		  # n is the number of points across the x range
		  rx <- range(x, na.rm=TRUE)
		  seq(rx[1], rx[2], length.out=n)
		}

		stackage <- function(x, id, n=ns) {
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

#	plot fitted distance or velocity curves by subject
		for (o in c('D', 'V')) {
		  oV <- as.numeric(o == 'V')
  		if (grepl(o, opt)) {
  		  newdata <- stackage(x[subset], id[subset])
  		  newdata <- cbind(newdata, y=predict(model, newdata=newdata, xfun=xfun, yfun=yfun, deriv=oV))
  		  ARG$ylab <- labels[2 + oV]
# adjust ARG=id to ARG=newid
  		  for (i in 1:length(ARG)) {
  		    arg <- ARG[[i]]
  		    if (length(arg) == length(id) &&
  		      length(lvl <- unlist(tapply(id, arg, unique))) == nlevels(id))
  		        ARG[[i]] <- lvl[newdata[, 2]]
  		  }
  		  do.call("mplot", c(list(x=xfun(newdata[, 1]), y=newdata[, 3], id=newdata[, 2],
  		                          data=newdata, add=add), ARG))
  		  add <- TRUE
  		}
		}

#	plot fitted distance and velocity curves
		if (grepl("d", opt) || grepl("v", opt) || apv) {
			xt <- xseq(x[subset])
  		newdata <- data.frame(x=xt)
# if subset, flag for predict
      if (!identical(subset, rep_len(TRUE, nf))) attr(newdata, 'subset') <- subset

# distance and velocity
			xt <- xfun(xt)
			yt <- yfun(predict(object=model, newdata=newdata, level=0, abc=abc))
			vt <- predict(object=model, newdata=newdata, level=0, deriv=1, abc=abc, xfun=xfun, yfun=yfun)

# plot curve(s)
			if (grepl("d", opt) && grepl("v", opt)) {
				xy <- do.call("y2plot", c(list(x=xt, y1=yt, y2=vt, labels=labels, add=add, xy=xy), ARG))
			} else
			if (grepl("d", opt)) {
				xy <- do.call("y2plot", c(list(x=xt, y1=yt, add=add, xy=xy), ARG))
			} else
			if (grepl("v", opt)) {
				ARG$ylab <- labels[3]
				xy <- do.call("y2plot", c(list(x=xt, y1=vt, labels=labels[c(1,3)], add=add, xy=xy), ARG))
			}
#	plot vertical line at age of peak velocity
  		if (apv) {
  			xy$apv <- getPeakTrough(xt, vt)
  			if (!is.na(opt)) print(signif(xy$apv, 4))
				if (is.null(ARG$y2par$lty)) ARG$y2par$lty <- 3
				do.call('abline', c(list(v=xy$apv[1]), ARG$y2par))
  		}
			add <- TRUE
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
			yt <- xyadj(model)
			xt <- yt$x
			yt <- yt$y
    	do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG))
			add <- TRUE
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

