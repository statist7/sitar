#	global variables
	if (getRversion() >= "3.0.0") utils::globalVariables(c('.par.usr2', 'fitnlme', 'fitenv'))

#############################
#
#	print.sitar
#
#############################





#' Print SITAR model
#'
#' Print method for \code{sitar} objects, based on \code{print.lme}.
#'
#'
#' @param x an object inheriting class \code{sitar}.
#' @param \dots other optional arguments.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @export
	print.sitar <- function (x, ...)
{
    dd <- x$dims
    cat("SITAR nonlinear mixed-effects model fit by ")
    cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
    cat("  Call:", deparse(x$call.sitar), fill=TRUE)
    cat("  Log-", ifelse(x$method == "REML", "restricted-", ""),
        "likelihood: ", format(x$logLik), "\n", sep = "")
    fixF <- x$call$fixed
    if (inherits(fixF, "formula") || is.call(fixF) || is.name(fixF)) {
        cat("  Fixed:", deparse(x$call$fixed), "\n")
    }
    else {
        cat("  Fixed:", deparse(lapply(fixF, function(el) as.name(deparse(el)))),
            "\n")
    }
    print(fixef(x))
    cat("\n")
    print(summary(x$modelStruct), sigma = x$sigma)
    cat("Number of Observations:", dd[["N"]])
    cat("\nNumber of Groups: ")
    Ngrps <- dd$ngrps[1:dd$Q]
    if ((lNgrps <- length(Ngrps)) == 1) {
        cat(Ngrps, "\n")
    }
    else {
        sNgrps <- 1:lNgrps
        aux <- rep(names(Ngrps), sNgrps)
        aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps,
            lNgrps))[!lower.tri(diag(lNgrps))])
        names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
        cat("\n")
        print(rev(Ngrps))
    }
    invisible(x)
}

#############################
#
#	summary.sitar
#
#############################





#' Create summary of SITAR model
#'
#' A \code{summary} method for \code{sitar} objects based on
#' \code{\link{summary.lme}}.
#'
#'
#' @param object object inheriting from class \code{sitar}.
#' @param adjustSigma optional logical (see \code{\link{summary.lme}}).
#' @param verbose optional logical to control the amount of output in
#' \code{print.summary.sitar}.
#' @param \dots some methods for this generic require additional arguments.
#' None are used in this method.
#' @return an object inheriting from class \code{summary.sitar} with all
#' components included in \code{object} (see \code{\link{lmeObject}} for a full
#' description of the components) plus the components for
#' \code{\link{summary.lme}} and the following components: \item{x.adj}{vector
#' of length \code{x} in \code{object} with \code{x} values adjusted for
#' subject-specific random effects b and c.} \item{y.adj}{vector of length
#' \code{y} in \code{object} with \code{y} values adjusted for subject-specific
#' random effects a.} \item{apv}{length 2 vector giving respectively age at
#' peak velocity and peak velocity based on the fitted distance curve (using
#' transformed \code{x} and \code{y} where specified).}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @export
	summary.sitar <- function (object, adjustSigma = TRUE, verbose = FALSE, ...)
{
	  class(object) <- class(object)[-1]
    object <- summary(object, adjustSigma=adjustSigma, verbose=verbose, ...)
    class(object) <- c("summary.sitar", "sitar", class(object))
#	save age at peak velocity
    x <- getCovariate(object)
    y <- fitted(object, level=0)
    y <- predict(smooth.spline(unique(cbind(x, y))), x, deriv=1)$y
    object$apv <- getPeakTrough(x, y)
    object
}

#############################
#
#	print.summary.sitar
#
#############################





#' Print summary of SITAR model
#'
#' A \code{print.summary} method for \code{sitar} objects.
#'
#'
#' @param x an object inheriting from class \code{summary.sitar}.
#' @param verbose a logical to control the amount of output.
#' @param \dots to specify extra arguments.
#' @return A formatted summary of the object.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @method print summary.sitar
#' @export
	print.summary.sitar <- function (x, verbose = FALSE, ...)
{
    dd <- x$dims
    verbose <- verbose || attr(x, "verbose")
    cat("SITAR nonlinear mixed-effects model fit by ")
    cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
    cat("  Call:", deparse(x$call.sitar), fill=TRUE)
    print(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = c(x$logLik),
        row.names = " "))
    if (verbose) cat("\nNumber of iterations:", x$numIter, "\n")
    cat("\n")
    print(summary(x$modelStruct), sigma = x$sigma, reEstimates = x$coef$random,
        verbose = verbose)
    cat("Fixed effects: ")
    fixF <- x$call$fixed
    if (inherits(fixF, "formula") || is.call(fixF)) {
        cat(deparse(x$call$fixed), "\n")
    }
    else {
        cat(deparse(lapply(fixF, function(el) as.name(deparse(el)))),
            "\n")
    }
    xtTab <- as.data.frame(x$tTable)
    wchPval <- match("p-value", names(xtTab))
    for (i in names(xtTab)[-wchPval]) {
        xtTab[, i] <- format(zapsmall(xtTab[, i]))
    }
    xtTab[, wchPval] <- format(round(xtTab[, wchPval], 4))
    if (any(wchLv <- (as.double(levels(xtTab[, wchPval])) ==
        0))) {
        levels(xtTab[, wchPval])[wchLv] <- "<.0001"
    }
    row.names(xtTab) <- dimnames(x$tTable)[[1]]
    print(xtTab)
    if (nrow(x$tTable) > 1) {
        corr <- x$corFixed
        class(corr) <- "correlation"
        print(corr, title = " Correlation:", ...)
    }
    cat("\nStandardized Within-Group Residuals:\n")
    print(x$residuals)
    cat("\nNumber of Observations:", x$dims[["N"]])
    cat("\nNumber of Groups: ")
    Ngrps <- dd$ngrps[1:dd$Q]
    if ((lNgrps <- length(Ngrps)) == 1) {
        cat(Ngrps, "\n")
    }
    else {
        sNgrps <- 1:lNgrps
        aux <- rep(names(Ngrps), sNgrps)
        aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps,
            lNgrps))[!lower.tri(diag(lNgrps))])
        names(Ngrps) <- unlist(lapply(aux, paste, collapse = " %in% "))
        cat("\n")
        print(rev(Ngrps))
    }
    invisible(x)
}

#############################
#
#	y2plot
#
#############################





#' Plot with two y axes
#'
#' Function to plot two y variables, y1 and y2, against a single x variable,
#' with the y1 and y2 axes on the left and right of the plot respectively.
#'
#' y2plot draws up to two superimposed plots, one with the y axis on the left
#' and the other on the right, with suitable adjustment for \code{par('mar')}
#' and including a legend.  The format for y1 is controlled by par arguments,
#' and that for y2 by the list y2par.
#'
#' @param x vector of ages.
#' @param y1 vector of measurements for plotting on left y axis.
#' @param y2 optional vector of measurements for plotting on right y axis.
#' @param labels character vector containing labels for x, y1 and y2.
#' @param y2par optional named list of par arguments to format the y2 axis.
#' @param add logical flag to specify if a new plot (with axes etc) is to be
#' drawn (FALSE) or an existing plot is to be added to (TRUE).
#' @param xy optional list to pass \code{usr} and \code{usr2} (see Value).
#' @param xlegend position for legend.
#' @param inset inset for legend.
#' @param \dots optional \code{par} arguments.
#' @return Returns the list \item{usr}{par('usr') for y1 axis}
#' \item{usr2}{par('usr') for y2 axis} In addition the variable
#' \code{.par.usr2}, equal to \code{usr2}, is created in globalenv().
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{plot.sitar}}
#' @examples
#'
#' ## plot boys median height and weight on the UK 1990 reference
#' with(uk90[uk90$sex == 1,],
#'   y2plot(x=years, y1=M.ht, y2=M.wt, y2par=list(col='red'))
#' )
#'
#' @export y2plot
	y2plot <- function(x, y1, y2=NULL, labels, y2par=list(), add=FALSE, xy=NULL, xlegend="topleft", inset=0.04, ...)
#	plot x versus y1 and y2 using left/right axes for the two y's
#	returns par() for y1 or list(par(), par("usr")) with y2
#	labels are labels for x, y1 and y2 (xlab and ylab override x and y1, y2par['ylab'] overrides y2)
#	y2par is optional list of par args for y2 axis
#	add suppresses plot axes
#	xy set uses par() from previous call
#	xlegend and inset place legend (NULL suppresses)
{
#	save par args
	opar <- par(no.readonly=TRUE)
#	args for y1 and y2
	ypar <- list(...)
	y2par <- as.list(y2par)
# default dotted line for velocity curve
	if (is.null(y2par$lty)) y2par$lty <- 2
#	if a new plot draw axes
	if (!add) {
#	if mar specified then set it
		if (!is.null(ypar$mar)) {
			mar <- ypar$mar
		}
#	otherwise reset mar
		else {
			mar <- opar$mar
#	if y2 set increase mar[4]
			if (!missing(y2)) mar[4] <- mar[2]
		}
		par(mar=mar)
# get axis labels
  	if (missing(labels))
  	  labels <- c(deparse(substitute(x)), deparse(substitute(y1)), deparse(substitute(y2)))
		if (is.null(ypar$xlab))
		  ypar$xlab <- labels[1]
		if (is.null(ypar$ylab))
		  ypar$ylab <- labels[2]
#	plot x & y1 axes and y1 ~ x
		do.call('plot', c(list(x=quote(x), y=quote(y1), type='l'), ypar))
#	save x & y1 axis limits
		xy$usr <- par('usr')
#	optionally plot y2 and add right axis, with y2par as ...
		if (!missing(y2)) {
			par(new=TRUE)
#	ensure x axis same for y2 as y1
		  if (!is.null(ypar$xlim) && is.null(y2par$xlim)) y2par$xlim <- ypar$xlim
			do.call('plot', c(list(x=quote(x), y=quote(y2), ann=FALSE, axes=FALSE, type="l"), y2par))
#	save y2 axis limits
			xy$usr2 <- par('usr')
			eval(parse(text=".par.usr2 <<- par('usr')"))
# local functions
			localdots <- function(..., col, bg, pch, lty, lwd) list(...)
			localaxis <- function(..., col, bg, pch, lty, lwd) axis(...)
			localmtext <- function(..., col, bg, pch, lty, lwd, las) mtext(...)
# copy relevant args for y2 axis
			y2par <- c(y2par, do.call('localdots', c(ypar)))
#	add y2 axis
			do.call('localaxis', c(list(side=4), y2par))
#	add y2 axis label
			mgp <-  if (!is.null(ypar$mgp))
			  ypar$mgp
			else
			  par('mgp')
			do.call('localmtext', c(list(text=labels[3], side=4, line=mgp[1]), y2par))
#	add legend
			if (!is.null(xlegend) && !is.null(inset)) {
  			parlu <- function(...) par(do.call('par', list(...)))
  			llc <- sapply(c('lty', 'lwd', 'col'), function(i) {
  			  p12 <- rep(par()[[i]], 2)
  			  if (!is.null(ypar[[i]]))
  			    p12[1] <- parlu(ypar[i])[[1]]
  			  if (!is.null(y2par[[i]]))
  			    p12[2] <- parlu(y2par[i])[[1]]
  			  p12
  			})
			  legend(xlegend, legend=labels[2:3], bty="o",
			         lty=llc[, 'lty'], lwd=llc[, 'lwd'], col=llc[, 'col'], inset=inset)
			}
#	reset axis limits
			par(usr=xy$usr)
		}
	}
	else {
#	plot y1
		do.call('lines', c(list(x, y1), ypar))
#	optionally plot y2
		if (!missing(y2)) {
			if (exists('xy$usr2')) {
				par(usr=xy$usr2)
				eval(parse(text=".par.usr2 <<- par('usr')"))
			}
			else if (exists('.par.usr2')) {
				par(usr=.par.usr2)
				xy$usr2 <- par('usr')
			}
			else stop("Error in y2plot: second y axis requires previous call to set it\n", call.=FALSE)
			do.call('lines', c(list(x=x, y=y2), y2par))
			par(usr=opar$usr)
		}
	}
	invisible(xy)
}

#############################
#
#	funcall
#
#############################





#' Function call with optional inverse
#'
#' Applies an expression to vector v, optionally inverting the expression
#' first. For example if the expression is log, funcall returns log(v) if
#' inverse is FALSE, and exp(v) if inverse is TRUE.
#'
#' Inverse covers functions log, exp, sqrt, ^, *, /, +, -.
#'
#' @param v vector
#' @param vcall expression
#' @param inverse logical
#' @return Returns a vector of length v.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @export funcall
	funcall <- function(v, vcall, inverse=FALSE)
#	returns v transformed according to vcall
#	v - vector
#	vcall - expression
#	assumes name  of v is vcall or vcall[[2]]
#		excludes e.g. log10
#	if inverse TRUE vcall is inverted first
{
	if (length(vcall) == 1) vcall <- substitute(v) else {
		if (inverse) {
			fun <- vcall[[1]]
			if (fun == 'log') vcall[[1]] <- as.symbol('exp') else
			if (fun == 'exp') vcall[[1]] <- as.symbol('log') else
			if (fun == 'sqrt') {
				vcall[[1]] <- as.symbol('^')
				vcall[[3]]  <- 2
			} else
			if (fun == '^') vcall[[3]] <- 1/vcall[[3]] else
			if (fun == '*') vcall[[1]] <- as.symbol('/') else
			if (fun == '/') vcall[[1]] <- as.symbol('*') else
			if (fun == '+') vcall[[1]] <- as.symbol('-') else
			if (fun == '-') vcall[[1]] <- as.symbol('+') else
			stop ('funcall: unknown function')
		}
		vcall[[2]] <- substitute(v)
	}
	eval(vcall)
}


#############################
#
#	bupdate
#
#############################

#	update value of bstart to minimise b-c correlation




#' Update the b fixed effect to minimise the b-c random effect correlation
#'
#' A function to update the value of \code{bstart}, the starting value for the
#' b fixed effect, to minimise the correlation between the random effects b and
#' c.
#'
#'
#' @param x a \code{sitar} object.
#' @return Returns an updated value of the b fixed effect, based on the random
#' effect covariance matrix.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @keywords regression
#' @examples
#'
#' ## fit sitar model with b fixed effect starting value defaulting to 'mean'
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)
#' print(fixef(m1)['b'])
#'
#' ## refit with starting value chosen to minimise b-c correlation and df increased
#' m2 <- update(m1, bstart=bupdate(m1), df=6)
#' print(fixef(m2)['b'])
#'
#' @export bupdate
	bupdate <- function(x) {
	cov <- getVarCov(x)
	fixef(x)['b'] + x$xoffset + cov[2,3] / cov[3,3]
	}

#############################
#
#	velout
#
#############################





#' Identify outliers with abnormal velocity in growth curves
#'
#' Quickly identifies putative outliers in a large number of growth curves.
#'
#' The algorithm works by viewing serial measurements in each growth curve as
#' triplets (A-B-C) and comparing the velocities between them. Velocity is
#' calculated as
#'
#' diff(y, lag = lag) / diff(x, lag = lag) ^ velpower
#'
#' Missing values for x or y are ignored. If any of the AB, BC or AC velocities
#' are abnormal (more than \code{limit} SDs in absolute value from the median
#' for the dataset) the code for B is non-zero.
#'
#' @param x age vector.
#' @param y outcome vector, typically weight or height.
#' @param id factor identifying each subject.
#' @param data data frame containing x, y and id.
#' @param lag lag between measurements for defining growth velocity.
#' @param velpower a value, typically between 0 and 1, defining the power of
#' delta x to use when calculating velocity as delta(y)/delta(x)^velpower. The
#' default of 0.5 is midway between velocity and increment.
#' @param limit the number of standard deviations beyond which a velocity is
#' deemed to be an outlier.
#' @param linearise if TRUE y is converted to a residual about the median curve
#' of y versus x.
#' @return Returns a data frame with columns: id, x, y (from the call), code
#' (as described below), vel1, vel2 and vel3 (corresponding to the velocities
#' AB, BC and AC above). The 'data' attribute contains the name of 'data'.
#'
#' Code is a factor taking values between 0 and 8, with 0 normal (see table
#' below). Values 1-6 depend on the pattern of abnormal velocities, while 7 and
#' 8 indicate a duplicate age (7 for the first in an individual and 8 for later
#' ones). Edge outliers, i.e. first or last for an individual, have just one
#' velocity. Code 4 indicates a conventional outlier, with both AB and BC
#' abnormal and AC normal. Code 6 is an edge outlier. Other codes are not
#' necessarily outliers, e.g. codes 1 or 3 may be adjacent to a code 4. Use
#' \code{codeplot} to look at individual curves, and \code{zapvelout} to delete
#' outliers. \tabular{cccl}{ code \tab AB+BC \tab AC \tab interpretation\cr 0
#' \tab 0 \tab 0 \tab no outlier\cr 0 \tab 0 \tab NA \tab no outlier\cr 1 \tab
#' 0 \tab 1 \tab rare pattern\cr 2 \tab 1 \tab 0 \tab complicated - look at
#' curve\cr 3 \tab 1 \tab 1 \tab adjacent to simple outlier\cr 4 \tab 2 \tab 0
#' \tab single outlier\cr 5 \tab 2 \tab 1 \tab double outlier\cr 6 \tab 1 \tab
#' NA \tab edge outlier\cr 7 \tab - \tab - \tab first duplicate age\cr 8 \tab -
#' \tab - \tab later duplicate age }
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{codeplot}}, \code{\link{zapvelout}}
#' @examples
#'
#' outliers <- velout(age, height, id, heights, limit=3)
#'
#' @export velout
	velout <- function(x, y, id, data, lag=1, velpower=0.5, limit=5, linearise=FALSE)
#	identifies velocity patterns in growth curve outliers
#	returns id, x and code for errors, plus coded velocities relative to neighbours
#	lag (default 1) to calculate velocity
#	velpower is power p in velocity = dy / (dx)^p
#		default 0.5 is halfway between velocity and increment
#	limit (default 5 SD) to identify outliers
#	linearise straightens out growth curves first
{
	mcall <- match.call()
	data <- eval(mcall$data)
#	create data frame of x, y and id with no NAs
	dc <- data[, sapply(mcall[c('x', 'y', 'id')], deparse)]
	nrow <- nrow(data)
	dc <- na.omit(cbind(dc, count=1:nrow))
#	sort into id/x order
	dc <- dc[order(dc[,3], dc[,1]), ]
	if (linearise) {
#	fit spline curve to convert y to residual adjusted for x
# fails if package quantreg not available
    spline.lm <- quantreg::rq(dc[,2] ~ bs(dc[,1], df=5))
		dc[,2] <- residuals(spline.lm)
	}
#	calculate velocity between successive measurements
	dt1 <-  diff(dc[,1], lag=lag)
	vel1 <- diff(dc[,2], lag=lag) / dt1 ^ velpower
	dt2 <-  diff(dc[,1], lag=lag * 2)
	vel2 <- diff(dc[,2], lag=lag * 2) / dt2 ^ velpower
#	flag duplicate ages
	dt1 <- dt1 == 0
#	set missing where ids differ
	idlev <- as.numeric(dc[,3])
	dt1[diff(idlev, lag=lag) != 0] <- FALSE
	vel1[diff(idlev, lag=lag) != 0] <- NA
	vel2[diff(idlev, lag=lag * 2) != 0] <- NA
#	standardise using robust SD
	vel1 <- trunc(vel1 / mad(vel1, na.rm=TRUE)) # or IQR * 1.349
	vel2 <- trunc(vel2 / mad(vel2, na.rm=TRUE))
#	insert NAs to make up length
	vel3 <- c(rep(NA, lag), vel2, rep(NA, lag))
	vel2 <- c(vel1, rep(NA, lag))
	vel1 <- c(rep(NA, lag), vel1)
#	derive outlier code
	code <- (as.numeric(abs(vel1) >= limit) + as.numeric(abs(vel2) >= limit)) * 2 + as.numeric(abs(vel3) >= limit)
	dt2 <- c(dt1, rep(FALSE, lag))
	dt1 <- c(rep(FALSE, lag), dt1)
	code[dt2 | dt1] <- 8
	code[dt2 & !dt1] <- 7
	t <- is.na(vel3) & !(dt1 | dt2)
	code[t] <-
		(as.numeric(!is.na(vel1[t]) & abs(vel1[t]) >= limit) +
		 as.numeric(!is.na(vel2[t]) & abs(vel2[t]) >= limit)) * 6

#  code	v1+v2	v3		interpretation
#	0		0	0		no outlier
#	0		0	NA		no outlier
#	1		0	1		rare pattern
#	2		1	0		complicated - look at curve
#	3		1	1		adjacent to simple outlier
#	4		2	0		single outlier
#	5		2	1		double outlier
#	6		1	NA		edge outlier
#	7		-	-		first duplicate age
#	8		-	-		later duplicate age

	dc <- cbind(dc[, c(3, 1, 2, 4)], code, vel1, vel2, vel3)
	mat <- as.data.frame(matrix(nrow=nrow, ncol=dim(dc)[2],
		dimnames=list(row.names(data), dimnames(dc)[[2]])))
	attr(mat, 'data') <- deparse(mcall$data)
	mat[dc$count, ] <- dc
	if (is.factor(dc[, 1])) {
		mat[, 1] <- as.factor(mat[, 1])
		levels(mat[, 1]) <- levels(dc[, 1])
	}
	mat$count <- NULL
	mat$code <- factor(mat$code)
	cat('code frequencies\n')
	print(summary(mat$code))
	invisible(mat)
}

#############################
#
#	codeplot
#
#############################





#' Plot and zap velocity outliers in growth curves
#'
#' Handles output from \code{velout} function to display growth curves with
#' outlying points, either plotting or zapping the outliers.
#'
#' The function \code{velout} identifies putative outliers for \code{y} in
#' \code{data}, \code{codeplot} plots them, and \code{zapvelout} sets missing
#' those confirmed as outliers.  Codes range from 0 (normal) to 8, where 4 and
#' 6 are conventional outliers (see \code{\link{velout}}).
#'
#' @aliases codeplot zapvelout
#' @param outliers Data frame returned from velout.
#' @param icode The code number(s) defining the subset of curves to be
#' displayed or zapped (between 1 and 6).
#' @param \dots Optional plot parameters.
#' @param print Option to print as well as plot information on each curve.
#' @return \code{codeplot} returns summary information on each curve with an
#' outlier of the relevant code, and optionally plots the curve.
#' \code{zapvelout} sets to NA values of \code{y} whose code is contained in
#' \code{icode}, and returns the modified data frame.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{velout}}
#' @examples
#'
#' ## identify outliers
#' outliers <- velout(age, height, id, heights, limit=2)
#'
#' ## plot outliers with code 4 or 6
#' codeplot(outliers, icode=c(4,6))
#'
#' ## set the 8 outliers missing
#' newheights <- zapvelout(outliers, icode=6)
#'
#' @export codeplot
	codeplot <- function(outliers, icode=4, ..., print=TRUE)
{
#	plots growth curves with outliers identified by velout
#	outliers	data frame returned by velout
#				id, x, y, code, vel1, vel2, vel3
#	icode		type(s) of outlier, from velout
#	...			for plot commands
#	print		TRUE to print case summaries
	cat('click or ESC in plot window to progress through cases', '\nESC in console then ESC in plot window to quit\n')
	id <- names(outliers)[1]
	x <- names(outliers)[2]
	y <- names(outliers)[3]
	outliers$res <- outliers$code %in% icode
	ids <- unique(outliers[outliers$res, 1]) # ids with outliers
	nids <- length(ids)
	dat <- outliers[outliers[,1] %in% ids,]
	dat <- within(dat, {
		vel1[is.na(vel1)] <- 999
		vel2[is.na(vel2)] <- 999
		vel3[is.na(vel3)] <- 999
	} )
	dat <- na.omit(dat)
	dat <- within(dat, {
		vel1[vel1 == 999] <- NA
		vel2[vel2 == 999] <- NA
		vel3[vel3 == 999] <- NA
	} )
	el <- new.env()
	assign('kount', 0, envir=el)
	by(dat, factor(dat[, 1]), function(z) {
		kount <- get('kount', envir=el) + 1
		assign('kount', kount, envir=el)
		if (print) cat('\ncase', kount, 'of', nids, '-', id, as.character(z[1, 1]), 'with', nrow(z), 'points\n')
		inc <- z[z$res, ]
		if (print) print(inc[, 1:7])
		ox <- order(z[, 2])
		dots <- list(...)
		if (is.null(dots$type)) dots$type <- 'b'
		if (is.null(dots$xlab)) dots$xlab <- x
		if (is.null(dots$ylab)) dots$ylab <- y
		pcht <- dots$pch
		dots$pch <- NULL
		do.call('plot', c(list(z[, 2][ox], z[, 3][ox], pch=46), dots))
		dots$type <- NULL
		if (is.null(pcht)) pcht <- 8
		do.call('points', c(list(inc[, 2], inc[, 3], pch=pcht), dots))
		title(paste(id, z[1, 1], 'code',
			paste(unique(inc$code), collapse=' ')))
		locator(1)
	} )
	invisible()
}

#############################
#
#	zapvelout
#
#############################

#' @rdname codeplot
#' @export
	zapvelout <- function(outliers, icode)
{
#	set outliers missing for specified code(s)
	data <- get(attr(outliers, 'data'))
	if (nrow(data) != nrow(outliers)) stop('numbers of rows in data and outliers differ')
	y <- names(outliers)[3]
	zap <- outliers$code %in% icode
	cat(sum(zap), y, 'values set missing\n')
	data[outliers$code %in% icode, y] <- NA
	invisible(data)
}

################################################################
#
#	recalib to recalibrate x,y data using SITAR random effects
#
################################################################





#' Recalibrate x, y data using SITAR random effects
#'
#' A function to recalibrate x,y data using SITAR random effects
#'
#' \code{recalib} recalibrates the values of \code{xc} and \code{yc} based on
#' \code{model}. \code{xc} values are changed to:
#'
#' (xc-c(coef[from,'b']))*exp(coef[from,'c']-coef[to,'c'])+coef[to,'b'].
#'
#' \code{yc} values are changed to: \code{yc-coef[from,'a']+coef[to,'a']}.
#'
#' @param xc character vector defining column name(s) of \code{x} data to be
#' recalibrated.
#' @param yc character vector defining column name(s) of \code{y} data to be
#' recalibrated.
#' @param id factor defining \code{from} and \code{to} rows. If \code{NULL}
#' then recalibrate all rows.
#' @param data dataframe containing \code{xc}, \code{yc} and \code{id}.
#' @param xcnew column names for replacement columns \code{xc}. If default
#' \code{NULL} then use names xcnew1... .
#' @param ycnew column names for replacement columns \code{yc}. If default
#' \code{NULL} then use names ycnew1... .
#' @param model \code{sitar} model defining the random effects to be used for
#' recalibration.
#' @param from level of \code{id} defining existing data (must be a single row
#' in \code{coef{model}}).
#' @param to level of \code{id} defining data to be recalibrated (a single row
#' in \code{coef{model}}).
#' @return Returns the dataframe \code{data} with the \code{from} rows of
#' \code{xc} and \code{yc} recalibrated.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @export recalib
	recalib <- function(xc, yc, id=NULL, data, xcnew=NULL, ycnew=NULL, model, from, to)
#	xc - column names of x data to be recalibrated
#	yc - column names of y data to be recalibrated
#	id - column of data containing to/from
#	(default NULL means include all, else just 'from' rows)
#	data - dataframe containing x and y data to be recalibrated
#	xcnew - column names for replacement x (default NULL=xnew1...)
#	ycnew - column names for replacement y (default NULL=ynew1...)
#	model - SITAR model defining random effects
#	from - id value of existing data (must be row in random effects)
#	to - id value to recalibrate to (ditto)

#	returns dataframe data, with x and y values changed
{
	xcall <- model$call.sitar$x
	ycall <- model$call.sitar$y
	mux <- model$mux
	coef <- random.effects(model)
	dcoef <- coef[from,] - coef[to,]
	include <- with(data, {
		if (is.null(id)) rep(TRUE, dim(data)[[1]]) else as.character(id) %in% from
	} )
	nxc <- length(xc)
	if (is.null(xcnew)) xcnew <- paste('xnew', 1:nxc, sep='')
	for (i in 1:nxc) {
		data[, xcnew[i]] <- data[, xc[i]]
		x <- c(data[include, xcnew[i]])
		x <-funcall(x, xcall)
		x <- (x - mux - c(coef[from, 'b'])) * c(exp(dcoef$c)) + c(coef[to, 'b']) + mux
		# print(signif(x, 2))
		data[include, xcnew[i]] <- funcall(x, xcall, TRUE)
		# print(signif(funcall(x, xcall, TRUE), 2))
	}
	nyc <- length(yc)
	if (is.null(ycnew)) ycnew <- paste('ynew', 1:nyc, sep='')
	for (i in 1:nyc) {
		data[, ycnew[i]] <- data[, yc[i]]
		y <- c(data[include, ycnew[i]])
		y <-funcall(y, ycall)
		y <- y - dcoef$a
		data[include, ycnew[i]] <- funcall(y, ycall, TRUE)
	}
	invisible(data)
}

###################
#	xaxsd  yaxsd  #
###################





#' Par args xaxs and yaxs option d
#'
#' Implements par('xaxs') and par('yaxs') option 'd'.
#'
#' Implements par('xaxs') and par('yaxs') option 'd', i.e. uses previous axis
#' scales in a new plot.
#'
#' @aliases xaxsd yaxsd
#' @param usr a length-2 vector defining the length of the x-axis or y-axis.
#' @return By default returns xlim/ylim args to match current setting of
#' par()$usr, i.e. previous plot scales.  Specifying \code{usr} gives scales
#' with the usr args at the extremes. If par('xlog') or par('ylog') are set the
#' returned limits are antilogged (to base 10).
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#'
#' ## generate and plot 100 data points
#' x <- rnorm(100)
#' y <- rnorm(100)
#' plot(x, y, pch=19)
#'
#' ## generate and plot 10 more
#' ## constraining axis scales to be as before
#' x <- rnorm(10)
#' y <- rnorm(10)
#' plot(x, y, pch=19, xlim=xaxsd(), ylim=yaxsd())
#'
#' ## force axis extremes to be -3 and 3
#' plot(x, y, pch=19, xlim=xaxsd(c(-3,3)), ylim=yaxsd(c(-3,3)))
#'
#' @export xaxsd
xaxsd <- function(usr=par()$usr[1:2]) {
  if (!missing(usr) && par('xlog')) usr <- log10(usr)
  usr <- (usr + mean(usr) * 0.08) / 1.08
	if (par('xlog')) 10 ^ usr else usr
}
#' @rdname xaxsd
#' @export yaxsd
yaxsd <- function(usr=par()$usr[3:4]) {
  if (!missing(usr) && par('ylog')) usr <- log10(usr)
	usr <- (usr + mean(usr) * 0.08) / 1.08
	if (par('ylog')) 10 ^ usr else usr
}

#	implements xaxs/yaxs option 'd'
#	by default returns xlim/ylim args to match current setting of par()$usr, i.e. previous plot scales
#	e.g. plot(..., xlim=xaxsd(), ylim=yaxsd())
#	specifying usr gives scale with usr args at extremes
# fails with plot.sitar: par([xy]log=TRUE) needed before plot(*, log='[xy]')

##########################
#	sampling for subset  #
##########################





#' Sample from SITAR dataset
#'
#' A function to sample from a SITAR dataset for experimental design purposes.
#' Two different sampling schemes are offered, based on the values of \code{id}
#' and \code{x}.
#'
#' With the first sampling scheme \code{xlim} is set to \code{NULL} (default),
#' and rows of \code{data} are sampled with probability \code{prob} without
#' replacement.  With the second sampling scheme \code{xlim} is set to a range
#' within \code{range(x)}.  Subjects \code{id} are then sampled with
#' probability \code{prob} without replacement, and all their rows where
#' \code{x} is within \code{xlim} are selected.  The second scheme is useful
#' for testing the power of the model to predict later growth when data only up
#' to a certain age are available. Setting \code{xlim} to \code{range(x)}
#' allows data to be sampled by subject. The returned value can be used as the
#' \code{subset} argument in \code{sitar} or \code{update.sitar}.
#'
#' @param x vector of age.
#' @param id factor of subject identifiers.
#' @param data dataframe containing \code{x} and \code{id}.
#' @param prob scalar defining sampling probability. See Details.
#' @param xlim length 2 vector defining range of \code{x} to be selected. See
#' Details.
#' @return Returns a logical the length of \code{x} where \code{TRUE} indicates
#' a sampled value.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{sitar}}
#' @examples
#'
#' ## draw 50% random sample
#' s50 <- subsample(age, id, heights, prob=0.5)
#'
#' ## truncate age range to 7-12 for 50% of subjects
#' t50 <- subsample(age, id, heights, prob=0.5, xlim=c(7, 12))
#'
#' @export subsample
	subsample <- function(x, id, data, prob = 1, xlim = NULL)
{
#	selects subset for sitar model design
#	x - age
#	id - subject id
#	data - data frame containing x and id
#	prob - probability of being selected
#		if xlim is NULL all points are sampled
#		otherwise IDs are sampled and their ages limited to xlim
#	xlim - age range to be selected (default NULL)
#	returns subset as logical vector of rows to be included
	mcall <- match.call()
	data <- eval(mcall$data)
	x <- eval(mcall$x, data)
	nx <- length(x)
	subset <- rep(TRUE, nx)
#	simple sampling of whole dataset
	if (is.null(xlim)) {
		ss <- sample(nx, prob * nx)
		subset[-ss] <- FALSE
	}
	else {
#	sampling of individuals and restricting their age range to xlim
		id <- eval(mcall$id, data)
		nid <- nlevels(factor(id))
		lid <- levels(factor(id))
		sid <- sample(lid, prob * nid)
		subset[id %in% sid & (x < min(xlim) | x > max(xlim))] <- FALSE
	}
	subset
}

################
#	wishlist   #
################

#	plot.sitar
#		add acceleration curve option
#		improve ylim for dau plots
