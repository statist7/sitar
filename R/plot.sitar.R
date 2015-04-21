	plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun, yfun, subset=NULL, abc=c(a=0, b=0, c=0), add=FALSE, nlme=FALSE, ...)
#	plot curves from sitar model
#	opt:
#		d = fitted distance curve (labels[1] = x, labels[2] = y)
#		v = fitted velocity curve (labels[3] = y)
#		e = fitted fixed effects distance curve (labels[1] = x, labels[2] = y)
#		f = fitted distance curves by subject
#		u = unadjusted y vs t curves by subject
#		a = adjusted y vs adjusted t curves by subject
#
#		multiple options all plot on same graph
#			(though axes not necessarily optimised)
#
#		labels are for opt dv - particularly v
#		or use xlab and ylab and y2par
#
#		apv TRUE draws vertical line at age of peak velocity
#		and returns apv/pv (respecting xfun/yfun settings)
#
#		xfun/yfun are functions to apply to x/y before plotting
#
#		subset is subset of values
#
#		abc is a set of named sitar parameters for opt dv e.g.
#			abc=list(a=1, b=0.1, c=-0.1)
#		or a single id level whose abc values are to be used
#
#		add TRUE overwrites previous graph (or use lines)
#
#		nlme TRUE plots model as nlme object
{
  xseq <- function(x, nx=101) {
    rx <- range(x, na.rm=TRUE)
    seq(rx[1], rx[2], length.out=nx)
  }

	model <- x
	rm(x)
	if (nlme) {
		mcall <- match.call()[-1]
		names(mcall)[[1]] <- 'x'
		do.call('plot.lme', as.list(mcall))
	}
	else {
		mcall <- model$call.sitar
		data <- eval(mcall$data, parent.frame())
#	subset used to fit model
		subset <- eval(mcall$subset, data)
		if (!is.null(subset)) data <- data[subset, ]
		x <- eval(mcall$x, data)
		y <- eval(mcall$y, data)
		id <- factor(eval(mcall$id, data)) # factor added 23/4/13
		nf <- length(fitted(model))
		if (nf != length(y)) stop(paste('model (length=', nf, ') incompatible with data (rows=', length(y), ')', sep=''))

#	extract list(...)
		ccall <- match.call()[-1]
#	subset to plot model
		subset <- eval(ccall$subset, data, parent.frame())
		if (is.null(subset)) subset <- rep(TRUE, nf)
		cnames <- names(ccall)
		dots <- cnames[!cnames %in% names(formals(plot.sitar))]
		if (length(dots) > 0) ARG <- lapply(ccall[dots], eval, data, parent.frame())
			else ARG <- NULL
#	if xlab not specified replace with label or else x name
		if (!"xlab" %in% names(ARG)) {
			xl <- ifelse(missing(labels), deparse(mcall$x), labels[1])
			ARG <- c(ARG, list(xlab=xl))
		}
		else xl <- ARG$xlab
#	if ylab not specified replace with label or else y name
		if (!"ylab" %in% names(ARG)) {
			yl <- ifelse(missing(labels), deparse(mcall$y), labels[2])
			ARG <- c(ARG, list(ylab=yl))
		}
		else yl <- ARG$ylab

#	create output list
		xy <- list()

#	plot fitted distance and velocity curves
		if (grepl("d", opt) || grepl("v", opt) || apv) {
			xt <- xseq(x[subset])
  		newdata <- data.frame(x=xt)
# if subset, estimate mean values for covariates
      if (!identical(subset, rep(TRUE, nf))) {
  			argnames <- names(formals(model$fitnlme))
  			xtra <- argnames[!argnames %in% c('x', names(fixef(model)))]
        if (length(xtra) > 0) {
          df <- update(model, returndata=TRUE)[subset, xtra]
          xtra <- unlist(lapply(df, mean, na.rm=TRUE))
     			newdata <- data.frame(newdata, t(xtra))
        }
     }
			yt <- predict(object=model, newdata=newdata, level=0)
#	adjust for abc
			if (!missing(abc)) {
#	if abc is named ensure 3 values match model random effects
				if (!is.null(names(abc))) {
					random <- dimnames(ranef(model))[[2]]
					for (l in letters[1:3])
						if (is.na(abc[l]) || !l %in% random) abc[l] <- 0
				}
				else
#	else abc is id level
				if (length(abc) == 1) {
					idabc <- dimnames(ranef(model))[[1]] %in% abc
					if (sum(idabc) == 0) stop(paste('id', abc, 'not found'))
					abc <- ranef(model)[idabc,]
				}
				else stop('abc should be either single id level or up to three named random effect values')
				xoffset <- model$xoffset
				if (is.null(xoffset)) xoffset <- 0
				if (!is.na(fixef(model)['b'])) xoffset <- xoffset + fixef(model)['b']
				xt <- (xt - xoffset) / exp(abc[['c']]) + xoffset + abc[['b']]
				yt <- yt + abc[['a']]
			}
#	derive cubic smoothing spline curve
			xy$ss <- ss <- makess(xt, yt, xfun=xfun, yfun=yfun)
			ss1 <- predict(ss, deriv=1)
			ss2 <- predict(ss, deriv=2)
			if (missing(labels)) labels <- c(xl, yl, paste(yl, 'velocity'))
			# if (missing(labels)) labels <- c(xl, yl, ifelse(typeof(yl) == 'expression', expression(paste(as.character(yl), '~~velocity')), paste(yl, 'velocity')))
			if (grepl("d", opt) && grepl("v", opt)) {
				xy <- do.call("y2plot", c(list(x=ss$x, y1=ss$y, y2=ss1$y, labels=labels, add=add, xy=xy), ARG))
				add <- TRUE
			} else
			if (grepl("d", opt)) {
				xy <- do.call("y2plot", c(list(x=ss$x, y1=ss$y, add=add, xy=xy), ARG))
				add <- TRUE
			} else
			if (grepl("v", opt)) {
				ARG$ylab <- labels[3]
				xy <- do.call("y2plot", c(list(x=ss$x, y1=ss1$y, add=add, xy=xy), ARG))
				add <- TRUE
			}
		}

#	plot fixed effects distance curve
		if (grepl("e", opt)) {
  		xt <- xseq(x[subset])
      yt <- predict(model$ns, newdata=data.frame(x=xt))
			if (!missing(xfun)) xt <- xfun(xt)
			if (!missing(yfun)) yt <- yfun(yt)
			ox <- order(xt)
			xy <- do.call("y2plot", c(list(x=xt[ox], y1=yt[ox], add=add, xy=xy), ARG))
			add <- TRUE
		}

#	plot y vs t by subject
		if (grepl("u", opt)) {
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
  		do.call("mplot", c(list(x=x, y=y, id=id, subset=subset, add=add), ARG))
			add <- TRUE
		}

#	plot fitted curves by subject
		if (grepl("f", opt)) {
			y <- fitted(model, level=1)
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
    	do.call("mplot", c(list(x=x, y=y, id=id, subset=subset, add=add), ARG))
			add <- TRUE
		}

#	plot adjusted y vs adjusted t by subject
		if (grepl("a", opt)) {
			fred <- summary(model)
			x <- fred$x.adj
			y <- fred$y.adj
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
    	do.call("mplot", c(list(x=x, y=y, id=id, subset=subset, add=add), ARG))
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
