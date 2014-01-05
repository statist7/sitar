#	global variables
	if (getRversion() >= "3.0.0") utils::globalVariables('.par.usr2')

#############################
#
#	print.sitar
#
#############################

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

	summary.sitar <- function (object, adjustSigma = TRUE, verbose = FALSE, ...) 
{
    class(object) <- class(object)[-1]
    object <- summary(object, adjustSigma=adjustSigma, verbose=verbose, ...)
    
#	adjust x and y for random effects a, b and c
	random <- as.character(object$call$random)[2]
	mcall <- object$call.sitar
	data <- eval(mcall$data)
	if (!is.null(mcall$subset)) data <- data[eval(mcall$subset),]
	x <- eval(mcall$x, data)
	y <- eval(mcall$y, data)
	id <- eval(mcall$id, data)
	nf <- length(fitted(object))
	if (nf != length(y)) stop(paste('model (length=', nf, ') incompatible with data (rows=', length(y), ')', sep=''))
	xoffset <- object$xoffset
	if (is.null(xoffset)) xoffset <- 0
	if (!is.na(fixef(object)['b'])) xoffset <- xoffset + fixef(object)['b'] # added 23/4/13
	# if (is.null(xoffset)) xoffset <- fixef(object)['b'] # was mean(x) # commented out 23/4/13
	# if (is.na(xoffset)) xoffset <- 0 # b not fixed effect # commented out 23/4/13
	x.adj <- x - xoffset
	if (grepl('b', random)) x.adj <- x.adj - ranef(object)[factor(id),'b'] # changed from coef 1/8/12
	if (grepl('c', random)) x.adj <- x.adj * exp(ranef(object)[factor(id),'c'])
	object$x.adj <- x.adj + xoffset # + bfe removed 1/8/12
	y.adj <- y
	if (grepl('a', random)) y.adj <- y.adj - ranef(object)[factor(id),'a']
	object$y.adj <- y.adj

#	save age at peak velocity
	object$apv <- makess(x, fitted(object, level=0))$apv
	
    class(object) <- c("summary.sitar", class(object))
    object
}

#############################
#
#	print.summary.sitar
#
#############################

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
#	plot.sitar
#
#############################

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
#		abc is a set of labelled sitar parameters for opt dv e.g.
#			abc=list(a=1, b=0.1, c=-0.1)
#		or a single id level whose abc values are to be used
#
#		add TRUE overwrites previous graph (or use lines)
#
#		nlme TRUE plots model as nlme object
{
	model <- x
	rm(x)
	if (nlme) {
		mcall <- match.call()[-1]
		names(mcall)[[1]] <- 'x'
		do.call('plot.lme', as.list(mcall))
	} 
	else {
		mcall <- model$call.sitar
		data <- eval(mcall$data)
#	subset used to fit model
		if (!is.null(mcall$subset)) data <- data[eval(mcall$subset),]
		x <- eval(mcall$x, data)
		y <- eval(mcall$y, data)
		id <- factor(eval(mcall$id, data)) # factor added 23/4/13
		nf <- length(fitted(model))
		if (nf != length(y)) stop(paste('model (length=', nf, ') incompatible with data (rows=', length(y), ')', sep=''))

#	extract list(...)
		ccall <- match.call()[-1]
#	subset to plot model
		subset <- eval(ccall$subset, data)
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
			xt <- x[subset]
			yt <- fitted(model, level=0)[subset]
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
				xt <- (xt - xoffset + abc[['b']]) / exp(abc[['c']]) + xoffset
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
				xy <- do.call("y2plot", c(list(x=quote(ss$x), y1=quote(ss$y), add=add, xy=xy), ARG))
				add <- TRUE
			} else
			if (grepl("v", opt)) {
				ARG$ylab <- labels[3]
				xy <- do.call("y2plot", c(list(x=quote(ss$x), y1=quote(ss1$y), add=add, xy=xy), ARG))
				add <- TRUE
			}
		}
		
#	plot fixed effects distance curve
		if (grepl("e", opt)) {
			xt <- x[subset]
			yt <- model$ns$fitted[subset]
			if (!missing(xfun)) xt <- xfun(xt)
			if (!missing(yfun)) yt <- yfun(yt)
			ox <- order(xt)
			xy <- do.call("y2plot", c(list(x=quote(xt[ox]), y1=quote(yt[ox]), add=add, xy=xy), ARG))
			add <- TRUE
		}

#	plot y vs t by subject
		if (grepl("u", opt)) {
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
			do.call("mplot", c(list(x=quote(x), y=quote(y), id=quote(id), subset=quote(subset), add=add), ARG))
			add <- TRUE
		}
	
#	plot fitted curves by subject
		if (grepl("f", opt)) {
			y <- fitted(model, level=1)
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
			do.call("mplot", c(list(x=quote(x), y=quote(y), id=quote(id), subset=quote(subset), add=add), ARG))
			add <- TRUE
		}
	
#	plot adjusted y vs adjusted t by subject
		if (grepl("a", opt)) {
			fred <- summary(model)
			x <- fred$x.adj
			y <- fred$y.adj
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
			do.call("mplot", c(list(x=quote(x), y=quote(y), id=quote(id), subset=quote(subset), add=add), ARG))
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

	lines.sitar <- function (x, ...) 
{
	mcall <- match.call()[-1]
	if (!"add" %in% names(mcall)) mcall <- c(as.list(mcall), list(add=TRUE))
	do.call('plot.sitar', as.list(mcall))
}
		
#############################
#
#	update.sitar
#
#############################

	update.sitar <- function (object, ..., evaluate = TRUE) 
{
	mcall <- object$call.sitar
	if (is.null(mcall)) 
		stop("need an object with call.sitar component")
	extras <- as.list(match.call(expand.dots = FALSE)$...)
#	drop null arguments
	for (i in names(extras)) 
		if (is.null(extras[[i]])) 
			mcall[[i]] <- extras[[i]] <- NULL
	start. <- list(fixed=fixef(object), random=ranef(object))
#	check if xoffset used
	bstart <- start.$fixed['b']
	if (is.na(bstart)) bstart <- 0 else names(bstart) <- NULL
	if (is.numeric(object$xoffset)) {
		mcall$xoffset <- NULL
		bstart <- bstart + object$xoffset
		# bstart <- object$xoffset # alternative
		if (is.null(extras$bstart)) extras$bstart <- bstart
	}		
	if (length(extras) > 0) {
#	update formulae
		all.pars <- c(as.list(mcall), formals(sitar))[-1]
		for (a in letters[1:3]) {
			pos.e <- pmatch(paste(a, 'form', sep='.'), names(extras), 0)
			if (pos.e) {
				pos.p <- pmatch(paste(a, 'formula', sep='.'), names(all.pars), 0)
				if (pos.p) extras[[pos.e]] <- update.formula(all.pars[[pos.p]], extras[[pos.e]])
			}
		}
#	update existing arguments
		existing <- pmatch(names(extras), names(mcall))
		if (sum(existing, na.rm=TRUE))
			for (a in 1:length(existing)) 
				mcall[existing[a]] <- extras[a]
#	add new arguments
		existing <- !is.na(existing)
		if (any(!existing))
			mcall <- as.call(c(as.list(mcall), extras[!existing]))
	}
#	check if can use previous start values
	if (sum(pmatch(names(extras), c("x", "y", "id", "data", "fixed", "random", "a.formula", "b.formula", "c.formula", "start", "subset", "returndata")), na.rm=TRUE) == 0) {
#	update start random effects if dataframe changed
		data <- eval(mcall$data)
		id <- factor(eval(mcall$id, data))
		if (nlevels(id) != dim(ranef(object))[1]) {
#	omit random effects for missing levels in id
			idcheck <- rownames(ranef(object)) %in% levels(id)
			start.$random <- ranef(object)[idcheck,]
			cat(dim(ranef(object))[1] - sum(idcheck), 'subjects omitted\n')
#	add zero random effects for new levels in id
			newid <-!levels(id) %in% rownames(ranef(object))
			if (sum(newid) > 0) {
				newre <- matrix(0, nrow=sum(newid), ncol=dim(ranef(object))[2], dimnames=list(levels(id)[newid], dimnames(ranef(object))[[2]]))
				start.$random <- rbind(start.$random, newre)
				cat(sum(newid), 'subjects added\n')
			}
		}	
#	update start fixed effects if df, knots, bounds or bstart updated
		if (sum(pmatch(names(extras), c("df", "knots", "bounds", "bstart")), na.rm=TRUE) > 0) {
			x <- eval(mcall$x, data)
			if (!is.null(object$bstart)) bstart <- object$bstart
			knots <- attr(object$ns$model$ns, 'knots') + bstart
			bounds <- attr(object$ns$model$ns, 'Boundary.knots') + bstart
			df <- object$ns$rank - 1
			if (length(fixef(object)) > df + 1) fixed.extra <- (df+2):length(fixef(object))
				else fixed.extra <- NULL
			if (!is.null(extras$knots)) {
				knots <- extras$knots
				df <- length(knots) + 1
			}
			else if (!is.null(extras$df)) {
				df <- extras$df
				knots <- quantile(x, (1:(df-1))/df)
			}
			if (!is.null(extras$bounds)) {
				bounds <- extras$bounds
				if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x))
			}
			if (!is.null(extras$bstart)) {
				bstart <- eval(extras$bstart)
				if (is.character(bstart)) bstart <- mean(x)
				start.$fixed['b'] <- bstart
			}
			start.$fixed['c'] <- 0
			knots <- knots - bstart
			bounds <- bounds - bstart
#	get spline start values
			spline.lm <- lm(fitted(object, level=0) ~ ns(x - bstart, knots=knots, Bound=bounds))
			start.$fixed <- c(coef(spline.lm)[c(2:(df+1), 1)], start.$fixed[fixed.extra])
		}
#	save start. object
		assign('start.', start., parent.frame())
		if (!'start' %in% names(mcall)) 
			mcall <- as.call(c(as.list(mcall), start=quote(start.)))
	}
	else mcall$start <- NULL
	if (evaluate) 
		eval(mcall, parent.frame())
	else mcall
}

#############################
#
#	y2plot
#
#############################

	y2plot <- function(x, y1, y2=NULL, labels, y2par=NULL, add=FALSE, xy=NULL, xlegend="topleft", inset=0.04, ...)
#	plot x versus y1 and y2 using left/right axes for the two y's
#	returns par() for y1 or list(par(), par("usr")) with y2
#	labels are labels for x, y1 and y2 (xlab and ylab override x and y1, y2par['ylab'] overrides y2)
#	y2par is optional list of par args for y2 axis
#	add suppresses plot axes
#	xy set uses par() from previous call
#	xlegend and inset place legend (NULL suppresses)
{
# get axis labels
	if (missing(labels)) labels <- c(deparse(substitute(x)), deparse(substitute(y1)), deparse(substitute(y2)))
#	plot y1
	ypar <- list(...)
#	save y1 par args
	opar <- par(no.readonly=TRUE)
	lty <- 1:2
	lwd <- rep(par('lwd'), 2)
	col <- rep(par('col'), 2)
	for (i in c('lty', 'lwd', 'col')) {
		j <- get(i)
		if (i %in% names(ypar)) {
			j[1] <- ypar[i]
			assign(i, unlist(j))
		}
		else ypar[i] <- j[1]
	}
#	save y2 par args
	if (!missing(y2)) {
		for (i in c('lty', 'lwd', 'col')) {
			j <- get(i)
			if (i %in% names(y2par)) {
				j[2] <- y2par[i]
				assign(i, unlist(j))
			}
			else y2par[i] <- j[2]
		}
	}
#	if a new plot draw axes
	if (!add) {
#	if mar specified then set it
		if (!is.null(ypar$mar)) {
			mar <- ypar$mar
		}
		else {
#	otherwise reset mar
			mar <- c(5,4,4,2) + 0.1
#	if y2 set increase mar[4]
			if (!missing(y2)) mar[4] <- 4.1
		}
		par(mar=mar)
#	plot x & y1 axes and y1 ~ x
		if (is.null(ypar$xlab)) ypar$xlab <- labels[1]
		if (is.null(ypar$ylab)) ypar$ylab <- labels[2]
		do.call('plot', c(list(x=quote(x), y=quote(y1), type='l'), ypar))
#	save x & y1 axis limits
		xy$usr <- par('usr')
#	optionally plot y2 and add right axis, with y2par as ...
		if (!missing(y2)) {
			par(new=TRUE)
#	ensure x axis same for y2 as y1
			if (!is.null(ypar$xlim) && is.null(y2par$xlim)) y2par$xlim <- ypar$xlim
			if (is.null(y2par$ylab)) y2par$ylab <- labels[3]
			do.call('plot', c(list(x=quote(x), y=quote(y2), ann=FALSE, bty="n", xaxt="n", yaxt="n", type="l"), y2par))
#	save y2 axis limits
			xy$usr2 <- par('usr')
			eval(parse(text=".par.usr2 <<- par('usr')"))
			# .par.usr2 <<- par('usr')
#	add y2 axis 
			if (par('mar')[4] >= 2) axis(4)
#	unset col
			y2par$col <- NULL
#	add y2 axis label
			if (par('mar')[4] >= 3) do.call('mtext', c(list(text=labels[3], side=4, line=3, cex=par()$cex), y2par))
#	add legend
			if (!is.null(xlegend) && !is.null(inset)) legend(xlegend, legend=labels[2:3], bty="o", lty=lty, lwd=lwd, col=col, inset=inset)
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
				# .par.usr2 <<- par('usr')
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
#	plotclean
#
#############################

	plotclean <- function(x, y=NULL, id=NULL, data=parent.frame(), n=length(x), pch=20, ...)
#	plot growth curves to identify outlying points and curves
#	plot y ~ x by id with data
#	identify up to n points, with pch as plot character
{	if (!is.null(data)) {
		if (!deparse(substitute(data)) %in% search()) {
			on.exit(detach(data))
			attach(data)	}	}
	xt <- substitute(x); xlab <- deparse(xt)
	yt <- substitute(y); ylab <- deparse(yt)
	if (is.null(id)) {
		plot(x, y, xlab=xlab, ylab=ylab, pch=46, ...)
		idlab <- NULL
	}
	else {
#	mplot(x, y, id, data, col="gray", type="o", pch=46)
		idlab <- deparse(substitute(id))
		id <- factor(id)
		plot(x, y, type="n", ...)
		tt <- by(data, id, function (z) lines(eval(yt, z) ~ eval(xt, z), col="gray", type="o", pch=46, ...))
	}
	xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
	sel <- rep(FALSE, length(x)); res <- integer(0)
	title('click on outliers in plot - then right-click to escape')
	while (sum(sel) < n) {
		ans <- identify(x[!sel], y[!sel], n=1, plot=FALSE, ...)
		if (!length(ans)) break
		ans <- which(!sel)[ans]
		points(x[ans], y[ans], pch = pch)
		if (!is.null(id)) {
			text(x[ans], y[ans], label=id[ans], adj=c(0.5, 0.1))
			idt <- id==id[ans]
			ox <- order(x[idt])
			lines(x[idt][ox], y[idt][ox])
		}	
		sel[ans] <- TRUE
		res <- c(res, ans)
	}
	res <- res[order(res)]
	if (is.null(data)) res else list(rows=res, data=data[res, c(idlab, xlab, ylab)])
}
	
#############################
#
#	do.call...
#
#############################

	do.call... <- function(what, args, ..., quote=FALSE, envir=parent.frame()) {
#	do.call avoiding clashes between what args and ... args
#	... args override what args
#	what	function name
#	args	list of args for what
#	quote and envir	from do.call
#	...		args from parent call

	ddd <- list(...)
	if (!is.null(ddd)) {
		for (i in names(ddd)) {
			if (i %in% names(args)) args[i] <- ddd[i]
			else {
				args <- c(args, i=ddd[i])
				names(args)[length(args)] <- i
			}
		}
	}
	do.call(what, args, quote=quote, envir=envir)
}
		
#############################
#
#	varexp
#
#############################

	varexp <- function(..., pattern=NULL)
#	returns % of variance explained by sitar model(s)
{	ARG <- list(...)
	if (!is.null(pattern)) {
		pattern <- ls(envir=parent.frame(), pattern=pattern)
		ARG <- c(ARG, lapply(as.list(pattern), get))
	}
	pc <- lapply(ARG, function(obj) {
		if (!'ns' %in% names(obj)) NA else
		100 * (1 - (obj$sigma / summary(obj$ns)$sigma)^2)		
	})
	names(pc) <- c(match.call(expand.dots=FALSE)$..., pattern)
	pc <- unlist(pc[!is.na(pc)])
	if (length(pc) == 0) return(invisible())
	round(pc[rev(order(pc))], 2)
}

#############################
#
#	funcall
#
#############################

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
#	makess
#
#############################

	makess <- function(x, y, xfun, yfun)
#	summarises y ~ x as smooth.spline
#	returns smooth.spline including apv/pv (peak velocity)
#		or amv/mv (max velocity)
{
	if(!missing(xfun)) x <- xfun(x)		
	if(!missing(yfun)) y <- yfun(y)		
	o <- order(x)
	af <- unique(cbind(x[o], y[o]))
	ss <- smooth.spline(af)
	ss1 <- predict(ss, deriv=1) # 1st derivative
	ss2 <- predict(ss, deriv=2) # 2nd derivative
	apv <- vector(length=2)
	ts <- as.ts(ss2$y) # 2nd derivative as time series
	ts <- ts * lag(ts) < 0 # 2nd derivative changes sign
	if (any(ts, na.rm=TRUE)) { # turning point(s) found
		ts <- c(NA, ts) # expand
		iapv <- which.max(ss1$y * ts) # point of peak velocity
		o <- (iapv-2):(iapv+2) # region of peak
		mn <- mean(ss$x[o])
		mc <- lm(ss1$y[o] ~ poly(ss$x[o] - mn, 2, raw=TRUE))$coef # velocity as quadratic in age
		apv[1] <- - mc[2] / 2 / mc[3] + mn # age at peak velocity
		apv[2] <- mc[1] - mc[2]^2 / 4 / mc[3] # peak velocity
		names(apv) <- c('apv', 'pv')
		if (any(is.na(apv))) warning('apv undefined')
	}
	else {
		iapv <- which.max(ss1$y) # max not peak velocity
		apv[1] <- ss$x[iapv] # age at max velocity
		apv[2] <- ss1$y[iapv] # max velocity
		names(apv) <- c('amv', 'mv')
	}
	ss$apv <- apv
	invisible(ss)
}

#############################
#
#	bupdate
#
#############################

#	update value of bstart to minimise b-c correlation
	bupdate <- function(x) {
	cov <- as.matrix(x$modelStruct[[1]][[1]])
	fixef(x)['b'] + cov[2,3] / cov[3,3]
	}

#############################
#
#	velout
#
#############################

	velout <- function(x, y, id, data, lag=1, velpower=0.5, limit=5, linearise=FALSE)
#	identifies velocity patterns in growth curve outliers
#	returns id, x and code for errors, plus coded velocities relative to neighbours
#	lag (default 1) to calculate velocity
#	velpower is power p in velocity = dy / (dx)^p
#		default 0.5 is halfway between velocity and increment
#	limit (default 5 SD) to identify outliers
#	linearise straightens out growth curves first
{	
	if (missing(x) | missing(y) | missing(id) | missing(data)) stop('x, y, id and data must all be specified')
	data.name <- deparse(substitute(data))
	if (data.name %in% search()) stop('data already attached - needs detaching')
#	create data frame of x, y and id with no NAs
	nrow <- dim(data)[1]
	xyid <- c(deparse(substitute(x)), deparse(substitute(y)), deparse(substitute(id)))
	dc <- data[, xyid]
	dc <- na.omit(cbind(dc, count=1:nrow))
#	sort into id/x order
	dc <- dc[order(dc[,3], dc[,1]), ]
	if (linearise) {
#	fit spline curve to convert y to residual adjusted for x
		require(quantreg)
		spline.lm <- rq(dc[,2] ~ bs(dc[,1], df=5))
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

	dc <- cbind(dc[,c(3,1,2,4)], code, vel1, vel2, vel3)
	mat <- as.data.frame(matrix(nrow=nrow, ncol=dim(dc)[2],
		dimnames=list(row.names(data), dimnames(dc)[[2]])))
	attr(mat, 'data') <- data.name
	mat[dc$count,] <- dc
	if (is.factor(dc[,1])) {
		mat[,1] <- as.factor(mat[,1])
		levels(mat[,1]) <- levels(dc[,1])
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
	by(dat, factor(dat[,1]), function(z) {
		kount <- get('kount', envir=el) + 1
		assign('kount', kount, envir=el)
		if (print) cat('\ncase', kount, 'of', nids, '-', id, as.character(z[1, 1]), 'with', dim(z)[[1]], 'points\n')
		inc <- z[z$res,]
		if (print) print(inc[,1:7])
		ox <- order(z[,2])
		do.call...('plot', list(z[,2][ox], z[,3][ox], type='b', xlab=x, ylab=y, pch=46), ...)
		do.call...('points', list(inc[,2], inc[,3], pch=8), ...)
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

	zapvelout <- function(outliers, icode)
{
#	set outliers missing for specified code(s)
	data <- get(attr(outliers, 'data'))
	if (dim(data)[1] != dim(outliers)[1]) stop('numbers of rows in data and outliers differ')
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

################################################################
#
#	nk2df to convert nk to df (from pre-December 2012)
#
################################################################

	nk2df <- function(obj) {
#	obj is sitar object with arg nk
#	returns obj with nk changed to df
		if (class(obj)[1] != 'sitar') 
			stop(paste("'", deparse(substitute(obj)), "' not a sitar object", sep=""))
		nc <- match('nk', names(obj$call.sitar))
		if (is.na(nc)) stop("'nk' arg not found")
		names(obj$call.sitar)[nc] <- 'df'
		invisible(obj)
	}

################################################################
#
#	revise to update sitar file format (from pre-August 2011)
#
################################################################

	revise <- function(obj) {
#	obj is old-style sitar object
#	returns new-style sitar object
		objname <- deparse(substitute(obj))
		if (grepl('[', objname, fixed=TRUE)) stop('object name must not be indexed')
		if (length(class(obj)) == 1 && class(obj) == 'by') {
			new <- vector('list', length(obj))
			class(new) <- class(obj)
			for (i in 1:length(obj)) {
				if (is.null(obj[[i]]$sitar)) stop('not an old-style sitar object')
				new[[i]] <- obj[[i]]$sitar
				if (!is.null(obj[[i]]$ns)) new[[i]]$ns <- obj[[i]]$ns
				if (!is.null(obj[[i]]$call)) new[[i]]$call.sitar <- obj[[i]]$call
				if (!is.null(obj[[i]]$xoffset)) new[[i]]$xoffset <- obj[[i]]$mux
				if (!'sitar' %in% class(new[[i]])[1]) class(new[[i]]) <- c('sitar', class(new[[i]]))				
				nc <- match('nk', names(new[[i]]$call.sitar))
				if (is.na(nc)) stop("'nk' arg not found")
				names(new[[i]]$call.sitar)[nc] <- 'df'
			}
		}
		else {
			if (is.null(obj$sitar)) stop('not an old-style sitar object')
			new <- obj$sitar
			if (!is.null(obj$ns)) new$ns <- obj$ns
			if (!is.null(obj$call)) new$call.sitar <- obj$call
			if (!is.null(obj$xoffset)) new$xoffset <- obj$mux
			if (!'sitar' %in% class(new)[1]) class(new) <- c('sitar', class(new))
			nc <- match('nk', names(new$call.sitar))
			if (is.na(nc)) stop("'nk' arg not found")
			names(new$call.sitar)[nc] <- 'df'
		}
		new
	}

###################
#	xaxsd  yaxsd  #
###################

	xaxsd <- function(usr=par()$usr[1:2]) (usr + mean(usr) * 0.08) / 1.08
	yaxsd <- function(usr=par()$usr[3:4]) (usr + mean(usr) * 0.08) / 1.08
#	implements xaxs/yaxs option 'd'
#	by default returns xlim/ylim args to match current setting of par()$usr, i.e. previous plot scales
#	e.g. plot(..., xlim=xaxsd(), ylim=yaxsd())
#	specifying usr gives scale with usr args at extremes

##########################
#	sampling for subset  #
##########################

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
	if (!missing(data)) {
		on.exit(detach(data))
		attach(data)
	}
	nx <- length(x)
	subset <- rep(TRUE, nx)
#	simple sampling of whole dataset
	if (is.null(xlim)) {
		ss <- sample(nx, prob * nx)
		subset[-ss] <- FALSE
	} 
	else {
#	sampling of individuals and restricting their age range to xlim
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
#		back-transform transformed x / y scales
#		derive new fitted values for spline curves
#		add acceleration curve option
#		improve ylim for dau plots