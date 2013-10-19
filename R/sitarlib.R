#	global variables
	if (getRversion() >= "3.0.0") utils::globalVariables(c('.par.usr2'))

#############################
#
#	sitar
#
#############################

	sitar <- function(x, y, id, data, df, knots, fixed=random, random='a+b+c', a.formula=~1, b.formula=~1, c.formula=~1, bounds=0.04, start, bstart='mean', xoffset='mean', file='.sitar.temp.R', returndata=FALSE, verbose=FALSE, correlation=NULL, weights=NULL, subset=NULL, method='ML', na.action=na.fail, control = nlmeControl(returnObject=TRUE), newform=TRUE)
#
#	fit growth curves of y ~ ns(f(x)) by id
#	df - number of knots (or length of knots + 1)
#	knots - default x quantiles
#	fixed - default random
#	random - random effects, default "a+b+c"
#	a.formula etc -  fixed effect formulae, default ~ 1
#			to omit fixed effects use ~ -1
#	bounds - span of ns, default x-range ±4%
#	start - starting values - default estimated
#			requires spline coefficients, any missing zeroes added
#	bstart - starting value for b, default 'mean', or 'apv' or value 
#		(subsumes xoffset)
#	xoffset - offset for x, default 'mean', alternatives 'apv' or value 
#	file - filename for temporary source code
#	returndata - if TRUE returns nlme data frame not nlme model
#	verbose - show code in file
#	verbose etc - arguments passed to nlme

#	newform - FALSE uses xoffset, TRUE uses bstart

#	returns call, ns, nlme.out if returndata FALSE
#			else fulldata
{
	b.origin <- function(b) {
		if (b == 'mean') return(mean(x))
		if (b == 'apv') {
			spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))
			return(makess(x, fitted(spline.lm))$apv[1])
		}
		if (!is.numeric(b)) stop('must be "mean" or "apv" or numeric')
		b
	}
	mcall <- match.call()
	if (!missing(data)) {
		if (!deparse(substitute(data), width.cutoff=99) %in% search()) {
			on.exit(detach(data))
			attach(data)
		}
	}
	if (missing(df) & missing(knots)) stop("either df or knots must be specified")
	if (!missing(df) & !missing(knots)) cat("both df and knots specified - df redefined from knots\n")
	if (missing(knots)) knots <- quantile(x, (1:(df-1))/df) 
		else df <- length(knots) + 1
#	define bounds, default x range ±4% 
	if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x)) 
	if (length(bounds) != 2) stop("bounds should be length 1 or 2")
#	derive x variable offset
	if (newform || !missing(bstart)) { # using bstart
		newform <- TRUE
		mcall$xoffset <- NULL
		if (b.formula == as.formula('~ -1') || b.formula == as.formula('~ 1-1') || !grepl('b', fixed)) bstart <- 0 
		else bstart <- b.origin(bstart)
		# else if (bstart == 'mean') bstart <- mean(x)
		# else if (bstart == 'apv') {
			# spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))
			# bstart <- makess(x, fitted(spline.lm))$apv[1]
		# }
		# else if (!is.numeric(bstart)) stop('bstart must be "mean" or "apv" or numeric')
		# xc <- x
		knots <- knots - bstart
		bounds <- bounds - bstart
#	get spline start values
		spline.lm <- lm(y ~ ns(x - bstart, knots=knots, Bound=bounds))
	} 
	else { # using xoffset
		xoffset <- b.origin(xoffset)
		# if (xoffset == 'mean') xoffset <- mean(x)
		# else if (xoffset == 'apv') {
			# spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))
			# xoffset <- makess(x, fitted(spline.lm))$apv[1]
		# }
		# else if (!is.numeric(xoffset)) stop('xoffset must be "mean" or "apv" or numeric')
#	apply xoffset	
		x <- x - xoffset
		knots <- knots - xoffset
		bounds <- bounds - xoffset
#	get spline start values
		spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))
	}
#	get starting values for ss and a
	if (missing(start)) start <- coef(spline.lm)[c(2:(df+1), 1)]
#	force fixed effect for a
	fix <- fixed
	if (!grepl('a', fix)) fix <- paste('a', fix, sep='+')
	fixed <- ss <- paste('s', 1:df, sep='')
	pars <- c('x', ss)
#	set up model elements for a, b and c
	names(model) <- model <- letters[1:3]
	constant <- mm.formula <- as.formula('~ 1')
	# mm.intercept <- TRUE
	if (is.null(subset)) subset <- 1:length(x)
	fulldata <- data.frame(x, y, id, subset)
	for (l in model) {
		if (!grepl(l, fix) && !grepl(l, random)) {
			model[l] <- NA
			next
		}
		pars <- c(pars, l)
		formula <- get(paste(l, 'formula', sep='.'))
		if (formula == as.formula('~ -1') || formula == as.formula('~ 1-1') || !grepl(l, fix)) next
		if (formula == constant) mm.intercept <- TRUE
		else {
			if (formula != mm.formula) {
				mm <- model.matrix(formula, data)
				if (dim(mm)[[1]] < length(y)) 
					stop('Missing values in data')
				mm.formula <- formula
				#	convert spaces to underlines
				colnames(mm) <- gsub(' ', '_', colnames(mm))
				#	convert colons to dots (from interaction)
				colnames(mm) <- gsub(':', '.', colnames(mm), fixed=TRUE)
				mm.intercept <- colnames(mm)[1] == '(Intercept)'
				# omit constant columns and centre others
				cmm <- apply(mm, 2, function(x) sd(x, na.rm=TRUE) > 0)
				mm.names <- colnames(mm)[cmm]
				mm <- as.matrix(mm[, mm.names])
				colnames(mm) <- mm.names
				mm <- scale(mm, scale=FALSE)
			}
			if (exists('mm')) for (i in 1:dim(mm)[[2]]) {
				var <- colnames(mm)[i]
				rc <- paste(l, var, sep='.')
				pars <- c(pars, rc)
				fixed <- c(fixed, rc)
				if (!is.list(start)) start <- c(start, 0)
				model[l] <- paste(model[l], '+', rc, '*', var, sep='')
				if (!var %in% pars) {
					pars <- c(pars, var)
					fulldata <- cbind(fulldata, mm[,i])
					names(fulldata)[dim(fulldata)[2]] <- var
					# assign(var, mm[,i])
				}					
			}		
		}
		if (mm.intercept) {
			fixed <- c(fixed, l)
			if (!is.list(start)) {
				if (newform) {
					if (l == 'b') start <- c(start, bstart)
					else if (l == 'c') start <- c(start, 0)					
				}
				else if (l != 'a') start <- c(start, 0)
			}
		}
	}
	if (returndata) invisible(fulldata) else {
		pars <- paste(pars, collapse=',')
		fixed <- paste(fixed, collapse='+')
		sscomma <- paste(ss, collapse=',')
	
	#	combine model elements
		nsd <- paste(model['a'], '+')
		nsf <- paste('(x', ifelse(!is.na(model['b']), paste('- (', model['b'], '))'), ')'))
		if (!is.na(model['c'])) nsf <- paste(nsf, '* exp(', model['c'], ')')
	
	#	print values
		if (verbose) {
			cat('df', df, 'bstart', bstart, 'xoffset', xoffset, '\nknots\n', knots, '\nbounds\n', bounds)
			if (is.list(start)) {
				cat('\nstarting values\n  fixed effects\n', start$fixed) 
				if (!is.null(start$random)) {
					cat('\n  random effects\n')
					print(start$random)
				}
			}
			else cat('\nstarting values\n', start, '\n')
		}
	
	#	create file to source		
		code <- c(
	"	fitnlme <- function($pars) {",
	"#		** created by sitar - edits will be lost **",
	"		splinecoefs <- as.matrix(cbind($sscomma))",
	"		as.vector( $nsd",
	"		(splinecoefs * as.matrix(ns($nsf,",
	"			knots=knots, Boundary.knots=bounds))) %*% ",
	"			matrix(rep(1,df), ncol=1))",
	"	}",
	"	assign('fitnlme', fitnlme, globalenv())",
	"	nlme.out <- nlme(y ~ fitnlme($pars),",
	"	fixed = $fixed ~ 1,",
	"	random = $random ~ 1 | id,",
	"	data = fulldata,",
	"#	groups, (unimplemented nlme argument)",
	"	start = start, correlation = correlation,",
	"	weights = weights, subset = subset, method = method,",
	"	na.action = na.action, control = control, verbose = verbose)")
	
		for (i in c('random', 'pars', 'fixed', 'sscomma', 'nsd', 'nsf', 'file')) {
			# if (verbose) cat(paste(i,get(i), '\n'))
			code <- gsub(paste('$', i, sep=''), get(i), code, fixed=TRUE)	}
		
	#	better to write to and source from clipboard, but machine-dependent
		write(code, file=file, sep='\n')
		if (verbose) {
			cat('file', file, '\n')
			cat(code, sep='\n')
		}
		elapsed <- system.time(source(file, local=TRUE))
	
		nlme.out$call.sitar <- mcall
		if (newform) nlme.out$bstart <- bstart
			else nlme.out$xoffset <- xoffset
		nlme.out$ns <- spline.lm
		nlme.out$elapsed <- elapsed
		if (class(nlme.out)[1] == 'nlme') class(nlme.out) <- c('sitar', class(nlme.out))
		nlme.out
	}
}

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
	on.exit(detach(data))
	attach(data, warn.conflicts=FALSE)
	y <- eval(mcall$y)
	nf <- length(fitted(object))
	if (nf != length(y)) stop(paste('model (length=', nf, ') incompatible with data (rows=', length(y), ')', sep=''))
	x <- eval(mcall$x)
	id <- eval(mcall$id)
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

	plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun, yfun, subset=NULL, abc=c(a=0, b=0, c=0), add=FALSE, xy=NULL, nlme=FALSE, ...)
#	plot curves from sitar model
#	opt: 
#		d = fitted distance curve (labels[1] = x, labels[2] = y)
#		v = fitted velocity curve (labels[3] = y)
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
		if (!is.null(mcall$subset)) data <- data[eval(mcall$subset),]
		on.exit(detach(data))
		attach(data)
#	rm(xy) needed to eval xy$apv at top level
		# rm(xy)
		x <- eval(mcall$x)
		y <- eval(mcall$y)
		id <- factor(eval(mcall$id)) # factor added 23/4/13
		nf <- length(fitted(model))
		if (nf != length(y)) stop(paste('model (length=', nf, ') incompatible with data (rows=', length(y), ')', sep=''))
		if (is.null(subset)) subset <- rep(TRUE, nf)
	
		ARG <- list(...)
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
		
#	plot y vs t by subject
		if (grepl("u", opt)) {
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
			do.call("mplot", c(list(x=quote(x), y=quote(y), id=quote(id), subset=subset, add=add), ARG))
			add <- TRUE
		}
	
#	plot fitted curves by subject
		if (grepl("f", opt)) {
			y <- fitted(model, level=1)
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
			do.call("mplot", c(list(x=quote(x), y=quote(y), id=quote(id), subset=subset, add=add), ARG))
			add <- TRUE
		}
	
#	plot adjusted y vs adjusted t by subject
		if (grepl("a", opt)) {
			fred <- summary(model)
			x <- fred$x.adj
			y <- fred$y.adj
			if (!missing(xfun)) x <- xfun(x)
			if (!missing(yfun)) y <- yfun(y)
			do.call("mplot", c(list(x=quote(x), y=quote(y), id=quote(id), subset=subset, add=add), ARG))
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
#	need to add "subset"
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
			c <- start.$fixed['c']
			if (is.na(c)) c <- 0
			x <- (x - bstart) * exp(c)
			knots <- knots - bstart
			bounds <- bounds - bstart
#	get spline start values
			spline.lm <- lm(fitted(object, level=0) ~ ns(x, knots=knots, Bound=bounds))
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
#	mplot
#
#############################

	mplot <- function(x, y, id, data=NULL, subset=NULL, add=FALSE, ...)
#	plots y ~ x by id with data
#	subset defines a subset of rows
#	add TRUE suppresses plot axes
#	... parameters where col, lty, lwd, pch can depend on id
{	
#	save names
	xl <- deparse(substitute(x))
	yl <- deparse(substitute(y))
	
#	attach or create data
	if (!is.null(data)) {
		data <- data[, c(xl, yl, deparse(substitute(id)))]
		on.exit(detach(data))
		attach(data)
		xv <- substitute(x)
		yv <- substitute(y)
	} 
	else {
#	save values
		xv <- substitute(x, globalenv())
		yv <- substitute(y, globalenv())
		data <- as.data.frame(cbind(x, y, id))
	}
#	extract and save vector ... args: col lty lwd pch
	ARG <- list(...)
	cnames <- names(ARG)
	if (!is.null(cnames)) {
		cnames <- cnames[lapply(ARG, length) == dim(data)[[1]]]
		data[, cnames] <- ARG[cnames]
		ARG[cnames] <- NULL
	}
#	handle subset
	if (is.null(subset)) subset <- rep(TRUE, length(id))
	subset <- ifelse(is.na(x) | is.na(y), FALSE, subset)
	if (sum(subset, na.rm=TRUE) == 0) stop("no data to plot")

#	if new graph plot axes
	if (!add) do.call...("plot", list(x=eval(x)[subset], y=eval(y)[subset], type='n', xlab=xl, ylab=yl), ...)

#	draw growth curves
	tt <- by(data[subset,], id[subset], function(z) {
#	sort by x
		zs <- z[order(eval(xv, z)),]
		xvt <- eval(xv, zs)
		yvt <- eval(yv, zs)
#	restore vector ... args
		if (!is.null(cnames)) ARG[cnames] <- as.list(as.data.frame(zs[, cnames]))

#	lines(xvt, yvt, ...)
		do.call("lines", c(list(x=xvt, y=yvt), ARG))
	}	)
}

#############################
#
#	y2plot
#
#############################

	y2plot <- function(x, y1, y2=NULL, labels, y2par=NULL, add=FALSE, xy, xlegend="topleft", inset=0.05, ...)
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
	if (is.null(y2par$ylab)) y2par$ylab <- labels[3]
#	plot y1
	ypar <- list(...)
#	save y1 par args
	opar <- par(no.readonly=TRUE)
	lty <- 1:2
	lwd <- c(1, 1)
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
		do.call('plot', c(list(x=x, y=y1, type='l'), ypar))
#	save x & y1 axis limits
		xy$usr <- par('usr')
#	optionally plot y2 and add right axis, with y2par as ...
		if (!missing(y2)) {
			par(new=TRUE)
#	ensure x axis same for y2 as y1
			if (!is.null(ypar$xlim) && is.null(y2par$xlim)) y2par$xlim <- ypar$xlim
			do.call('plot', c(list(x=x, y=y2, ann=FALSE, bty="n", xaxt="n", yaxt="n", type="l"), y2par))
#	save y2 axis limits
			xy$usr2 <- par('usr')
			assign('.par.usr2', par('usr'), globalenv())
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
				assign('.par.usr2', par('usr'), globalenv())
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

	plotclean <- function(x, y=NULL, id=NULL, data=NULL, n=length(x), pch=20, ...)
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
			lines(x[id==id[ans]], y[id==id[ans]])
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
#	devadj
#
#############################

	devadj <- function(..., pattern=NULL)
#	input SITAR model(s)
#	returns deviance for power transformed y
#	looks for log, sqrt or ^ - assumes no multipliers
{	ARG <- list(...)
	if (!is.null(pattern)) {
		pattern <- ls(envir=parent.frame(), pattern=pattern)
		ARG <- c(ARG, lapply(as.list(pattern), get))
	}
	dev <- lapply(ARG, function(obj) {
		if (!'call' %in% names(obj)) return(NA)
#	check for call.sitar$y or else call$y or else call$formula
		if (!is.null(obj$call.sitar)) obj$call <- obj$call.sitar
#	check for call$y or call$model
		if (!is.null(obj$call$y)) ycall <- obj$call$y 
			else if (!is.null(obj$call$model)) ycall <- obj$call$model[[2]] 
			else if (!is.null(obj$call$formula)) ycall <- obj$call$formula[[2]] 
				else return(NA)
		lambda <- 1
		if (length(ycall) == 1) yt <- ycall else {
			yt <- as.symbol(ycall[[2]])
			fun <- ycall[1]
			if (fun == "log()") lambda <- 0 else
			if (fun == "sqrt()") lambda <- 0.5 else 
			if (fun == "`^`()") lambda <- eval(ycall[[3]])
		}
		y <- eval(yt, eval(obj$call$data, sys.parent()))
		sly <- ifelse(lambda == 1, 0, sum(log(y)))
		-2 * (obj$logLik + (lambda-1) * sly + length(y) * log(abs(lambda) + (lambda == 0)))		
	})
	names(dev) <- c(match.call(expand.dots=FALSE)$..., pattern)
	dev <- unlist(dev[!is.na(dev)])
	if (length(dev) == 0) return(invisible())
	dev[order(dev)]
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
		mc <- lm(ss1$y[o] ~ poly(ss$x[o], 2, raw=TRUE))$coef # velocity as quadratic in age
		apv[1] <- - mc[2] / 2 / mc[3] # age at peak velocity
		apv[2] <- mc[1] - mc[2]^2 / 4 / mc[3] # peak velocity
		names(apv) <- c('apv', 'pv')
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
		do.call...('plot', list(z[,2], z[,3], type='b', xlab=x, ylab=y, pch=46), ...)
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

	zLMS <- function(x, L, M, S, data=NULL) {
		with(data, {
			L0 <- L + 1e-7 * (L == 0)
			( (x / M) ^ L0 - 1) / L0 / S
		} )
	}
	
	cLMS <- function(z, L, M, S, data=NULL) {
		with(data, {
			L0 <- L + 1e-7 * (L == 0)
			c <- M * (1 + L0 * S %o% z) ^ (1 / L0)
			if (length(z) == 1) as.numeric(c)
				else c
		} )
	}
	
	lms2z <- function(x, y, sex, data=NULL, measure, ref, toz=TRUE) {
#	converts measurement y to/from z-score adjusted for x & sex
#		using LMS reference 'ref' for 'measure'
#	x		age
#	y		measurement (or z-score if toz FALSE)
#	sex		sex variable (male=1, female=2)
#	data	source of x, y and sex
#	measure label for measurement, one of:
#		'ht' 'wt' 'bmi' 'head' 'sitht' 'leglen' 'waist' 'bfat'
#	ref		name of reference, one of: 'uk90' 'who06'
#	toz		if TRUE returns measurement converted to z-score using ref
#			if FALSE returns z-score converted to measurement using ref
	on.exit(detach(data))
	attach(data)
	lms <- paste(c('L','M','S'), measure, sep='.')
	lmsout <- matrix(nrow=length(x), ncol=3)
	ref <- get(ref)
	for (i in 1:3) {
		for (ix in 1:2) {
			if (sum(as.numeric(sex) == ix, na.rm=TRUE) > 0)
			lmsout[as.numeric(sex) == ix, i] <- spline(ref$years[ref$sex == ix], ref[ref$sex == ix, lms[i]], method='natural', xout=x[as.numeric(sex) == ix])$y
		}
	}
	if (toz) zLMS(y, lmsout[,1], lmsout[,2],lmsout[,3])
	else cLMS(y, lmsout[,1], lmsout[,2],lmsout[,3])
}	
	z2cent <- function(z)
#	z is z-score
#	returns corresponding centile as label
{
	np <- ifelse(abs(z) < 2.33, 0, 1)
	ct <- round(pnorm(z) * 100, np)
	mod10 <- ifelse(np == 1, 0, floor(ct %% 10))
	th <- ifelse(mod10 == 0 | mod10 > 4, 4, mod10)
	th <- paste(ct, c('st','nd','rd','th')[th], sep='')
	th[th == '0th'] <- paste('SDS', round(z[th == '0th'], 1), sep='')
	th[th == '100th'] <- paste('SDS', round(z[th == '100th'], 1), sep='+')
	th
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
#		omit long line segments from adjusted plot
#		derive new fitted values for spline curves
#		sort out mar problem
#		add acceleration curve option
#		improve ylim for dau plots
#		provide lines.sitar 									# added 04/01/13
#		change newplot (default TRUE) to add (default FALSE)	# added 05/04/13