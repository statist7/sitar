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
		subset <- eval(mcall$subset, data)
		if (!is.null(subset)) data <- data[subset, ]
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
				knots <- eval(extras$knots)
				df <- length(knots) + 1
			}
			else if (!is.null(extras$df)) {
				df <- eval(extras$df)
				knots <- quantile(x, (1:(df-1))/df)
			}
			if (!is.null(extras$bounds)) {
				bounds <- eval(extras$bounds)
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
