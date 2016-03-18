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
#	update formulae
	if (length(extras) > 0) {
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
#	check if can use start
	if (sum(pmatch(names(extras), c("x", "y", "id", "fixed", "random", "a.formula", "b.formula", "c.formula", "start", "returndata")), na.rm=TRUE) == 0) {
  	start. <- list(fixed=fixef(object), random=ranef(object))
#	args data and subset
		data <- eval(mcall$data)
		subset <- eval(mcall$subset, data)
		if (!is.null(subset)) data <- data[subset, ]
		id <- factor(eval(mcall$id, data))
		levels.obj <- levels(getGroups(object))
		if (!identical(levels(id), levels.obj)) {
#	omit random effects for missing levels in id
			idcheck <- levels.obj %in% levels(id)
			start.$random <- ranef(object)[idcheck,]
			cat(length(levels.obj) - sum(idcheck), 'subjects omitted\n')
#	add zero random effects for new levels in id
			newid <- !levels(id) %in% levels.obj
			if (sum(newid) > 0) {
				newre <- matrix(0, nrow=sum(newid), ncol=dim(ranef(object))[2], dimnames=list(levels(id)[newid], dimnames(ranef(object))[[2]]))
				start.$random <- rbind(start.$random, newre)
				cat(sum(newid), 'subjects added\n')
			}
		}
#	update start fixed effects if df, knots, bounds or bstart updated
		if (sum(pmatch(names(extras), c("df", "knots", "bounds", "bstart")), na.rm=TRUE) > 0) {
			x <- eval(mcall$x, data)
			df <- object$ns$rank - 1
			knots <- attr(object$ns$model$ns, 'knots')
			bounds <- attr(object$ns$model$ns, 'Boundary.knots')
			if (length(fixef(object)) > df + 1) fixed.extra <- (df+2):length(fixef(object))
				else fixed.extra <- NULL
# arg knots
			if (!is.null(extras$knots)) {
				knots <- eval(extras$knots) - mean(x)
				df <- length(knots) + 1
				mcall$df <- NULL
			}
# arg df
			else if (!is.null(extras$df)) {
				df <- eval(extras$df)
				knots <- quantile(x, (1:(df-1))/df) - mean(x)
				mcall$knots <- NULL
			}
# arg bounds
			if (!is.null(extras$bounds)) {
				bounds <- eval(extras$bounds)
				if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x))
				bounds <- bounds - mean(x)
			}
# arg bstart
			if (!is.null(extras$bstart)) {
				bstart <- eval(extras$bstart)
				if (is.character(bstart)) bstart <- mean(x)
				start.$fixed['b'] <- bstart
			}
#	get spline start values
			spline.lm <- lm(predict(object, data, level=0) ~ ns(x - mean(x), knots=knots, Bound=bounds))
			start.$fixed <- c(coef(spline.lm)[c(2:(df+1), 1)], start.$fixed[fixed.extra])
		}
#	save start. object
		assign('start.', start., parent.frame())
		mcall <- as.call(c(as.list(mcall), start=quote(start.)))
	}
	if (evaluate)
		eval(mcall, parent.frame())
	else mcall
}
