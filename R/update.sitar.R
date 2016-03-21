	update.sitar <- function (object, ..., evaluate = TRUE)
{
	mcall <- object$call.sitar
	if (is.null(mcall))
		stop("need an object with call.sitar component")
	extras <- as.list(match.call(expand.dots = FALSE)$...)
#	drop null args
	mcall$start <- NULL
	for (i in names(extras))
		if (is.null(extras[[i]]))
			mcall[[i]] <- extras[[i]] <- NULL
# update args
	mcall[names(extras)] <- extras
#	add start arg if none of these args specified
	if (!sum(pmatch(names(extras), c("x", "y", "id", "fixed", "random", "a.formula", "b.formula",
	                                "c.formula", "start", "returndata"), 0))) {
  	start. <- list(fixed=fixef(object), random=ranef(object))
# update start if any of these args specified
  	if (sum(pmatch(names(extras), c('data', 'subset', 'df', 'knots', 'bounds', 'bstart'), 0))) {
# get data etc
  		data <- eval(mcall$data)
  		subset <- eval(mcall$subset, data)
  		if (!is.null(subset)) data <- data[subset, ]
  		x <- eval(mcall$x, data)
  		df <- object$ns$rank - 1
  		knots <- attr(object$ns$model$ns, 'knots')
  		bounds <- attr(object$ns$model$ns, 'Boundary.knots')
# update random effects
  		if (!is.null(extras$data) || !is.null(extras$subset)) {
  		  id <- factor(eval(mcall$id, data))
  		  levels.obj <- levels(getGroups(object))
  		  if (!identical(levels(id), levels.obj)) {
#	omit random effects for missing levels in id
  		    start.$random <- start.$random[idcheck <- levels.obj %in% levels(id), ]
  		    cat(length(levels.obj) - sum(idcheck), 'subjects omitted\n')
#	add zero random effects for new levels in id
  		    newid <- !levels(id) %in% levels.obj
  		    if (sum(newid) > 0) {
  		      newre <- matrix(0, nrow=sum(newid), ncol=dim(ranef(object))[2], dimnames=list(levels(id)[newid], dimnames(ranef(object))[[2]]))
  		      start.$random <- rbind(start.$random, newre)
  		      cat(sum(newid), 'subjects added\n')
  		    }
  		  }
  		}
#	update fixed effects
  		if (length(fixef(object)) > df + 1) fixed.extra <- (df+2):length(fixef(object))
  			else fixed.extra <- NULL
# new arg knots
  		if (!is.null(extras$knots)) {
  			knots <- eval(extras$knots) - mean(x)
  			df <- length(knots) + 1
  			mcall$df <- NULL
  		}
# new arg df
  		else if (!is.null(extras$df)) {
  			df <- eval(extras$df)
  			knots <- quantile(x, (1:(df-1))/df) - mean(x)
  			mcall$knots <- NULL
  		}
# new arg bounds
  		if (!is.null(extras$bounds)) {
  			bounds <- eval(extras$bounds)
  			if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x))
  			bounds <- bounds - mean(x)
  		}
# new arg bstart
  		if (!is.null(extras$bstart) && !is.null(start.$fixed['b'])) {
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
