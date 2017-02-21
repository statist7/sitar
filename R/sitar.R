#' Fit SITAR growth curve model
#'
#' SITAR is a method of growth curve analysis, based on \pkg{nlme}, that
#' summarises a set of growth curves with a mean growth curve as a regression
#' spline, plus a set of up to three fixed and random effects (a, b and c)
#' defining how individual growth curves differ from the mean curve.
#'
#' \code{xoffset} allows the origin of \code{x} to be varied, while
#' \code{bstart} specifies the starting value for \code{b}, both of which can
#' affect the model fit and particularly \code{b}. The values of \code{bstart},
#' \code{knots} and \code{bounds} are offset by \code{xoffset} for fitting
#' purposes, and similarly for fixed effect \code{b}.
#'
#' The formulae \code{a.formula}, \code{b.formula} and \code{c.formula} can
#' include functions and interactions, but \code{\link{make.names}} is used to
#' ensure that the names of the corresponding model terms are valid. The
#' modified not the original names need to be specified in \code{predict.sitar}.
#'
#' \code{update} updates the model by taking the \code{object} call, adding any
#' new parameters and replacing changed ones. Where feasible the fixed and
#' random effects of the model being updated are suitably modified and passed
#' via the \code{start} argument.
#'
#' @aliases sitar update.sitar
#' @param x vector of ages.
#' @param y vector of measurements.
#' @param id factor of subject identifiers.
#' @param data data frame containing variables \code{x}, \code{y} and
#' \code{id}.
#' @param df degrees of freedom for cubic regression spline (2 or more).
#' @param knots vector of values for knots (default \code{df} quantiles of
#' \code{x} distribution).
#' @param fixed character string specifying a, b, c fixed effects (default
#' \code{random}).
#' @param random character string specifying a, b, c random effects (default
#' \code{"a+b+c"}).
#' @param a.formula formula for fixed effect a (default \code{~ 1}).
#' @param b.formula formula for fixed effect b (default \code{~ 1}).
#' @param c.formula formula for fixed effect c (default \code{~ 1}).
#' @param bounds span of \code{x} for regression spline, or fractional
#' extension of range (default 0.04).
#' @param start optional numeric vector of initial estimates for the fixed
#' effects, or list of initial estimates for the fixed and random effects (see
#' \code{\link{nlme}}).
#' @param xoffset optional value of offset for \code{x} (either "mean"
#' (default), "apv" or value).
#' @param bstart optional starting value for fixed effect \code{b} (either
#' "mean", "apv" or value (default \code{xoffset})).
#' @param returndata logical which if TRUE causes the model matrix to be
#' returned, or if FALSE (default) the fitted model. Setting returndata TRUE is
#' useful in conjunction with \code{subset} and \code{\link{subsample}} for
#' simulation purposes.
#' @param verbose optional logical value to print information on the evolution
#' of the iterative algorithm (see \code{\link{nlme}}).
#' @param correlation optional \code{corStruct} object describing the
#' within-group correlation structure (see \code{\link{nlme}}).
#' @param weights optional \code{varFunc} object or one-sided formula
#' describing the within-group heteroscedasticity structure (see
#' \code{\link{nlme}}).
#' @param subset optional expression indicating the subset of the rows of data
#' that should be used in the fit (see \code{\link{nlme}}).
#' @param method character string, either "REML" or "ML" (default) (see
#' \code{\link{nlme}}).
#' @param na.action function for when the data contain NAs (see
#' \code{\link{nlme}}).
#' @param control list of control values for the estimation algorithm (see
#' \code{\link{nlme}}) (default {nlmeControl(returnObject=TRUE)}).
#' @param object object of class \code{sitar}.
#' @param \dots further parameters for \code{update} consisting of any of the
#' above \code{sitar} parameters.
#' @param evaluate logical to control evaluation.  If TRUE (default) the
#' expanded \code{update} call is passed to \code{sitar} for evaluation, while
#' if FALSE the expanded call itself is returned.
#' @return An object inheriting from class \code{sitar} representing the
#' nonlinear mixed-effects model fit, with all the components returned by
#' \code{nlme} (see \code{\link{nlmeObject}} for a full description) plus the
#' following components:
#' \item{fitnlme}{the function returning the predicted value of \code{y}.}
#' \item{call.sitar}{the internal \code{sitar} call that produced the object.}
#' \item{xoffset}{the value of \code{xoffset}.}
#' \item{ns}{the \code{lm} object providing starting values for the B-spline curve.}
#'
#' Generic functions such as \code{print}, \code{plot}, \code{anova} and
#' \code{summary} have methods to show the results of the fit. The functions
#' \code{resid}, \code{coef}, \code{fitted}, \code{fixed.effects},
#' \code{random.effects}, \code{predict}, \code{getData}, \code{getGroups},
#' \code{getCovariate} and \code{getVarCov} can be used to extract some of its
#' components.
#'
#' Note that versions of \code{sitar} prior to 1.0.4 did not return
#' \code{fitnlme}. Both \code{plot} and \code{predict} may require it, in which
#' case they \code{update} the SITAR object on the fly, with a message. Also
#' version 1.0.5 altered the defaults for \code{xoffset} and \code{bstart}.
#' Models fitted with versions prior to 1.0.5 need refitting.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @keywords package nonlinear regression models
#' @examples
#'
#' data(heights)
#' ##  fit simple model
#' (m1 <- sitar(x=age, y=height, id=id, data=heights, df=5))
#'
#' ##  relate random effects to age at menarche (with censored values +ve)
#' ##  both a (size) and b (tempo) are positively associated with age at menarche
#' amen <- abs(heights$men)
#' (m2 <- update(m1, a.form=~amen, b.form=~amen, c.form=~amen))
#' @import nlme splines
#' @importFrom stats AIC BIC as.formula coef cor fitted lag lm logLik
#' mad model.frame model.matrix na.fail na.omit pnorm predict qnorm quantile
#' residuals sd setNames smooth.spline spline update update.formula
#' @export sitar
	sitar <- function(x, y, id, data, df, knots, fixed=random, random='a+b+c',
	                  a.formula=~1, b.formula=~1, c.formula=~1, bounds=0.04, start, xoffset='mean', bstart=xoffset,
	                  returndata=FALSE, verbose=FALSE, correlation=NULL, weights=NULL, subset=NULL, method='ML',
	                  na.action=na.fail, control = nlmeControl(returnObject=TRUE))
{
	b.origin <- function(b) {
		if (b == 'mean') return(mean(x))
		if (b == 'apv') {
			spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))
			return(getPeakTrough(x, predict(smooth.spline(x, fitted(spline.lm)), x, deriv=1)$y)[1])
		}
		if (!is.numeric(b)) stop('must be "mean" or "apv" or numeric')
		b
	}

# get data
	mcall <- match.call()
	data <- eval(mcall$data, parent.frame())
	subset <- eval(mcall$subset, data)
	if (!is.null(subset))
	  data <- data[subset, ]
	x <- eval(mcall$x, data)
	y <- eval(mcall$y, data)

# get df, knots and bounds
	if (missing(df) & missing(knots)) stop("either df or knots must be specified")
	if (!missing(df) & !missing(knots)) cat("both df and knots specified - df redefined from knots\n")
	if (missing(knots)) {
	  if (df < 2) stop("df must be 2 or more")
	  knots <- quantile(x, (1:(df-1))/df)
	} else {
	  df <- length(knots) + 1
	}
	if (nrow(data) <= df) stop("too few data to fit spline curve")
	if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x))
	else if (length(bounds) != 2) stop("bounds should be length 1 or 2")

	xoffset <- b.origin(xoffset)
	bstart <- b.origin(bstart) - xoffset
	x <- x - xoffset
	knots <- knots - xoffset
	bounds <- bounds - xoffset

# get spline start values
  spline.lm <- lm(y ~ ns(x, knots=knots, Bound=bounds))

#	if start missing get start values for ss and a
  if (nostart <- missing(start)) start <- coef(spline.lm)[c(2:(df+1), 1)]

#	force fixed effect for a
	fix <- fixed
	if (!grepl('a', fix)) fix <- paste('a', fix, sep='+')

# set up args for fitnlme
	fixed <- ss <- paste0('s', 1:df)
	pars <- c('x', ss)

# if subsetted restore data
  if (!is.null(subset)) {
    data <- eval(mcall$data, parent.frame())
    x <- eval(mcall$x, data) - xoffset
    y <- eval(mcall$y, data)
  } else {
    subset <- 1:length(x)
  }
  id <- eval(mcall$id, data)
  fulldata <- data.frame(x, y, id, subset)

#	set up model elements for a, b and c
	names(model) <- model <- letters[1:3]
	constant <- mm.formula <- as.formula('~ 1')
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
				if (nrow(mm) < length(y))
					stop('Missing values in data')
				mm.formula <- formula
				mm.intercept <- colnames(mm)[1] == '(Intercept)'
# ensure names are valid
      	colnames(mm) <- make.names(colnames(mm), unique=TRUE)
# omit constant columns
				mm <- mm[, apply(mm, 2, function(x) max(x) > min(x)), drop=FALSE]
# centre columns
				mm <- scale(mm, scale=FALSE)
			}
			if (exists('mm')) for (i in 1:ncol(mm)) {
				var <- colnames(mm)[i]
				rc <- paste(l, var, sep='.')
				pars <- c(pars, rc)
				fixed <- c(fixed, rc)
				if (nostart) start <- c(start, 0)
				model[l] <- paste0(model[l], '+', rc, '*', var)
				if (!var %in% pars) {
					pars <- c(pars, var)
					fulldata <- cbind(fulldata, mm[, i, drop=FALSE])
				}
			}
		}
		if (mm.intercept) {
			fixed <- c(fixed, l)
			if (nostart) {
			  if (l == 'b') start <- c(start, bstart)
			    else if (l == 'c') start <- c(start, 0)
			}
		}
	}

# 	if (!is.null(weights)) {
#     if (is.list(weights)) form <- asOneFormula(lapply(weights, function(z) attr(z, 'formula')))
#       else form <- attr(weights, 'formula')
#     if (!is.null(form)) {
#       wt <- as.data.frame(model.frame(form, data))
#       wt.names <- names(wt)[!names(wt) %in% names(fulldata)]
#       wt <- as.data.frame(wt[, wt.names])
#       names(wt) <- wt.names
#       fulldata <- cbind(fulldata, wt)
# 	  }
# 	}

	if (returndata) invisible(fulldata) else {
		pars <- paste(pars, collapse=',')
		fixed <- paste(fixed, collapse='+')
		sscomma <- paste(ss, collapse=',')

#	combine model elements
		nsd <- paste(model['a'], '+')
		nsf <- paste('(x', ifelse(!is.na(model['b']), paste('- (', model['b'], '))'), ')'))
		if (!is.na(model['c'])) nsf <- paste(nsf, '* exp(', model['c'], ')')

#	code to parse
  	fitcode <- c(
	"fitenv <- new.env()",
	"fitenv$fitnlme <- function($pars) {",
	"as.vector( $nsd",
	"(as.matrix(cbind($sscomma)) * as.matrix(ns($nsf,",
	"knots=knots, Boundary.knots=bounds))) %*%",
	"matrix(rep(1,df), ncol=1))",
	"}",
	"on.exit(detach(fitenv))",
	"attach(fitenv)",
	"nlme(y ~ fitnlme($pars),",
	"fixed = $fixed ~ 1,",
	"random = $random ~ 1 | id,",
	"data = fulldata,",
	"start = start, correlation = correlation,",
	"weights = weights, subset = subset, method = method,",
	"na.action = na.action, control = control, verbose = verbose)")

		for (i in c('random', 'pars', 'fixed', 'sscomma', 'nsd', 'nsf'))
  		fitcode <- gsub(paste0('$', i), get(i), fitcode, fixed=TRUE)

#	print values
		if (verbose) {
			cat('\nconstructed code', fitcode, sep='\n')
			cat('\ndf', df, 'bstart', bstart, 'xoffset', xoffset, '\nknots\n', knots, '\nbounds\n', bounds)
			if (is.list(start)) {
				cat('\nstarting values\n  fixed effects\n', start$fixed)
				if (!is.null(start$random)) {
					cat('\n  random effects\n')
					print(start$random)
				}
			}
			else cat('\nstarting values\n', start, '\n')
		}

#	save fitted model
    nlme.out <- eval(parse(text=fitcode))
#     if (exists('start.')) rm(start., inherits=TRUE)
    nlme.out$fitnlme <- fitenv$fitnlme
		nlme.out$call.sitar <- mcall
		nlme.out$xoffset <- xoffset
		nlme.out$ns <- spline.lm
		if (!'sitar' %in% class(nlme.out)) class(nlme.out) <- c('sitar', class(nlme.out))
		nlme.out
	}
}

#' @rdname sitar
#' @method update sitar
#' @export
update.sitar <- function (object, ..., evaluate = TRUE)
{
  mcall <- object$call.sitar
  if (is.null(mcall))
    stop("need an object with call.sitar component")
  extras <- as.list(match.call(expand.dots = FALSE)$...)
  mcall$start <- NULL
  #	expand formulae
  if (any(grep('formula', names(extras)))) {
    for (n in paste0(letters[1:3], '.formula')) {
      if (!is.null(extras[[n]]) && !is.null(mcall[[n]]))
        extras[[n]] <- update.formula(mcall[[n]], extras[[n]])
    }
  }
  # update args
  mcall[names(extras)] <- extras
  #	add start arg if none of these args specified
  if (!sum(pmatch(names(extras), c("x", "y", "id", "fixed", "random", "a.formula", "b.formula",
                                   "c.formula", "start", "returndata"), 0))) {
    start. <- list(fixed=fixef(object), random=ranef(object))
    # update start if any of these args specified
    if (sum(pmatch(names(extras), c('data', 'subset', 'df', 'knots', 'bounds', 'xoffset', 'bstart'), 0))) {
      # get data etc
      data <- eval(mcall$data)
      subset <- eval(mcall$subset, data)
      if (!is.null(subset)) data <- data[subset, ]
      x <- eval(mcall$x, data)
      xoffset <- object$xoffset
      if (is.null(xoffset)) xoffset <- mean(x)
      x <- x - xoffset
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
            newre <- matrix(0, nrow=sum(newid), ncol=dim(ranef(object))[2],
                            dimnames=list(levels(id)[newid], dimnames(ranef(object))[[2]]))
            start.$random <- rbind(start.$random, newre)
            cat(sum(newid), 'subjects added\n')
          }
        }
      }
      #	update fixed effects
      if (length(fixef(object)) > df + 1) fixed.extra <- (df+2):length(fixef(object))
      else fixed.extra <- NULL
      # new arg xoffset
      if (!is.null(extras$xoffset)) {
        xoffset.t <- xoffset
        xoffset <- eval(extras$xoffset)
        xoffset.t <- xoffset - xoffset.t
        x <- x - xoffset.t
        knots <- knots - xoffset.t
        bounds <- bounds - xoffset.t
      }
      # new arg knots
      if (!is.null(extras$knots)) {
        knots <- eval(extras$knots) - xoffset
        df <- length(knots) + 1
        mcall$df <- NULL
      }
      # new arg df
      else if (!is.null(extras$df)) {
        df <- eval(extras$df)
        knots <- quantile(x, (1:(df-1))/df)
        mcall$knots <- NULL
      }
      # new arg bounds
      if (!is.null(extras$bounds)) {
        bounds <- eval(extras$bounds)
        if (length(bounds) == 1) bounds <- range(x) + abs(bounds) * c(-1,1) * diff(range(x))
        else bounds <- bounds - xoffset
      }
      #	get spline start values
      spline.lm <- lm(predict(object, data, level=0) ~ ns(x, knots=knots, Bound=bounds))
      start.$fixed <- c(coef(spline.lm)[c(2:(df+1), 1)], start.$fixed[fixed.extra])
      # new arg bstart
      if (!is.null(extras$bstart) && !is.null(start.$fixed['b'])) {
        bstart <- eval(extras$bstart)
        if (bstart == 'mean') bstart <- mean(x)
        else bstart <- bstart - xoffset
        start.$fixed['b'] <- bstart
      }
    }
    #	save start. object
    assign('start.', start., parent.frame())
    mcall[['start']] <- quote(start.)
  }
  if (evaluate)
    eval(mcall, parent.frame())
  else mcall
}
