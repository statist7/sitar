#' Ways to compare SITAR models for fit
#'
#' \code{BICadj} and \code{AICadj} calculate the BIC and AIC for SITAR models,
#' adjusting the likelihood for Box-Cox transformed y variables. \code{varexp}
#' calculates the variance explained by SITAR models, compared to the
#' corresponding fixed effect models. \code{getL} is used by \code{[AB]ICadj} to
#' find the power that the y variable is raised to.
#'
#' The deviance is adjusted if the y variable is power-transformed, using the
#' formula
#' \deqn{adjusted deviance = deviance - 2n ( (\lambda-1) * log(gm) + %
#' log(abs(\lambda)) )}{%
#' deviance - 2n ( (lambda-1) * log(gm) + log(abs(lambda)) )}
#' where \eqn{\lambda}{lambda} is the power transform, and \eqn{n} and
#' \eqn{gm} are the length and geometric mean of \code{y}.
#'
#' The variance explained is given by \deqn{\% explained = 100 * (1 -%
#' (\sigma_2/\sigma_1)^2)}{\% explained = 100 * (1 - (sigma2/sigma1)^2)} where
#' \eqn{\sigma_1}{sigma1} is the fixed effects RSD and \eqn{\sigma_2}{sigma2}
#' the SITAR random effects RSD.
#'
#' \code{BICadj} and \code{AICadj} accept non-\code{sitar} models with a
#' \code{logLik} class. \code{varexp} ignores objects not of class
#' \code{sitar}.
#'
#' @aliases BICadj AICadj varexp getL
#' @param \dots one or more SITAR models.
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is
#' the classical AIC.
#' @param pattern regular expression defining names of models.
#' @param expr quoted or unquoted expression containing a single variable name.
#' @param data data frame to evaluate \code{expr}.
#' @return For \code{BICadj} and \code{AICadj} a named vector of deviances in
#' increasing order.  For \code{varexp} a named vector of percentages in
#' decreasing order. For \code{getL} the power the variable in \code{expr}
#' is raised to.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{BIC}}, \code{\link{AIC}}
#' @keywords regression
#' @examples
#' data(heights)
#' ## fit sitar model for height
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)
#'
#' ## update it for log(height)
#' m2 <- update(m1, y=sqrt(height))
#'
#' ## compare variance explained in the two models
#' varexp(m1, m2)
#'
#' ## compare BIC adjusting for sqrt transform
#' ## the pattern matches names starting with "m" followed by a digit
#' BICadj(pattern="^m[0-9]")
#'
#' ## find which power height is raised to
#' getL(quote(sqrt(sqrt(height))), heights)
#' @export BICadj
	BICadj <- function(..., pattern=NULL)
{	ARG <- match.call(expand.dots=FALSE)$...
	if (!is.null(pattern)) pattern <- ls(envir=parent.frame(), pattern=pattern)
	ARG <- unique(c(unlist(sapply(ARG, deparse)), pattern))
	dev <- sapply(ARG, function(obj) {
		obj <- eval(parse(text=obj))
		if (is.character(try(ll <- logLik(obj), TRUE))) return(NA)
#	check for call.sitar$y or else call$y or else call$formula
		if (!is.null(obj$call.sitar)) obj$call <- obj$call.sitar
#	check for call$y or call$model
		if (!is.null(obj$call$y)) ycall <- obj$call$y
			else if (!is.null(obj$call$model)) ycall <- obj$call$model[[2]]
			  else if (!is.null(obj$call$formula)) ycall <- obj$call$formula[[2]]
			    else return(NA)
		data <- eval(obj$call$data)
		if (!is.null(obj$call$subset))
		  data <- data[eval(obj$call$subset, data), ]
		y <- eval(parse(text=all.vars(ycall)), data)
		lambda <- getL(ycall, data)
		sly <- ifelse(lambda == 1, 0, sum(log(y)))
		BIC(ll) - 2 * ((lambda - 1) * sly + length(y) * log(abs(lambda) + (lambda == 0)))
	})
	dev <- dev[!is.na(dev)]
	if (length(dev) == 0) return(invisible())
	round(dev[order(dev)], 1)
}

#' @rdname BICadj
#' @export
	AICadj <- function(..., k=2, pattern=NULL)
{	ARG <- match.call(expand.dots=FALSE)$...
	if (!is.null(pattern)) pattern <- ls(envir=parent.frame(), pattern=pattern)
	ARG <- unique(c(unlist(sapply(ARG, deparse)), pattern))
	dev <- sapply(ARG, function(obj) {
	  obj <- eval(parse(text=obj))
	  if (is.character(try(ll <- logLik(obj), TRUE))) return(NA)
#	check for call.sitar$y or else call$y or else call$formula
		if (!is.null(obj$call.sitar)) obj$call <- obj$call.sitar
#	check for call$y or call$model
		if (!is.null(obj$call$y)) ycall <- obj$call$y
			else if (!is.null(obj$call$model)) ycall <- obj$call$model[[2]]
			  else if (!is.null(obj$call$formula)) ycall <- obj$call$formula[[2]]
			    else return(NA)
	  data <- eval(obj$call$data)
	  if (!is.null(obj$call$subset))
	    data <- data[eval(obj$call$subset, data), ]
	  y <- eval(parse(text=all.vars(ycall)), data)
	  lambda <- getL(ycall, data)
		sly <- ifelse(lambda == 1, 0, sum(log(y)))
		AIC(ll, k=k) - 2 * ((lambda - 1) * sly + length(y) * log(abs(lambda) + (lambda == 0)))
	})
	dev <- dev[!is.na(dev)]
	if (length(dev) == 0) return(invisible())
	round(dev[order(dev)], 1)
}

#############################
#
#	varexp
#
#############################

#' @rdname BICadj
#' @export
	varexp <- function(..., pattern=NULL)
	  #	returns % of variance explained by sitar model(s)
{	ARG <- list(...)
	if (!is.null(pattern)) {
	  pattern <- ls(envir=parent.frame(), pattern=pattern)
	  ARG <- c(ARG, lapply(as.list(pattern), get))
	}
	pc <- lapply(ARG, function(obj) {
	  if (!'ns' %in% names(obj))
	    NA
	  else
	    100 * (1 - (obj$sigma / summary(obj$ns)$sigma)^2)
	})
	names(pc) <- c(match.call(expand.dots=FALSE)$..., pattern)
	pc <- unlist(pc[!is.na(pc)])
	if (length(pc) == 0)
	  return(invisible())
	round(pc[rev(order(pc))], 2)
}

#############################
#
#	getL
#
#############################

#' @rdname BICadj
#' @export
	getL <- function(expr, data=parent.frame()) {
	  y <- all.vars(expr)
	  if (length(y) != 1) {
	    warning('y involves more than one variable')
	    return(NA)
	  }
	  . <- eval(expr, data)
	  . <- (eval(D(expr, y), data) * eval(parse(text=y), data) / .)[. != 0]
	  if (isTRUE(all.equal(min(.), max(.))))
	    median(.)
	  else	if ('log' %in% all.names(expr))
	    0
	  else
	    NA
	}
