\name{devadj}
\alias{devadj}
\title{
Compares the deviance of models including an adjustment for power-transformed y variables.}
\description{
\code{devadj} calculates the deviance for one or more SITAR models, adjusting the likelihood for Box-Cox transformed y variables.
}
\usage{
devadj(..., pattern)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{one or more SITAR model names}
  \item{pattern}{a character string defining a search term for names of SITAR models. It is passed to \code{ls(pattern=pattern)}, and  any objects so identified are concatenated with  explicitly given model names.}
}
\details{
The deviance for each model is adjusted if the y variable is power-transformed, using the formula 
        adjusted deviance = deviance - 2n {(lambda-1) * log(gm) + log(abs(lambda))}
where lambda is the power transform, and n and gm are the length and geometric mean of y.
}
\value{
A vector of deviances, in increasing order, each named by the corresponding model.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{varexp}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (..., pattern) 
{
	ARG <- list(...)
	if (!missing(pattern)) {
		pattern <- ls(env=parent.frame(), pattern=pattern)
		ARG <- c(ARG, pattern)
	}
	else pattern <- NULL
	dev <- numeric(length(ARG))
	names(dev) <- c(match.call(expand.dots=FALSE)$..., pattern)
	for (i in 1:length(ARG)) {
		obj <- ARG[[i]]
		if (is.character(obj)) obj <- get(obj)
		if (!is.null(obj$call.sitar)) obj$call <- obj$call.sitar
		if (!is.null(obj$call$y)) ycall <- obj$call$y else ycall <- obj$call$model[[2]]
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
		dev[i] <- -2 * (obj$logLik + (lambda-1) * sly + length(y) * log(abs(lambda) + (lambda == 0)))
	}
	dev[order(dev)]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
