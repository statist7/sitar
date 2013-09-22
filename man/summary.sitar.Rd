\name{summary.sitar}
\alias{summary.sitar}
\title{summary.sitar}
\description{
A \code{summary} method for \code{sitar} objects based on \code{\link{summary.lme}}.}
\usage{
\method{summary}{sitar}(object, adjustSigma = TRUE, verbose = FALSE, ...)
}
\arguments{
  \item{object}{an object inheriting from class \code{sitar}.}
  \item{adjustSigma}{an optional logical value (see \code{\link{summary.lme}}.}
  \item{verbose}{an optional logical value used to control the amount of output in the \code{print.summary.sitar} method.}
  \item{\dots}{some methods for this generic require additional arguments. None are used in this method.}
}
\value{
an object inheriting from class \code{summary.sitar} with all components included in \code{object} (see \code{\link{lmeObject}} for a full description of the components) plus the components detailed for \code{\link{summary.lme}} and the following components:
\item{x.adj}{vector of length \code{x} in \code{object} with \code{x} values adjusted for subject-specific random effects b and c.}
\item{y.adj}{vector of length \code{y} in \code{object} with \code{y} values adjusted for subject-specific random effects a.}
\item{apv}{length 2 vector giving respectively age at peak velocity and peak velocity based on the fitted distance curve (using transformed \code{x} and \code{y} where specified).}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, adjustSigma = TRUE, verbose = FALSE, ...) 
{
    class(object) <- class(object)[-1]
    object <- summary(object, adjustSigma = adjustSigma, verbose = verbose, 
        ...)
    random <- as.character(object$call$random)[2]
    mcall <- object$call.sitar
    data <- eval(mcall$data)
    attach(data, warn.conflicts = FALSE)
    on.exit(detach(data))
    y <- eval(mcall$y)
    nf <- length(fitted(object))
    if (nf != length(y)) 
        stop(paste("model (length=", nf, ") incompatible with data (rows=", 
            length(y), ")", sep = ""))
    x <- eval(mcall$x)
    id <- eval(mcall$id)
    xoffset <- object$xoffset
    if (is.null(xoffset)) 
        xoffset <- mean(x)
    x.adj <- x - xoffset
    if (grepl("b", random)) 
        x.adj <- x.adj - ranef(object)[factor(id), "b"]
    if (grepl("c", random)) 
        x.adj <- x.adj * exp(ranef(object)[factor(id), "c"])
    y.adj <- y
    if (grepl("a", random)) 
        y.adj <- y.adj - ranef(object)[factor(id), "a"]
    object$y.adj <- y.adj
    object$x.adj <- x.adj + xoffset
    object$apv <- makess(x, fitted(object, level = 0))$apv
    class(object) <- c("summary.sitar", class(object))
    object
  }
}
