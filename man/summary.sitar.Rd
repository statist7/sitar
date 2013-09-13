\name{summary.sitar}
\alias{summary.sitar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{summary.sitar}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.sitar(object, adjustSigma = TRUE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{adjustSigma}{
%%     ~~Describe \code{adjustSigma} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
