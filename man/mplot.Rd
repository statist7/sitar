\name{mplot}
\alias{mplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{mplot}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mplot(x, y, id, data = NULL, subset = NULL, newplot = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{newplot}{
%%     ~~Describe \code{newplot} here~~
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
function (x, y, id, data = NULL, subset = NULL, newplot = TRUE, 
    ...) 
{
    xl <- deparse(substitute(x))
    yl <- deparse(substitute(y))
    if (!is.null(data)) {
        data <- data[, c(xl, yl, deparse(substitute(id)))]
        on.exit(detach(data))
        attach(data)
        xv <- substitute(x)
        yv <- substitute(y)
    }
    else {
        xv <- substitute(x, globalenv())
        yv <- substitute(y, globalenv())
        data <- as.data.frame(cbind(x, y, id))
    }
    ARG <- list(...)
    cnames <- names(ARG)
    if (!is.null(cnames)) {
        cnames <- cnames[lapply(ARG, length) == dim(data)[[1]]]
        data[, cnames] <- ARG[cnames]
        ARG[cnames] <- NULL
    }
    if (is.null(subset)) 
        subset <- rep(TRUE, length(id))
    subset <- ifelse(is.na(x) | is.na(y), FALSE, subset)
    if (sum(subset, na.rm = TRUE) == 0) 
        stop("no data to plot")
    if (newplot) 
        do.call...("plot", list(x = eval(x)[subset], y = eval(y)[subset], 
            type = "n", xlab = xl, ylab = yl), ...)
    tt <- by(data[subset, ], id[subset], function(z) {
        zs <- z[order(eval(xv, z)), ]
        xvt <- eval(xv, zs)
        yvt <- eval(yv, zs)
        if (!is.null(cnames)) 
            ARG[cnames] <- as.list(as.data.frame(zs[, cnames]))
        do.call("lines", c(list(x = xvt, y = yvt), ARG))
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
