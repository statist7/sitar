\name{plot.clean}
\alias{plot.clean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.clean}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.clean(x, y = NULL, id = NULL, data = NULL, n = length(x), pch = 20, ...)
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
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
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
function (x, y = NULL, id = NULL, data = NULL, n = length(x), 
    pch = 20, ...) 
{
    if (!is.null(data)) {
        if (!deparse(substitute(data)) \%in\% search()) {
            on.exit(detach(data))
            attach(data)
        }
    }
    xt <- substitute(x)
    xlab <- deparse(xt)
    yt <- substitute(y)
    ylab <- deparse(yt)
    if (is.null(id)) {
        plot(x, y, xlab = xlab, ylab = ylab, pch = 46, ...)
        idlab <- NULL
    }
    else {
        idlab <- deparse(substitute(id))
        id <- factor(id)
        plot(x, y, type = "n", ...)
        tt <- by(data, id, function(z) lines(eval(yt, z) ~ eval(xt, 
            z), col = "gray", type = "o", pch = 46, ...))
    }
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    sel <- rep(FALSE, length(x))
    res <- integer(0)
    title("click on outliers in plot - then right-click to escape")
    while (sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, 
            ...)
        if (!length(ans)) 
            break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch)
        if (!is.null(id)) {
            text(x[ans], y[ans], label = id[ans], adj = c(0.5, 
                0.1))
            lines(x[id == id[ans]], y[id == id[ans]])
        }
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res <- res[order(res)]
    if (is.null(data)) 
        res
    else list(rows = res, data = data[res, c(idlab, xlab, ylab)])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
