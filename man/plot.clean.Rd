\name{plot.clean}
\alias{plot.clean}
\title{plot.clean}
\description{
A version of \code{mplot} to plot growth curves and identify outliers. Outliers are clicked on and the corresponding growth curve is highlighted.}
\usage{
plot.clean(x, y = NULL, id = NULL, data = NULL, n = length(x), pch = 20, ...)
}
\arguments{
  \item{x}{
vector of x coordinates.}
  \item{y}{
vector of y coordinates.}
  \item{id}{
factor of subject levels indexing each growth curve.}
  \item{data}{
optional dataframe containing \code{x}, \code{y} and \code{id}.}
  \item{n}{
maximum number of points to be identified.}
  \item{pch}{
plot character.}
  \item{\dots}{
Further graphical parameters (see \code{\link{par}}) may also be supplied as arguments, particularly line type, lty, line width, lwd and color, col.}
}
\value{
\code{plot.clean} returns either a vector \code{rows} or a list:
\item{rows}{a vector of row numbers corresponding to the selected points.}
\item{data}{a dataframe returned if \code{data} has been provided in the call, containing as columns \code{rows}, \code{id}, \code{x} and \code{y} for the selected points.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
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
