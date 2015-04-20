\name{plotclean}
\alias{plotclean}
\title{Plot multiple growth curves to identify outliers}
\description{
A version of \code{mplot} to plot growth curves and identify outliers. When outliers are
 clicked on, and if id is specified, the corresponding growth curve is highlighted.
 If id is not specified the selected point is highlighted. Use right-click to exit.}
\usage{
plotclean(x, y, id = NULL, data = parent.frame(), n = length(x),
par.out = list(pch = 20), ...)
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
  \item{par.out}{
list of optional graphical parameters to control appearance of selected outlying points and lines.}
  \item{\dots}{
Further graphical parameters (see \code{\link{par}}) may also be supplied as arguments
for lines and points, particularly line type, lty, line width, lwd and color, col.}
}
\value{
\code{plotclean} returns either a vector \code{rows} (if data is not specified) or a list:
\item{rows}{a vector of row numbers corresponding to the selected points.}
\item{data}{a subset of \code{data} consisting of rows \code{rows},
and columns \code{id}, \code{x} and \code{y}.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
if (interactive()) plotclean(age, height, id, heights)
}
