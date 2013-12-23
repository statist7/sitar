\name{plotclean}
\alias{plotclean}
\title{plotclean}
\description{
A version of \code{mplot} to plot growth curves and identify outliers. When outliers are
 clicked on the corresponding growth curve is highlighted.}
\usage{
plotclean(x, y = NULL, id = NULL, data = parent.frame(), n = length(x), pch = 20, ...)
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
Further graphical parameters (see \code{\link{par}}) may also be supplied as arguments,
 particularly line type, lty, line width, lwd and color, col.}
}
\value{
\code{plotclean} returns either a vector \code{rows} or a list:
\item{rows}{a vector of row numbers corresponding to the selected points.}
\item{data}{a dataframe returned if \code{data} has been provided in the call, containing
 as columns \code{rows}, \code{id}, \code{x} and \code{y} for the selected points.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
\donttest{plotclean(age, height, id, heights)}
}
