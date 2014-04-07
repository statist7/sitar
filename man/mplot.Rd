\name{mplot}
\alias{mplot}
\title{Plot multiple growth curves}
\description{
Function to plot multiple growth curves indexed by subject id.}
\usage{
mplot(x, y, id, data = parent.frame(), subset = NULL, add = FALSE, ...)
}
\arguments{
  \item{x}{
vector of x coordinates.}
  \item{y}{
vector of y coordinates.}
  \item{id}{
factor denoting subject levels.}
  \item{data}{
optional dataframe containing \code{x}, \code{y} and \code{id}.}
  \item{subset}{
optional logical defining a subset of rows in \code{data}.}
  \item{add}{
optional logical defining whether the plot is pre-existing (TRUE) or new (FALSE).}
  \item{\dots}{
Further graphical parameters (see \code{\link{par}}) may also be supplied as arguments, 
  particularly line type, lty, line width, lwd, color, col and character pch.}
}
\details{
The arguments \code{x}, \code{y} and \code{id} can be given as character strings. The 
  \code{\link{par}} parameters can be functions of vector variables, e.g. to 
  colour curves separately by \code{id} use: col = 1 + as.integer(id) \%\% 6.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{y2plot}}}
\examples{
mplot(age, height, id, heights)
}
