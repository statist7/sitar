\name{y2plot}
\alias{y2plot}
\title{Plot with two y axes}
\description{
Function to plot two y variables, y1 and y2, against a single x variable, with the y1
 and y2 axes on the left and right of the plot respectively.}
\usage{
y2plot(x, y1, y2 = NULL, labels, y2par = list(), add = FALSE,
  xy = NULL, xlegend = "topleft", inset = 0.04, ...)
}
\arguments{
  \item{x}{vector of ages.}
  \item{y1}{vector of measurements for plotting on left y axis.}
  \item{y2}{optional vector of measurements for plotting on right y axis.}
  \item{labels}{character vector containing labels for x, y1 and y2.}
  \item{y2par}{optional named list of par arguments to format the y2 axis.}
  \item{add}{logical flag to specify if a new plot (with axes etc) is to be drawn (FALSE)
   or an existing plot is to be added to (TRUE).}
  \item{xy}{optional list to pass \code{usr} and \code{usr2} (see Value).}
  \item{xlegend}{position for legend.}
  \item{inset}{inset for legend.}
  \item{\dots}{optional \code{par} arguments.}
}
\details{
y2plot draws up to two superimposed plots, one with the y axis on the left and the other
 on the right, with suitable adjustment for \code{par('mar')} and including a legend.
 The format for y1 is controlled by par arguments, and that for y2 by the list y2par.}
\value{
Returns the list
  \item{usr}{par('usr') for y1 axis}
  \item{usr2}{par('usr') for y2 axis}
In addition the variable \code{.par.usr2}, equal to \code{usr2}, is created in globalenv().
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{plot.sitar}}}
\examples{
## plot boys median height and weight on the UK 1990 reference
with(uk90[uk90$sex == 1,],
  y2plot(x=years, y1=M.ht, y2=M.wt, y2par=list(col='red'))
)
}
