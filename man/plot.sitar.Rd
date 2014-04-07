\name{plot.sitar}
\alias{plot.sitar}
\alias{lines.sitar}
\title{Plot SITAR model}
\description{
\code{plot} and \code{lines} methods for objects of class \code{sitar}, providing various flavours of plot of the fitted growth curves. }
\usage{
\method{plot}{sitar}(x, opt = "dv", labels, apv = FALSE, xfun, yfun, subset = NULL, 
 abc = c(a = 0, b = 0, c = 0), add = FALSE, nlme = FALSE, ...)
\method{lines}{sitar}(x, ...)
}
\arguments{
  \item{x}{object of class \code{sitar}.}
  \item{opt}{
character string containing a subset of letters corresponding to the options: 'd' for fitted Distance curve, 'v' for fitted Velocity curve, 
'e' for fitted fixed Effects distance curve,
'f' for individual Fitted growth curves, 'u' for Unadjusted individual growth curves, and 'a' for Adjusted individual growth curves. Options 'dvef' give spline curves, while 'ua' give curves made up as line segments. If both 'd' and 'v' are specified, the axis for the velocity curve appears on the right side of the plot (y2), and a legend identifying the distance and velocity curves is provided.}
  \item{labels}{
optional character vector containing plot labels for \code{x}, \code{y} and \code{y} velocity from the original SITAR model. The first two elements can alternatively be provided via \code{\link{par}} parameters \code{xlab} and \code{ylab}, and the third element via \code{y2par} (see Details). Default labels are the names of \code{x} and \code{y}, and "\code{y} velocity".}
  \item{apv}{
optional logical specifying whether or not to calculate the age at peak velocity from the velocity curve. If TRUE, age at peak velocity is calculated as the age when the second derivative of the fitted curve changes sign (after applying \code{xfun} and/or \code{yfun} if specified). Age at peak velocity is marked in the plot with a vertical dotted line, and its value, along with peak velocity, is printed and returned.}
  \item{xfun}{
optional function to be applied to \code{x} prior to plotting.}
  \item{yfun}{
optional function to be applied to \code{y} prior to plotting.}
  \item{subset}{
optional logical vector of length \code{x} defining a subset of \code{data} rows to be plotted.}
  \item{abc}{
character string defining an \code{id} level, or a list of values labelled a, b and c, e.g. abc=list(a=1, b=0.1, c=-0.1), where the values of a, b and c are random effects defining an individual growth curve. If a, b or c are missing their values are treated as zero.}
  \item{add}{
optional logical defining whether the plot is pre-existing (TRUE) or new (FALSE). TRUE is equivalent to using \code{lines}.}
  \item{nlme}{
optional logical which set TRUE plots the model as an \code{nlme} object, using \code{plot.nlme} arguments.}
  \item{\dots}{
Further graphical parameters (see \code{par}) may also be supplied as arguments, e.g. axis labels \code{xlab} and \code{ylab}, line type \code{lty}, line width \code{lwd}, and color \code{col}. For the velocity (y2) plot \code{y2par} can be used (see Details).}
}
\details{
For option 'dv' (the default) the velocity curve plot (with right axis) can be annotated with \code{par} parameters given as a named list called \code{y2par}. To suppress the
legend that comes with it set \code{xlegend = NULL}.
}
\value{
Returns a list of three objects:
\item{ss}{\code{smooth.spline} object corresponding to the fitted distance curve.}
\item{usr}{value of \code{par('usr')} for the main plot.}
\item{usr2}{the value of \code{par('usr')} for the velocity (y2) plot.}
\item{apv}{if argument \code{apv} is TRUE a named list giving the age at peak velocity (apv) and peak velocity (pv) from the fitted velocity curve.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{
\code{\link{lines.sitar}}, \code{\link{mplot}}, \code{\link{plotclean}}, \code{\link{y2plot}}
}
\examples{
## fit sitar model
m1 <- sitar(x=age, y=height, id=id, data=heights, df=7)

## draw fitted distance and velocity curves
## with velocity curve in blue
## adding age at peak velocity
plot(m1, y2par=list(col='blue'), apv=TRUE)

## draw individually coloured growth curves adjusted for random effects
## using same x-axis limits as for previous plot
plot(m1, opt='a', col=1 + as.integer(id) \%\% 5, xlim=xaxsd())

## add mean curve in red
lines(m1, opt='d', col='red', lwd=2)
}
\keyword{aplot}
