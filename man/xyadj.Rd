\name{xyadj}
\alias{xyadj}
\title{Adjust x and y variables for SITAR random effects}
\description{
\code{xyadj} Adjusts \code{x} and \code{y} values for subject-specific
random effects from a SITAR model.
}
\usage{
xyadj(x, y=NULL, id, object, abc=NULL, tomean=TRUE)
}
\arguments{
  \item{x}{a vector of x coordinates.}
  \item{y}{a vector of y coordinates (default NULL).}
  \item{id}{a factor denoting the subject levels corresponding to \code{x} and \code{y}.}
  \item{object}{a SITAR model.}
  \item{abc}{a data frame containing random effects for (a subset of) a, b and c (default ranef(object)).}
  \item{tomean}{a logical defining the direction of adjustment. TRUE (default) indicates that
individual curves are translated and rotated to match the mean curve, while FALSE
indicates the reverse, the mean curve being translated and rotated to match an individual curve.}
}
\details{
When tomean = TRUE the x and y values are adjusted to
  \deqn{(x - xoffset - b<fixed> - b<random>) * exp(c<random>) + xoffset + b<fixed>}
  \deqn{y - a<random>}
When tomean = FALSE they are adjusted to
  \deqn{(x - xoffset - b<fixed>) / exp(c<random>) + xoffset + b<fixed> + b<random>}
  \deqn{y + a<random>}
In each case missing values of \code{xoffset} or the fixed or random effects are treated as zero.

When \code{abc} consists of a single row, it is applied to all \code{x,y} pairs
irrespective of \code{id}.
}
\value{
The list of adjusted values:
  \item{x}{numeric vector.}
  \item{y}{numeric vector the same length as x, or NULL.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
data(heights)
## fit sitar model for height
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

## plot unadjusted data as growth curves
plot(m1, opt='u')

## overplot with adjusted data as points
with(heights, points(xyadj(age, height, id, m1), col='red', pch=19))
}
