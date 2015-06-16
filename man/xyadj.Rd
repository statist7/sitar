\name{xyadj}
\alias{xyadj}
\title{Adjust x and y variables for SITAR random effects}
\description{
\code{xyadj} takes raw \code{x} and \code{y} values for a SITAR model
and returns them adjusted for the subject-specific random effects.
}
\usage{
xyadj(object)
}
\arguments{
  \item{object}{a SITAR model.}
}
\details{
The x and y values are adjusted using the formulae
  \deqn{(x - xoffset - b<fixed> - b<random>) * exp(c<random>) + xoffset + b<fixed>}
  \deqn{y - a<random>}
where missing values of \code{xoffset} or the fixed or random effects are treated as zero.
}
\value{
The list
  \item{x}{numeric vector.}
  \item{y}{numeric vector the same length as x.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
data(heights)
## fit sitar model for height
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

## plot unadjusted data as growth curves
plot(m1, opt='u')

## overplot with adjusted data as points
points(xyadj(m1), col='red', pch=19)
}
