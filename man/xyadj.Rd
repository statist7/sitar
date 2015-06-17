\name{xyadj}
\alias{xyadj}
\title{Adjust x and y variables for SITAR random effects}
\description{
\code{xyadj} takes raw \code{x} and \code{y} values for a SITAR model
and returns them adjusted for the subject-specific random effects.
}
\usage{
xyadj(object, data)
}
\arguments{
  \item{object}{a SITAR model.}
  \item{data}{an optional data frame containing values for x and id,
and optionally y, the variable names matching those in \code{object}.}
}
\details{
The x and y values are adjusted using the formulae
  \deqn{(x - xoffset - b<fixed> - b<random>) * exp(c<random>) + xoffset + b<fixed>}
  \deqn{y - a<random>}
where missing values of \code{xoffset} or the fixed or random effects are set to zero.
}
\value{
The list
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
points(xyadj(m1), col='red', pch=19)
}
