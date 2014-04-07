\name{bupdate}
\alias{bupdate}
\title{Update the b fixed effect to minimise the b-c random effect correlation}
\description{
A function to update the value of \code{bstart}, the starting value for the b fixed effect,
 to minimise the correlation between the random effects b and c.}
\usage{bupdate(x)}
\arguments{
  \item{x}{a \code{sitar} object.}
}
\value{
Returns an updated value of the b fixed effect, based on the random effect covariance matrix.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
## fit sitar model with b fixed effect starting value defaulting to 'mean'
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)
print(fixef(m1)['b'])

## refit with starting value chosen to minimise b-c correlation
m2 <- update(m1, bstart=bupdate(m1))
print(fixef(m2)['b'])
}
\keyword{regression}
