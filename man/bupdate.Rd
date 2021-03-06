% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sitarlib.R
\name{bupdate}
\alias{bupdate}
\title{Update the b fixed effect to minimise the b-c random effect correlation}
\usage{
bupdate(x)
}
\arguments{
\item{x}{a \code{sitar} object.}
}
\value{
Returns an updated value of the b fixed effect, based on the random
effect covariance matrix.
}
\description{
A function to update the value of \code{bstart}, the starting value for the
b fixed effect, to minimise the correlation between the random effects b and
c.
}
\examples{

## fit sitar model with b fixed effect starting value defaulting to 'mean'
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)
print(fixef(m1)['b'])

## refit with starting value chosen to minimise b-c correlation and df increased
m2 <- update(m1, bstart=bupdate(m1), df=6)
print(fixef(m2)['b'])

}
\author{
Tim Cole \email{tim.cole@ucl.ac.uk}
}
\keyword{regression}
