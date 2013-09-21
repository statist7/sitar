\name{bupdate}
\alias{bupdate}
\title{bupdate}
\description{
A function to update the value of \code{bstart} when fitting a \code{sitar} model to minimise the correlation between the random effects b and c.}
\usage{bupdate(x)}
\arguments{
  \item{x}{a \code{sitar} object.}
}
\value{
Returns an updated value of the b fixed effect, based on the random effect covariance matrix.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
    cov <- as.matrix(x$modelStruct[[1]][[1]])
    fixef(x)["b"] + cov[2, 3]/cov[3, 3]
  }
}
