\name{zLMS}
\alias{zLMS}
\title{zLMS}
\description{Routines to handle references constructed with the
LMS method. Given a set of LMS values, the functions convert
z-scores to measurement centiles and vice versa.}
\usage{
zLMS(x, L, M, S, data = NULL)
}
\arguments{
  \item{x}{
x}
  \item{L}{
L}
  \item{M}{
L}
  \item{S}{
S}
  \item{data}{
data}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, L, M, S, data = NULL) 
{
    attach(data)
    on.exit(detach(data))
    ((x/M)^L - 1)/L/S
  }
}
