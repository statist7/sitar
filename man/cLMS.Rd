\name{cLMS}
\alias{cLMS}
\alias{zLMS}
\title{LMS conversion to and from z-scores}
\description{Routines to handle references constructed with the
LMS method. Given a set of LMS values, the functions convert
z-scores to measurement centiles and vice versa. }
\usage{
zLMS(x, L, M, S, data = NULL)
cLMS(z, L, M, S, data = NULL)
}
\arguments{
  \item{x}{vector of measurements to be converted to z-scores.}
  \item{z}{scalar or vector of z-scores to be converted to measurements.}
  \item{L}{vector of Box-Cox transformation lambda values, L in the LMS method.}
  \item{M}{vector of medians, M in the LMS method.}
  \item{S}{vector of coefficients of variation, S in the LMS method.}
  \item{data}{optional data frame containing the variables L, M and S.}
}
\details{The variables x, L, M and S should all be the same length. Variable \code{z}
can be scalar or vector.}
\value{ If \code{z} is scalar \code{cLMS} returns a vector the same length as \code{L}, \code{M} and \code{S}.  
If \code{z} is a vector \code{cLMS} returns a matrix where dim(matrix) = \code{c(length(z),length(M))}.
\code{zLMS} returns z-scores the same length as \code{L}, \code{M} and \code{S}.}
\references{
%% ~put references to the literature/web site here ~
%% }
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{z2cent}}, \code{\link{lms2z}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
##
## The function is currently defined as
function (z, L, M, S, data = NULL) 
{
    attach(data)
    on.exit(detach(data))
    as.numeric(M * (1 + L * S \%o\% z)^(1/L))
}
## The function is currently defined as
function (x, L, M, S, data = NULL) 
{
    attach(data)
    on.exit(detach(data))
    ((x/M)^L - 1)/L/S
  }
}
