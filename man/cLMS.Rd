\name{cLMS}
\alias{cLMS}
\alias{zLMS}
\title{LMS conversion to and from z-scores}
\description{Routines to handle references constructed with the
LMS method. Given a set of LMS values, the functions convert
z-scores to measurement centiles and vice versa.}
\usage{
zLMS(x, L, M, S, data = NULL)
cLMS(z, L, M, S, data = NULL)
}
\arguments{
  \item{x}{vector of measurements to be converted to z-scores.}
  \item{z}{scalar or vector of z-scores to be converted to measurements.}
  \item{L}{vector of Box-Cox transformation (lambda) values, L in the LMS method.}
  \item{M}{vector of medians (mu), M in the LMS method.}
  \item{S}{vector of coefficients of variation (sigma), S in the LMS method.}
  \item{data}{optional data frame containing the variables L, M and S.}
}
\details{x, L, M and S should all be the same length. \code{z}
can be scalar or vector. The formulae converting \code{x} to \code{z} and vice versa are:
\deqn{z = \frac{(x/M)^L - 1}{L S}}{%z = ((x/M)^L - 1)/L/S}
\deqn{x = M (1 + L S z)^{1/L})}{%
x = M (1 + L S z)^(1/L)}
where L has 10^-7 added if it is zero.}
\value{\code{zLMS} returns a vector the same length as \code{M}.
If \code{z} is scalar \code{cLMS} returns a vector the same length as \code{M}.  
If \code{z} is a vector \code{cLMS} returns a matrix with the number of columns equal to its length.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{z2cent}}, \code{\link{lms2z}}}
\examples{
cLMS(z = -2:2, L = 1:-1, M = 5:7, S = rep(0.1, 3))
zLMS(x = 6.5, L = 1:-1, M = 5:7, S = rep(0.1, 3))
}
\keyword{arith}
