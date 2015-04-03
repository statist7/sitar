\name{cLMS}
\alias{cLMS}
\alias{zLMS}
\title{LMS conversion to and from z-scores}
\description{Routines to handle references constructed with the
LMS method. Given a set of LMS values, the functions convert
z-scores to measurement centiles and vice versa.}
\usage{
zLMS(x, L = 1, M, S)
cLMS(z, L = 1, M, S)
}
\arguments{
  \item{x}{vector of measurements to be converted to z-scores.}
  \item{z}{scalar or vector of z-scores to be converted to measurements.}
  \item{L}{vector of Box-Cox transformation (lambda) values, L in the LMS method.}
  \item{M}{vector of medians (mu), M in the LMS method.}
  \item{S}{vector of coefficients of variation (sigma), S in the LMS method.}
}
\details{x, L, M and S should all be the same length. \code{z}
can be scalar or vector. The formulae converting \code{x} to \code{z} and vice versa are:
\deqn{z = \frac{(x/M)^L - 1}{L S}}{%
z = ((x/M)^L - 1)/L/S}

\deqn{x = M (1 + L S z)^{1/L})}{%
x = M (1 + L S z)^(1/L)}
where L has 10^-7 added if it is zero. The LMS method is the same as the BCCG family in
 the \code{gamlss} package.}
\value{\code{zLMS} and \code{cLMS} return a vector the same length as \code{M}.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{z2cent}}, \code{\link{LMS2z}}}
\examples{
cLMS(z = 0:2, L = 1:-1, M = 5:7, S = rep(0.1, 3))
zLMS(x = 6.5, L = 1:-1, M = 5:7, S = rep(0.1, 3))
}
\keyword{arith}
