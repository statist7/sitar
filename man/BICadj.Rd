\name{BICadj}
\alias{BICadj}
\alias{AICadj}
\alias{varexp}
\title{Ways to compare SITAR models for fit}
\description{
\code{BICadj} and \code{AICadj} calculate the BIC and AIC for SITAR models, adjusting 
the likelihood for Box-Cox transformed y variables. \code{varexp} calculates the variance
explained by SITAR models, compared to the corresponding fixed effect models.
}
\usage{
BICadj(..., pattern)
AICadj(..., k = 2, pattern)
varexp(..., pattern)
}
\arguments{
  \item{\dots}{one or more SITAR models.}
  \item{k}{numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.}
  \item{pattern}{regular expression defining names of SITAR models.}
}
\details{
The deviance is adjusted if the y variable is power-transformed, using the formula 
        \deqn{adjusted deviance = deviance - 2n ( (\lambda-1) * log(gm) + log(abs(\lambda)) )}{%
        deviance - 2n ( (lambda-1) * log(gm) + log(abs(lambda)) )}
where lambda is the power transform, and n and gm are the length and geometric mean of y.

The variance explained is given by
        \deqn{\% explained = 100 * (1 - (\sigma2/\sigma1)^2)}{%
        \% explained = 100 * (1 - (sigma2/sigma1)^2)}
where sigma1 is the fixed effects RSD and sigma2 the SITAR random effects RSD.

\code{BICadj} and \code{AICadj} accept non-\code{sitar} models with a \code{logLik} class.
\code{varexp} ignores objects not of class \code{sitar}.
}
\value{
For \code{BICadj} and \code{AICadj} a named vector of deviances in increasing order. 
For \code{varexp} a named vector of percentages in decreasing order.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{BIC}}, \code{\link{AIC}}}
\examples{
data(heights)
## fit sitar model for height
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

## update it for log(height)
m2 <- update(m1, y=sqrt(height))

## compare variance explained in the two models
varexp(m1, m2)

## compare BIC adjusting for sqrt transform
## the pattern matches names starting with "m" followed by a digit
BICadj(pattern="^m[0-9]")
}
\keyword{regression}
