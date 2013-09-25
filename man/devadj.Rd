\name{devadj}
\alias{devadj}
\alias{varexp}
\title{Two ways to compare SITAR models for fit.}
\description{
\code{devadj} calculates the deviance for SITAR models, adjusting 
the likelihood for Box-Cox transformed y variables. \code{varexp} 
calculates the variance explained by SITAR models, compared to 
the fixed effect model.
}
\usage{
devadj(..., pattern)
varexp(..., pattern)
}
\arguments{
  \item{\dots}{one or more SITAR models.}
  \item{pattern}{regular expression defining names of SITAR models.}
}
\details{
The deviance is adjusted if the y variable is power-transformed, using the formula 
        \deqn{adjusted deviance = deviance - 2n ( (\lambda-1) * log(gm) + log(abs(\lambda)) )}{%
        adjusted deviance = deviance - 2n ( (lambda-1) * log(gm) + log(abs(lambda)) )}
where lambda is the power transform, and n and gm are the length and geometric mean of y.

The variance explained is given by
        \deqn{% explained = 100 * (1 - (\sigma2/\sigma1)^2)}{%
        \% explained = 100 * (1 - (sigma2/sigma1)^2)}
where sigma1 is the cross-sectional RSD and sigma2 the SITAR RSD.

Objects not of class \code{sitar} are ignored by \code{pattern} in \code{varexp}. In \code{devadj} non-\code{sitar} objects containing \code{call} and\code{logLik} are included.
}
\value{
For \code{devadj} a vector of deviances in increasing order. For \code{varexp} a vector of percentages in decreasing order. Value names correspond to the models.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
data(heights)
## fit sitar model for height
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

## update it for log(height) 
m2 <- update(m1, y=sqrt(height))

## compared variance explained in the two models
varexp(m1, m2)

## compare their deviances, adjusting for sqrt transform
## the regular expression looks for model names starting "m"
## followed by a digit
devadj(pattern="^m[0-9]")
}
\keyword{regression}