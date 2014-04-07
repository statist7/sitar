\name{makess}
\alias{makess}
\title{Create cubic smoothing spline from data}
\description{
A wrapper to fit a cubic smoothing spline to the supplied data, allowing derivatives to be calculated.}
\usage{makess(x, y, xfun, yfun)}
\arguments{
  \item{x}{
vector giving the values of the predictor variable.}
  \item{y}{
vector of responses.}
  \item{xfun}{
optional function to be applied to \code{x} prior to fitting.}
  \item{yfun}{
optional function to be applied to \code{y} prior to fitting.}
}
\value{
An object of class \code{smooth.spline} with the extra component \code{apv}, a length-2 vector containing the age at peak velocity and peak velocity. If no peak is identified (based on the 2nd derivative changing sign) the values are for maximum rather than peak velocity.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
## create smooth.spline mean height curve
data(heights)
ss <- with(heights, makess(age, height))

## and plot it
plot(ss, type='l', xlab='age', ylab='height')

## age at peak velocity, and peak velocity
print(ss$apv)
}
