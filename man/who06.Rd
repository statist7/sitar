\name{who06}
\alias{who06}
\docType{data}
\title{The WHO 2006 growth standard}
\description{
The WHO growth standard (WHO 2006) for height, weight, body mass index, circumferences and skinfold thicknesses, fitted by the LMS method and summarised by values of L, M and S by sex from birth to 5 years.
}
\usage{data(who06)}
\format{
  A data frame with 150 observations on the following 23 variables:
  \describe{
    \item{\code{years}}{age from 0 to 5 years}
    \item{\code{L.ht}}{numeric vector}
    \item{\code{M.ht}}{numeric vector}
    \item{\code{S.ht}}{numeric vector}
    \item{\code{L.wt}}{numeric vector}
    \item{\code{M.wt}}{numeric vector}
    \item{\code{S.wt}}{numeric vector}
    \item{\code{L.bmi}}{numeric vector}
    \item{\code{M.bmi}}{numeric vector}
    \item{\code{S.bmi}}{numeric vector}
    \item{\code{L.head}}{numeric vector}
    \item{\code{M.head}}{numeric vector}
    \item{\code{S.head}}{numeric vector}
    \item{\code{L.arm}}{numeric vector}
    \item{\code{M.arm}}{numeric vector}
    \item{\code{S.arm}}{numeric vector}
    \item{\code{L.subscap}}{numeric vector}
    \item{\code{M.subscap}}{numeric vector}
    \item{\code{S.subscap}}{numeric vector}
    \item{\code{L.tricep}}{numeric vector}
    \item{\code{M.tricep}}{numeric vector}
    \item{\code{S.tricep}}{numeric vector}
    \item{\code{sex}}{two-level factor with level 1 male and level 2 female}
  }
}
\details{
The L, M and S values for each measurement correspond respectively to the Box-Cox power, median and coefficient of variation of the distribution by age and sex (Cole & Green 1992). The short names for each measurement (see \code{\link{lms2z}}) are as follows: height (ht), weight (wt), body mass index (bmi), head circumference (head), arm circumference (arm), subscapular skinfold (subscap), and tricep skinfold (tricep).
}
\source{
\url{http://www.who.int/childgrowth/en/}
}
\references{
World Health Organization. WHO Child Growth Standards: Methods and development: Length/height-for-age, weight-for-age, weight-for-length, weight-for-height and body mass index-for-age. Geneva: WHO; 2006.

Cole TJ, Green PJ. Smoothing reference centile curves: the LMS method and penalized likelihood. Stat Med 1992;11:1305-19.
}
\examples{
data(who06)
## calculate z-score for length 60 cm in boys at age 0:12 months
lms2z(x = 0:12/12, y = 60, sex = 1, data = NULL, measure = 'ht', ref = 'who06')
}
\keyword{datasets}
