\name{uk90}
\alias{uk90}
\docType{data}
\title{UK 1990 growth reference}
\description{
The UK 1990 growth reference (Freeman et al 1995, Cole et al 1998) for height, weight, body mass index, circumferences and percent body fat, fitted by the LMS method and summarised by values of L, M and S by sex from 23 weeks gestation to 23 years.
}
\usage{data(uk90)}
\format{
  A data frame with 588 observations on the following 26 variables:
  \describe{
    \item{\code{years}}{numeric vector}
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
    \item{\code{L.sitht}}{numeric vector}
    \item{\code{M.sitht}}{numeric vector}
    \item{\code{S.sitht}}{numeric vector}
    \item{\code{L.leglen}}{numeric vector}
    \item{\code{M.leglen}}{numeric vector}
    \item{\code{S.leglen}}{numeric vector}
    \item{\code{L.waist}}{numeric vector}
    \item{\code{M.waist}}{numeric vector}
    \item{\code{S.waist}}{numeric vector}
    \item{\code{L.bfat}}{numeric vector}
    \item{\code{M.bfat}}{numeric vector}
    \item{\code{S.bfat}}{numeric vector}
    \item{\code{sex}}{two-level factor with level 1 male and level 2 female}
  }
}
\details{
The L, M and S values for each measurement correspond respectively to the Box-Cox power, median and coefficient of variation of the distribution by age and sex (Cole & Green 1992). The short names for each measurement (see \code{\link{lms2z}}) are as follows: height (ht), weight (wt), body mass index (bmi), head circumference (head), sitting height (sitht), leg length (leglen), waist circumference (waist) and percent body fat (fat).
}
\source{
The values are tabulated in the Excel spreadsheet British1990.xls provided with the Excel add-in LMSgrowth from \url{http://www.healthforallchildren.com/?product_cat=software}.
}
\references{
Cole TJ, Green PJ. Smoothing reference centile curves: the LMS method and penalized likelihood. Stat Med 1992;11:1305-19.

Cole TJ, Freeman JV, Preece MA. British 1990 growth reference centiles for weight, height, body mass index and head circumference fitted by maximum penalized likelihood. Stat Med 1998;17:407-29.

Freeman JV, Cole TJ, Chinn S, et al. Cross sectional stature and weight reference curves for the UK, 1990. Arch Dis Child 1995;73:17-24.
}
\examples{
data(uk90)
## calculate median BMI in girls from birth to 10 years
lms2z(x = 0:10, y = 0, sex = 2, data = NULL, measure = 'bmi', ref = 'uk90', toz = FALSE)
}
\keyword{datasets}
