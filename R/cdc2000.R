#' The CDC 2000 growth reference
#'
#' The CDC growth reference (Kuczmarski et al 2000) for height,
#' weight, body mass index and head circumference, fitted by the LMS method and
#' summarised by values of L, M and S by sex from birth to 19 years.
#'
#' BMI starts at 2 years, and head circumference stops at 3 years.
#'
#' The L, M and S values for each measurement correspond respectively to the
#' Box-Cox power, median and coefficient of variation of the distribution by
#' age and sex (Cole & Green 1992). The short names and units for each measurement
#' (see \code{\link{LMS2z}}) are as follows: height (ht, cm), weight (wt, kg),
#' body mass index (bmi, kg/m2), head circumference (hc, cm).
#'
#' @name cdc2000
#' @docType data
#' @format A tibble with 484 observations on the following 14 variables:
#' \describe{ \item{years}{age from 0 to 19 years}
#' \item{L.ht}{numeric vector}
#' \item{M.ht}{numeric vector}
#' \item{S.ht}{numeric vector}
#' \item{L.wt}{numeric vector}
#' \item{M.wt}{numeric vector}
#' \item{S.wt}{numeric vector}
#' \item{L.bmi}{numeric vector}
#' \item{M.bmi}{numeric vector}
#' \item{S.bmi}{numeric vector}
#' \item{L.hc}{numeric vector}
#' \item{M.hc}{numeric vector}
#' \item{S.hc}{numeric vector}
#' \item{sex}{two-level factor with level 1 male and level 2 female} }
#' @references Cole TJ, Green PJ. Smoothing reference centile curves: the
#' LMS method and penalized likelihood. Stat Med 1992;11:1305-19.
#'
#' Kuczmarski RJ, Ogden CL, Guo SS, Grummer-Strawn LM, Flegal KM, Mei Z, Wei R,
#' Curtin LR, Roche AF, Johnson CL. 2000 CDC growth charts for the United States:
#' methods and development. Vital Health Stat, 2002, 11, 246, 1-190.

#' @keywords datasets
#' @examples
#' data(cdc2000)
#' ## calculate 98th centile for weight in girls from birth to 19 years
#' round(
#'   setNames(
#'     LMS2z(x = 0:19, y = 2, sex = 2, measure = 'wt', ref = 'cdc2000',
#'       toz = FALSE), 0:19), 1)
"cdc2000"
