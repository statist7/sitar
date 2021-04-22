#' The WHO 2006 growth standard and WHO 2007 growth reference
#'
#' The WHO growth standard (WHO 2006) and growth reference (2007) for height,
#' weight and body mass index, fitted by the LMS method and
#' summarised by values of L, M and S by sex from birth to 19 years.
#'
#' The L, M and S values for each measurement correspond respectively to the
#' Box-Cox power, median and coefficient of variation of the distribution by
#' age and sex (Cole & Green 1992). The short names and units for each measurement (see
#' \code{\link{LMS2z}}) are as follows: height (ht, cm), weight (wt, kg) and body mass
#' index (bmi, kg/m2).
#'
#' @name who0607
#' @docType data
#' @format A tibble with 486 observations on the following 11 variables:
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
#' \item{sex}{two-level factor with level 1 male and level 2 female} }
#' @references Cole TJ, Green PJ. Smoothing reference centile curves: the
#' LMS method and penalized likelihood. Stat Med 1992;11:1305-19.
#'
#' World Health Organization. WHO Child Growth Standards: Methods
#' and development: Length/height-for-age, weight-for-age, weight-for-length,
#' weight-for-height and body mass index-for-age. Geneva: WHO; 2006.
#'
#' de Onis M, Onyango AW, Borghi E, Siyam A, Nishida C, Siekmann J. Development
#' of a WHO growth reference for school-aged children and adolescents.
#' Bull WHO 2007;85:660-7.
#' @source \url{https://www.who.int/growthref/growthref_who_bull/en/}
#' @keywords datasets
#' @examples
#' data(who0607)
#' ## calculate 98th centile for BMI in girls from birth to 19 years
#' round(
#'   setNames(
#'     LMS2z(x = 0:19, y = 2, sex = 2, measure = 'bmi', ref = 'who0607',
#'       toz = FALSE), 0:19), 1)
"who0607"
