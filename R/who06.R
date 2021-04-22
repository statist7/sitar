#' The WHO 2006 growth standard
#'
#' The WHO growth standard (WHO 2006) for height, weight, body mass index,
#' circumferences and skinfold thicknesses, fitted by the LMS method and
#' summarised by values of L, M and S by sex from birth to 5 years.
#'
#' The L, M and S values for each measurement correspond respectively to the
#' Box-Cox power, median and coefficient of variation of the distribution by
#' age and sex (Cole & Green 1992). The short names and units for each measurement (see
#' \code{\link{LMS2z}}) are as follows: height (ht, cm), weight (wt, kg), body mass
#' index (bmi, kg/m2), head circumference (head, cm), arm circumference (arm, cm), subscapular
#' skinfold (subscap, mm), and tricep skinfold (tricep, mm).
#'
#' @name who06
#' @docType data
#' @format A tibble with 150 observations on the following 23 variables:
#' \describe{ \item{years}{age from 0 to 5 years}
#' \item{L.ht}{numeric vector}
#' \item{M.ht}{numeric vector}
#' \item{S.ht}{numeric vector}
#' \item{L.wt}{numeric vector}
#' \item{M.wt}{numeric vector}
#' \item{S.wt}{numeric vector}
#' \item{L.bmi}{numeric vector}
#' \item{M.bmi}{numeric vector}
#' \item{S.bmi}{numeric vector}
#' \item{L.head}{numeric vector}
#' \item{M.head}{numeric vector}
#' \item{S.head}{numeric vector}
#' \item{L.arm}{numeric vector}
#' \item{M.arm}{numeric vector}
#' \item{S.arm}{numeric vector}
#' \item{L.subscap}{numeric vector}
#' \item{M.subscap}{numeric vector}
#' \item{S.subscap}{numeric vector}
#' \item{L.tricep}{numeric vector}
#' \item{M.tricep}{numeric vector}
#' \item{S.tricep}{numeric vector}
#' \item{sex}{two-level factor with level 1 male and level 2 female} }
#' @references World Health Organization. WHO Child Growth Standards: Methods
#' and development: Length/height-for-age, weight-for-age, weight-for-length,
#' weight-for-height and body mass index-for-age. Geneva: WHO; 2006.
#'
#' Cole TJ, Green PJ. Smoothing reference centile curves: the LMS method and
#' penalized likelihood. Stat Med 1992;11:1305-19.
#' @source \url{https://www.who.int/toolkits/child-growth-standards}
#' @keywords datasets
#' @examples
#' data(who06)
#' ## calculate z-score for length 60 cm in boys at age 0:12 months
#' LMS2z(x = 0:12/12, y = 60, sex = 1, measure = 'ht', ref = 'who06')
"who06"
