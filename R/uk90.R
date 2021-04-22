#' UK 1990 growth reference
#'
#' The UK 1990 growth reference (Freeman et al 1995, Cole et al 1998) for
#' height, weight, body mass index, circumferences and percent body fat, fitted
#' by the LMS method and summarised by values of L, M and S by sex from 23
#' weeks gestation to 23 years.
#'
#' The L, M and S values for each measurement correspond respectively to the
#' Box-Cox power, median and coefficient of variation of the distribution by
#' age and sex (Cole & Green 1992). The short names and units for each measurement (see
#' \code{\link{LMS2z}}) are as follows: height (ht, cm), weight (wt, kg), body mass
#' index (bmi, kg/m2), head circumference (head, cm), sitting height (sitht, cm), leg length
#' (leglen, cm), waist circumference (waist, cm) and percent body fat (fat, %).
#'
#' @name uk90
#' @docType data
#' @format A tibble with 588 observations on the following 26 variables:
#' \describe{
#' \item{years}{numeric vector}
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
#' \item{L.sitht}{numeric vector}
#' \item{M.sitht}{numeric vector}
#' \item{S.sitht}{numeric vector}
#' \item{L.leglen}{numeric vector}
#' \item{M.leglen}{numeric vector}
#' \item{S.leglen}{numeric vector}
#' \item{L.waist}{numeric vector}
#' \item{M.waist}{numeric vector}
#' \item{S.waist}{numeric vector}
#' \item{L.bfat}{numeric vector}
#' \item{M.bfat}{numeric vector}
#' \item{S.bfat}{numeric vector}
#' \item{sex}{two-level factor with level 1 male and level 2 female} }
#' @references Cole TJ, Green PJ. Smoothing reference centile curves: the LMS
#' method and penalized likelihood. Stat Med 1992;11:1305-19.
#'
#' Cole TJ, Freeman JV, Preece MA. British 1990 growth reference centiles for
#' weight, height, body mass index and head circumference fitted by maximum
#' penalized likelihood. Stat Med 1998;17:407-29.
#'
#' Freeman JV, Cole TJ, Chinn S, et al. Cross sectional stature and weight
#' reference curves for the UK, 1990. Arch Dis Child 1995;73:17-24.
#' @source The values are tabulated in the spreadsheet British1990.xls provided
#' with the Excel add-in LMSgrowth from:
#'
#' \url{https://www.healthforallchildren.com/shop-base/software/lmsgrowth/}.
#' @keywords datasets
#' @examples
#' data(uk90)
#' ## calculate median BMI in girls from birth to 10 years
#' LMS2z(x = 0:10, y = 0, sex = 2, measure = 'bmi', ref = 'uk90', toz = FALSE)
"uk90"
