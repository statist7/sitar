#' UK-WHO growth reference omitting preterm data
#'
#' The UK-WHO growth reference for height, weight, BMI and head circumference
#' (see Wright et al 2010), fitted by the LMS method and summarised by values of
#' L, M and S by sex and postnatal age from term birth (see Details) to 20 years.
#'
#' The growth reference combines term birth data from the British 1990 growth
#' reference (Cole et al 2011), the WHO growth standard from 2 postnatal weeks
#' to 4 years, and the British 1990 reference from 4 to 20 years.
#'
#' Age is measured in years, and term birth corresponds to ages between 37 and
#' 42 weeks gestation, where 40 weeks gestation is 0 years. The conversion is:
#' \code{years = (weeks - 40) * 7 / 365.25}.
#'
#' The L, M and S values for each measurement correspond respectively to the
#' Box-Cox power, median and coefficient of variation of the distribution by
#' age and sex (Cole & Green 1992). The measurement short names and units (see
#' \code{\link{LMS2z}}) are as follows: height (ht, cm), weight (wt, kg),
#' BMI (bmi, kg/m2) and head circumference (head, cm).
#'
#' @name ukwhoterm
#' @docType data
#' @format A tibble with 512 observations on the following 15 variables:
#' \describe{
#' \item{years}{numeric vector - postnatal age in years}
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
#' \item{origin}{two-level factor indicating the provenance of the data, with
#' levels British1990 and WHO2006}
#' \item{sex}{two-level factor with level 1 male and level 2 female} }
#' @references Cole TJ, Green PJ. Smoothing reference centile curves: the LMS
#' method and penalized likelihood. Stat Med 1992;11:1305-19.
#'
#' Cole TJ, Williams AF, Wright CM, et al. Revised birth centiles for weight,
#' length and head circumference in the UK-WHO growth charts. Ann Hum Biol
#' 2011;38:7-11.
#'
#' Wright CM, Williams AF, Elliman D, et al. Using the new UK-WHO growth
#' charts. BMJ 2010;340:c1140.
#' @source The values are tabulated in the Excel spreadsheet UK_WHO_preterm.xls
#' provided with the Excel add-in LMSgrowth from
#' https://www.healthforallchildren.com/shop-base/software/lmsgrowth/
#' @keywords datasets
#' @examples
#' data(ukwhoterm)
#' ## calculate median weight (kg) in girls from 0 to 10 years
#' v <- LMS2z(x = 0:10, y = 0, sex = 2, measure = 'wt',
#'   ref = 'ukwhoterm', toz = FALSE)
#' setNames(v, 0:10)
"ukwhoterm"
