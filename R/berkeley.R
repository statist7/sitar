#' The Berkeley Child Guidance Study
#'
#' The Berkeley Child Guidance Study dataset contains longitudinal anthropometry
#' data for 136 children from birth to 21 years.
#'
#' The data are for 66 boys and 70 girls from Berkeley, California born in 1928-29
#' of north European ancestry, and followed from birth to 21 years. Measurements were
#' at ages 0, 0.085, 0.25 to 2 (3-monthly), 2 to 8 (annually), and 8 to 21 (6-monthly) years.
#'
#' The children were measured for height, weight (undressed), stem length, biacromial diameter,
#' bi-iliac diameter, leg circumference, and dynamometric strength.
#' The data were provided as an appendix to the book by Tuddenham and Snyder (1954),
#' and a few transcription errors are corrected here.
#' A further 19 errors in height and weight as reported in \code{sitar} issue #7 are also now corrected.
#' The \code{growth} dataset in the \code{fda} package uses heights from the same study.

#' @name berkeley
#' @docType data
#' @format A data frame with 4884 observations on the following 10 variables:
#' \describe{
#' \item{id}{factor with levels 201-278 male and 301-385 female}
#' \item{age}{years, numeric vector}
#' \item{height}{cm, numeric vector}
#' \item{weight}{kg, numeric vector}
#' \item{stem.length}{cm, numeric vector}
#' \item{bi.acromial}{cm, numeric vector}
#' \item{bi.iliac}{cm, numeric vector}
#' \item{leg.circ}{cm, numeric vector}
#' \item{strength}{lb, numeric vector}
#' \item{sex}{factor with level 1 male and level 2 female}
#' }
#' @references
#' Tuddenham RD, Snyder MM. Physical growth of California boys and girls from birth to eighteen years.
#' University of California Publications in Child Development 1954;1:183-364.
#'
#' @keywords datasets
#' @examples
#' data(berkeley)
#' ## frequencies of age of measurement for each variable
#' ## weight and length/height from birth, other variables from 6-8 years
#' ## few measurements after 18 years
#' . <- as.factor(berkeley$age)
#' plot(levels(.), summary(.), type='s', las=1,
#'   xlab='age of measurement (years)', ylab='frequency of measurements')
#' points(levels(.), levels(.) < 0, pch=15)
#' for (i in 3:9) {
#'   .. <- .[!is.na(berkeley[, names(berkeley)[i]])]
#'   lines(levels(..), summary(..), type='s', col=i)
#' }
#' legend('topright', names(berkeley)[c(3:9)], text.col=c(3:9), bty='n', inset=0.04)
"berkeley"
