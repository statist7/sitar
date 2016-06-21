#' Serial heights measured in 12 girls
#'
#' Heights of 12 girls from the Chard Growth Study measured twice a year
#' between 8 and 16 years of age.
#'
#'
#' @name heights
#' @docType data
#' @format A data frame with 124 observations on the following 4 variables:
#' \describe{
#' \item{id}{factor of subject ids (levels 1:12).}
#' \item{age}{vector of ages (years).}
#' \item{height}{vector of heights (cm).}
#' \item{men}{vector of ages at menarche (years), where negative values
#' are right censored.} }
#' @keywords datasets
#' @examples
#'   require(graphics)
#'   data(heights)
#'   coplot(height ~ age | id, data = heights, panel=panel.smooth,
#'     show.given=FALSE, xlab='age (years)', ylab='height (cm)', pch=19)
"heights"
