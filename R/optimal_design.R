#' Optimal design for growth reference centile studies
#'
#' Two functions for estimating optimal sample size and sample composition
#' when constructing growth reference centiles.
#'
#' Studies to construct growth reference centiles using \code{GAMLSS} need to
#' be of optimal size. Cole (SMMR, 2020) has shown that
#' the sample composition, i.e. the age distribution of the measurements,
#' needs to be optimised as well as the sample size. Sample composition is defined in terms of the age power
#' lambda which determines the degree of infant oversampling.
#'
#' There are two criteria that determine the optimal sample size and sample
#' composition: the centile of interest (as z-score z) and the required level
#' of precision for that centile (as the z-score standard error SEz).
#'
#' @aliases optimal_design n_agegp
#' @param z z-score on which to base the design, with default -2 which
#' equates to the 2nd centile. If NA, optimal z is calculated from lambda.
#' @param lambda power of age that defines the sample composition.
#' The default NA means calculate optimal lambda from z.
#' @param N total sample size per sex. The default NA means calculate from
#' z or lambda, and SEz if provided.
#' @param SEz target z-score standard error. The default NA means calculate
#' from z or lambda, and N if provided.
#' @param age age at which to calculate SEz. The default 10 returns mean SEz,
#' and if z or lambda are optimal SEz is independent of age.
#' @param minage youngest age (default 0).
#' @param maxage oldest age (default 20).
#' @param n_groups number of age groups (default 20).
#'
#' @return For optimal_design, a tibble with columns:
#' \item{z}{as above.}
#' \item{lambda}{as above.}
#' \item{N}{as above.}
#' \item{SEz}{as above.}
#' \item{age}{as above.}
#' \item{p}{the centile corresponding to z.}
#' \item{plo}{lower 95\% confidence interval for p.}
#' \item{phi}{upper 95\% confidence interval for p.}
#'
#' For n_agegp, a tibble giving the numbers of measurements to be collected
#' per equal width age group, with columns:
#' \item{n_varying}{numbers for equal width age groups.}
#' \item{age}{mean ages for equal width age groups.}
#' \item{n}{number for each unequal width age group (only for longitudinal studies).}
#' \item{age_varying}{target ages for unequal width age groups (only for longitudinal studies).}
#'
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{gamlss} to fit the centiles
#' with the \code{BCCG}, \code{BCT} or \code{BCPE} family.
#' @examples
#' ## estimate optimal sample composition lambda and precision SEz for 9 centiles
#' ## spaced 2/3 of a z-score apart, based on a sample of 10,000 children
#'
#' optimal_design(z = -4:4*2/3, N = 10000)
#'
#' ## calculate age group sizes optimised for centiles from the 50th to the 99.6th
#' ## (or equivalently from the 50th to the 0.4th)
#' ## with a sample of 10,000 children from 0 to 20 years in one-year groups
#'
#' purrr::map_dfc(0:4*2/3, ~{
#'   n_agegp(z = .x, N = 10000) %>%
#'       dplyr::select(!!z2cent(.x) := n_varying)
#'       }) %>%
#'         dplyr::bind_cols(tibble::tibble(age = paste(0:19, 1:20, sep='-')), .)
#' @importFrom dplyr select bind_cols
#' @importFrom purrr map_dfc
#' @importFrom rlang .data !! :=
#' @importFrom stats pnorm
#' @importFrom tibble tibble
#' @export optimal_design
optimal_design <- function(z = -2, lambda = NA, N = NA, SEz = NA,
                           age = 10) {
  stopifnot(!(any(is.na(z) & is.na(lambda))))
  N0 <- 6878
  coef <- c(`(Intercept)` = -3.0878266723299, age = 0.0231399511712584,
            fz2 = 0.536999348905583, lambda = -0.283052090744945,
            `age:fz2` = 0.0204776896422858, `age:lambda` = -0.0603665153280415)
  fz2 <- log(1 + z^2/2)
  if (all(is.na(lambda))) {
    lambda <- -(coef['age'] + coef['age:fz2'] * fz2) / coef['age:lambda']
  }
  if (all(is.na(z))) {
    fz2 <- -(coef['age'] + coef['age:lambda'] * lambda) / coef['age:fz2']
    z <- sqrt(2 * (ifelse(fz2 >= 0, exp(fz2) - 1, NA))) # inadmissible if lambda too small
  }
  SEz0 <- exp(coef['(Intercept)'] + coef['fz2'] * fz2 + coef['lambda'] * lambda +
                (age - 10) * (coef['age'] + coef['age:fz2'] * fz2 + coef['age:lambda'] * lambda))
  if (all(is.na(N))) {
    N <- N0
    if (all(is.na(SEz)))
      SEz <- SEz0
    else
      N <- N * (SEz0 / SEz)^1.85
  } else
    SEz <- SEz0 * (N0 / N)^(1/1.85)
  plo <- pnorm(z - 2 * SEz) * 100
  phi <- pnorm(z + 2 * SEz) * 100
  tibble(z = z, lambda = lambda, N = N, SEz = SEz, age = age, p = z2cent(z), plo = plo, phi = phi)
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "!!", ":="))

#' @importFrom dplyr select
#' @importFrom purrr map_dfc
#' @importFrom tibble tibble
#' @rdname optimal_design
#' @export
n_agegp <- function(z = -2, lambda = NA, N = NA, SEz = NA, minage = 0, maxage = 20, n_groups = 20) {
  results <- optimal_design(z, lambda, N, SEz)
  if (is.na(lambda))
    lambda <- results$lambda
  if (is.na(N))
    N <- results$N
  N <- round(N / n_groups) * n_groups
  table <- tibble(n = N / n_groups,
                  even = 1:(2 * n_groups + 1) %% 2 == 0, # TRUE design age, FALSE age group boundary
                  age = seq(minage, # equal width groups
                            maxage,
                            length = length(.data$even)),
                  age0 = seq(minage^lambda, # unequal width groups
                             maxage^lambda,
                             length = length(.data$even))^(1/lambda)) %>%
    mutate(n0 = .data$age^lambda - dplyr::lag(.data$age^lambda, 2),
           n0 = ifelse(.data$even, 0, .data$n0),
           n0 = N * .data$n0 / sum(.data$n0, na.rm=TRUE),
           n0 = round(dplyr::lead(.data$n0))) %>%
    filter(.data$even) %>% # age group means
    select(n_varying = .data$n0,
           age = .data$age,
           n = .data$n,
           age_varying = .data$age0)
  return(table)
}

