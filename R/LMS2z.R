#' Convert to/from measurement from/to z-score with growth reference
#'
#' A function to convert between measurements and z-scores using a growth
#' reference previously fitted by the LMS method.
#'
#' Growth references fitted by the LMS method consist of a table of L, M and S
#' values by age and sex. Vectors of L, M and S corresponding to \code{x} and
#' \code{sex} are extracted using cubic interpolation and passed to either
#' \code{\link{cLMS}} or \code{\link{zLMS}}, depending on \code{toz}.
#'
#' Disjunct references are supported, where there is a disjunction in the
#' centiles at a particular age. This may be because the measurement changes,
#' e.g. from length to height, or because two different references have been
#' joined together. The disjunction is flagged by including two rows at the
#' common age, but with different L, M and S values, and measurements at this
#' age are ascribed to the older reference. For example the \code{who06}
#' reference has a disjunction at 2 years reflecting the switch from length to
#' height. As a result height at just below and just above 2 years returns a
#' different z-score.
#'
#' @param x vector of ages in units of years.
#' @param y vector or one-column matrix of either measurements or z-scores,
#'   depending on the value of \code{toz}.
#' @param sex vector where 1/2 = males/females = boys/girls = TRUE/FALSE, based
#'   on the uppercased first character of the string.
#' @param measure unique measurement name, as character string, the choice depending on
#'   the choice of \code{ref} (see e.g. references \code{uk90}, \code{who06} and
#'   \code{ukwhopt}).
#' @param ref unique growth reference, either as name or character string, available as
#'   a \code{data} object or data frame (e.g. \code{uk90}, \code{who06} or
#'   \code{ukwhopt}).
#' @param toz logical set to TRUE for conversion from measurement to z-score, or
#'   FALSE for the reverse.
#' @param LMStable logical set to TRUE to return the associated LMS table as a
#'   data frame in attribute \code{LMStable}.
#' @return A vector or matrix containing the transformed values. If \code{y} is
#'   a vector then a vector of \code{length(x)} is returned, else if \code{y} is
#'   a one-column matrix then a matrix is returned, with \code{length(x)} rows
#'   and \code{length(y)} columns. The matrix row names are set to \code{x}, and
#'   the column names to either \code{y} or if \code{toz} is FALSE,
#'   \code{z2cent(y)}. If LMStable is TRUE the associated LMS table is returned
#'   as a data frame in attribute \code{LMStable}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{z2cent}}. The LMS method can be fitted to data using the
#'   package \code{gamlss} with the \code{BCCG} or \code{BCCGo} family, where nu
#'   (originally lambda), mu and sigma correspond to L, M and S respectively.
#' @keywords arith
#' @examples
#'
#' ## convert girls' heights data to UK 90 z-scores
#' data(heights)
#' data(uk90)
#' with(heights, LMS2z(age, height, sex = 2, measure = 'ht', ref = 'uk90'))
#'
#' ## construct table of boys' weight centiles by age for WHO standard
#' data(who06)
#' zs <- -4:4*2/3 # z-scores for 9 centiles
#' ages <- 0:20/4 # 3-month ages to 5 years
#' LMS2z(ages, as.matrix(zs), sex = 'm', measure = 'wt', ref = who06,
#'   toz = FALSE, LMStable = TRUE)
#'
#' @export
LMS2z <- function(x, y, sex, measure, ref, toz=TRUE, LMStable=FALSE) {
  # check sex coded correctly
  test_sex <- function(sex) {
    # check sex contains only 1/M/B or 2/F/G
    fsex <- toupper(substr(sex, 1, 1))
    fsex <- factor(fsex, levels = c(1:2, 'M', 'F', 'B', 'G', 'T'))
    levels(fsex) <- c(rep(1:2, 3), 1)
    as.numeric(fsex)
  }
  measure <- unique(measure)
  stopifnot('measure not unique' = length(measure) == 1L)
  xy <- data.frame(x, sex=test_sex(sex))
  x <- xy$x
  sex <- xy$sex
  stopifnot('x or sex missing' = any(!is.na(x)) & any(!is.na(sex)))
  LMS <- c('L', 'M', 'S')
  LMSnames <- paste(LMS, measure, sep='.')
  if (is.character(ref)) {
    ref <- unique(ref)
    if (length(ref) != 1)
      stop('ref not unique')
    ref <- get(ref)
  }
  x[x < min(ref$years) | x > max(ref$years)] <- NA
  for (ix in 1:2) {
    sexvar <- sex == ix
    if (any(sexvar)) {
      # unique omits duplicated rows # na.omit omits blank rows
      refx <- na.omit(unique(ref[ref$sex == ix, c('years', LMSnames)]))
      # check for location of any duplicated ages
      nref <- which(diff(refx$years) == 0)
      nref <- c(nref, nrow(refx))
      end <- 0L
      for (i in seq_along(nref)) {
        start <- end + 1L
        end <- nref[[i]]
        # age range from start to end
        # if age duplicated end value is overwritten by start value on next pass
        refrange <- x[sexvar] >= refx$years[[start]] & x[sexvar] <= refx$years[[end]]
        refrange[is.na(refrange)] <- FALSE
        refrange <- which(sexvar)[refrange]
        if (length(refrange) > 0L) {
        # cubic interpolation of L, M and S
          xy[refrange, LMS] <- lapply(LMSnames, function(lms) {
              with(refx[start:end, ], spline(years, get(lms), method='natural',
                                           xout=x[refrange])$y)
          })
        }
      }
    }
  }
  cz <- do.call(ifelse(toz, 'zLMS', 'cLMS'), c(list(y), xy[, LMS]))
  if (is.matrix(cz)) {
    dimnames(cz) <- if (toz)
      list(x, y)
    else
      list(x, z2cent(y))
  }
  if (LMStable) {
    names(xy)[[1L]] <- 'years'
    attr(cz, 'LMStable') <- xy
  }
  cz
}

#' LMS conversion to and from z-scores
#'
#' Routines to handle references constructed with the LMS method. Given a set
#' of LMS values, the functions convert z-scores to measurement centiles and
#' vice versa.
#'
#' L, M and S -- and if vectors then \code{x} and \code{z} --
#' should all be the same length, recycled if necessary.
#' The formulae converting \code{x} to \code{z} and vice versa are:
#' \deqn{z = \frac{(x/M)^L - 1}{L S}}{z = ((x/M)^L - 1)/L/S}
#'
#' \deqn{x = M (1 + L S z)^{1/L})}{x = M (1 + L S z)^(1/L)} where L is reset
#' to 10^-7 if it is zero. The LMS method is the same as the BCCG
#' family in the \code{gamlss} package, except that lambda in LMS is referred
#' to as nu in BCCG.
#'
#' @aliases cLMS zLMS
#' @param x vector or one-column matrix of measurements to be converted to z-scores.
#' @param z vector or one-column matrix of z-scores to be converted to measurements.
#' @param L vector of Box-Cox transformation (lambda) values, L in the LMS
#' method.
#' @param M vector of medians (mu), M in the LMS method.
#' @param S vector of coefficients of variation (sigma), S in the LMS method.
#' @return If \code{x} and \code{z} are vectors \code{zLMS} and \code{cLMS}
#' each return a vector, respectively of z-scores and measurement centiles, with length
#' matching the length of (the longest of) \code{x} or \code{z}, L, M and S.
#' If \code{x} or \code{z} are matrices \code{zLMS} and \code{cLMS} each return a matrix,
#' the number of rows matching the length of (the longest of) L, M and S,
#' and the number of columns matching the length of \code{x} or \code{z}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{z2cent}}, \code{\link{LMS2z}}, \code{\link{pdLMS}}
#' @keywords arith
#' @examples
#'
#' cLMS(z = as.matrix(-2:2), L = 1:-1, M = 5:7, S = rep(0.1, 3))
#' cLMS(z = 0:2, L = 1:-1, M = 7, S = 0.1)
#' cLMS(z = as.matrix(0:2), L = 1:-1, M = 7, S = 0.1)
#' zLMS(x = 6.5, L = 1:-1, M = 5:7, S = rep(0.1, 3))
#'
#' @export
cLMS <- function(z, L = 1, M, S) {
  L[L == 0] <- 1e-7
  if (is.vector(z))
    with(data.frame(z, L, M, S),
         M * (1 + L * S * z) ^ (1/L))
  else
    with(data.frame(L, M, S),
         M * (1 + L * outer(S, as.vector(z))) ^ (1/L))
}

#' @rdname cLMS
#' @export
zLMS <- function(x, L = 1, M, S) {
  L[L == 0] <- 1e-7
  if (is.vector(x))
    with(data.frame(x, L, M, S),
         ((x/M) ^ L - 1) / L / S)
  else
    with(data.frame(L, M, S),
         (t(outer(as.vector(x), M, `/`)) ^ L - 1) / L / S)
}

#' Express z-scores as centile character strings for plotting
#'
#' Converts z-scores, typically defining centiles in a growth chart, to
#' character strings that can be used to label the centile curves.
#'
#' @param z a scalar or vector of z-scores.
#' @return A character string is returned, the same length as z. Z-scores between
#' the 1st and 99th centile are converted to centiles with one or two significant
#' figures (lower tail) or to their complement (upper tail). For larger z-scores
#' in absolute value the character consists of "SDS" appended to the z-score
#' rounded to one decimal place.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{cLMS}}
#' @keywords character
#' @examples
#'
#' z2cent(-4:4)
#' z2cent(qnorm(0:100/100))
#'
#' @export
z2cent <- function(z) {
  np <- abs(z) > qnorm(0.99)
  ct <- round(pnorm(z) * 100, np)
  mod10 <- ifelse(np, 0, floor(ct %% 10))
  th <- ifelse(mod10 == 0 | mod10 > 4 | (ct > 10 & ct < 14), 4, mod10)
  th <- ifelse(is.na(th), 'NA', paste0(ct, c('st', 'nd', 'rd', 'th')[th]))
  th[th == '0th'] <- paste0('SDS', round(z[th == '0th'], 1))
  th[th == '100th'] <- paste('SDS', round(z[th == '100th'], 1), sep='+')
  th
}
