#' Convert to/from measurement from/to z-score with growth reference
#'
#' A function to convert between measurements and z-scores using a growth
#' reference previously fitted by the LMS method.
#'
#' Vectors of L, M and S corresponding to \code{x} and \code{sex} are extracted
#' and passed to either \code{\link{cLMS}} or \code{\link{zLMS}}, depending on
#' \code{toz}.
#'
#' @param x vector of ages.
#' @param y vector of either measurements or z-scores, depending on the value
#' of \code{toz}.
#' @param sex two-level factor with level 1 male and level 2 female, of length
#' one or \code{length(x)}.
#' @param measure character indicating the name of the measurement, depending
#' on the choice of \code{ref} (see e.g. references \code{uk90}, \code{who06}
#' and \code{ukwhopt}).
#' @param ref character indicating the name of the growth reference, available
#' as a \code{data} object.
#' @param toz logical set to TRUE for conversion from measurement to z-score,
#' or FALSE for the reverse.
#' @return A vector or matrix, depending on the lengths of \code{x} and
#' \code{y}, containing the transformed values. If the two lengths are the
#' same, or either is one, then a vector is returned. If they are different and
#' not one, then a matrix with \code{length(x)} rows and \code{length(y)}
#' columns is returned. In this latter case the row names are set to \code{x},
#' and if \code{toz} is false the column names are set using \code{z2cent}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{z2cent}}. The LMS method can be fitted to data using
#' the package \code{gamlss} with the \code{BCCG} family, where nu (originally
#' lambda), mu and sigma correspond to L, M and S respectively.
#' @keywords arith
#' @examples
#'
#' ## convert girls' heights data to UK 90 z-scores
#' data(heights)
#' data(uk90)
#' with(heights, LMS2z(age, height, sex = 2, measure = 'ht', ref = 'uk90'))
#'
#' ## construct table of boys weight centiles by age for WHO standard
#' data(who06)
#' zs <- -4:4*2/3 # z-scores for centiles
#' ages <- 0:12/4 # 3-month ages
#' v <- as.data.frame(lapply(as.list(zs), function(z) {
#'   v <- LMS2z(ages, z, sex = 1, measure = 'wt', ref = 'who06', toz = FALSE)
#' }))
#' names(v) <- z2cent(zs)
#' rownames(v) <- ages
#' round(v, 2)
#'
#' @export LMS2z
	LMS2z <- function(x, y, sex, measure, ref, toz=TRUE) {
#	converts measurement y to/from z-score adjusted for x & sex
#		using LMS reference 'ref' for 'measure'
#	x		age
#	y		measurement (or z-score if toz FALSE)
#	sex		sex variable (male=1, female=2)
#	measure label for measurement, one of:
#		'ht' 'wt' 'bmi' 'head' 'sitht' 'leglen' 'waist' 'bfat'
#	ref		name of reference, one of: 'uk90' 'who06'
#	toz		if TRUE returns measurement converted to z-score using ref
#		  	if FALSE returns z-score converted to measurement using ref
	if (!length(sex) %in% c(1, length(x))) stop('sex wrong length for x')
	lms <- paste(c('L', 'M', 'S'), measure, sep='.')
	v <- matrix(nrow=length(x), ncol=4)
	v[, 4] <- sex
	ref <- get(ref)
	for (i in 1:3) {
		for (ix in 1:2) {
			sexvar <- as.numeric(v[, 4]) == ix
			sexref <- as.numeric(ref$sex) == ix
			if (any(sexvar)) v[sexvar, i] <- spline(ref$years[sexref],
				ref[sexref, lms[i]], method='natural', xout=x[sexvar])$y
		}
	}
	cz <- if (toz) zLMS(y, v[, 1], v[, 2], v[, 3])
		else cLMS(y, v[, 1], v[, 2], v[, 3])
	if (!is.null(dim(cz))) {
	  cz <- as.data.frame(cz)
	  if (!toz) names(cz) <- z2cent(y)
	  rownames(cz) <- x
	}
	cz
}


#' LMS conversion to and from z-scores
#'
#' Routines to handle references constructed with the LMS method. Given a set
#' of LMS values, the functions convert z-scores to measurement centiles and
#' vice versa.
#'
#' L, M and S should all be the same length, recycled if necessary. The
#' formulae converting \code{x} to \code{z} and vice versa are:
#' \deqn{z = \frac{(x/M)^L - 1}{L S}}{z = ((x/M)^L - 1)/L/S}
#'
#' \deqn{x = M (1 + L S z)^{1/L})}{x = M (1 + L S z)^(1/L)} where L is reset
#' to 10^-7 if it is zero. \code{x} and \code{z} are usually the same length as
#' L M and S, but can be different. The LMS method is the same as the BCCG
#' family in the \code{gamlss} package, except that lambda in LMS is referred
#' to as nu in BCCG.
#'
#' @aliases cLMS zLMS
#' @param x vector of measurements to be converted to z-scores.
#' @param z vector of z-scores to be converted to measurements.
#' @param L vector of Box-Cox transformation (lambda) values, L in the LMS
#' method.
#' @param M vector of medians (mu), M in the LMS method.
#' @param S vector of coefficients of variation (sigma), S in the LMS method.
#' @return \code{zLMS} and \code{cLMS} each return a vector or matrix,
#' respectively of z-scores and measurement centiles, with the number of rows
#' matching the length of \code{x} or \code{z}, and the number of columns
#' matching the length of L, M and S. If the two lengths are the same, or if
#' either length is 1, a vector is returned.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{z2cent}}, \code{\link{LMS2z}}
#' @keywords arith
#' @examples
#'
#' cLMS(z = -2:2, L = 1:-1, M = 5:7, S = rep(0.1, 3))
#' cLMS(z = -2:2, L = 1:-1, M = 7, S = 0.1)
#' zLMS(x = 6.5, L = 1:-1, M = 5:7, S = rep(0.1, 3))
#'
#' @export cLMS
	cLMS <- function(z, L = 1, M, S) {
	  L0 <- L + 1e-7 * (L == 0)
	  LMS <- data.frame(cbind(L0, M, S))
	  if (length(z) == nrow(LMS) || min(length(z), nrow(LMS)) == 1)
	    drop(with(LMS, M * (1 + L0 * S * z) ^ (1/L0)))
	  else
	    drop(with(LMS, M * (1 + L0 * S %*% t(z)) ^ (1/L0)))
	}

#' @rdname cLMS
#' @export
	zLMS <- function(x, L = 1, M, S) {
	  L0 <- L + 1e-7 * (L == 0)
	  LMS <- data.frame(cbind(L0, M, S))
	  if (length(x) == nrow(LMS) || min(length(x), nrow(LMS)) == 1)
	    drop(with(LMS, ((x / M) ^ L0 - 1) / L0 / S))
	  else
	    drop(with(LMS, (t(x %*% t(1 / M)) ^ L0 - 1) / L0 / S))
	}



#' Express z-scores as centile character strings for plotting
#'
#' Converts z-scores, typically defining centiles in a growth chart, to
#' character strings that can be used to label the centile curves.
#'
#'
#' @param z a scalar or vector of z-scores.
#' @return A character string is returned, the same length as z. Z-scores in
#' the range +/- 3.3 are converted to centiles with one or two significant
#' figures (lower tail) or to their complement (upper tail). For z-scores
#' exceeding 3.3 in absolute value the character consists of "SDS" appended to
#' the z-score rounded to one decimal place.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{cLMS}}
#' @keywords character
#' @examples
#'
#' z2cent(-4:4)
#'
#' @export z2cent
	z2cent <- function(z) {
#	z is z-score
#	returns corresponding centile as label
	np <- ifelse(abs(z) < 2.33, 0, 1)
	ct <- round(pnorm(z) * 100, np)
	mod10 <- ifelse(np == 1, 0, floor(ct %% 10))
	th <- ifelse(mod10 == 0 | mod10 > 4, 4, mod10)
	th <- paste(ct, c('st','nd','rd','th')[th], sep='')
	th[th == '0th'] <- paste('SDS', round(z[th == '0th'], 1), sep='')
	th[th == '100th'] <- paste('SDS', round(z[th == '100th'], 1), sep='+')
	th
}
