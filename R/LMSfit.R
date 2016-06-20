#' Estimate LMS curves from tabulated growth reference centiles
#' 
#' A function to summarise an existing set of growth reference centiles as the
#' L, M and S curves of the LMS method.
#' 
#' At each age the optimal Box-Cox power Lopt is estimated to render the
#' centiles closest to Normal, and the corresponding median Mopt and
#' coefficient of variation Sopt are derived. The three sets of values are then
#' smoothed across age to give L, M and S.
#' 
#' @param x vector of tabulated ages.
#' @param y matrix of corresponding measurement centiles, e.g. of height or
#' weight, with \code{nrows = length(x)} and \code{ncols = length(centiles)}.
#' @param sex two-level factor where level 1 corresponds to male and level 2 to
#' female.
#' @param data optional data frame containing \code{x}, \code{y} and
#' \code{sex}.
#' @param centiles vector of centiles corresponding to the columns of \code{y},
#' default c(3, 10, 25, 50, 75, 90, 97).
#' @param df length-3 vector with the cubic smoothing spline equivalent degrees
#' of freedom (edf) for the L, M and S curves, default c(6, 10, 8).
#' @param L1 logical constraining the L curve to 1, i.e. a Normal distribution,
#' default FALSE.
#' @param plot logical to plot the estimated L, M and S curves, default TRUE.
#' @param \dots optional graphical parameters for the plots.
#' @return A list with the results: \describe{ \item{list("LMS")}{data frame of
#' sex, x, L, M, S, Lopt, Mopt, Sopt.} \item{list("ey")}{matrix of predicted
#' values of \code{y}.} \item{list("fit")}{matrix of summary statistics for
#' \code{ey}, giving for each column \code{cmean} the mean centile, \code{zSD}
#' the SD of the z-score difference between \code{ey} and \code{y}, and the
#' minimum and maximum z-scores \code{zmin} and \code{zmax}.} }
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{LMS2z}}, \code{\link{z2cent}}. The LMS method can be
#' fitted to data using the package \code{gamlss} with the \code{BCCG} family,
#' where nu (originally lambda), mu and sigma correspond to L, M and S
#' respectively.
#' @keywords arith
#' @examples
#' 
#' ## first construct table of boys weight centiles by age for WHO standard
#' data(who06)
#' zs <- -4:4*2/3 # z-scores for centiles
#' ages <- 0:12/4 # ages 0-3 years by 3 months
#' v <- LMS2z(ages, zs, sex = 1, measure = 'wt', ref = 'who06', toz = FALSE)
#' round(v, 2)
#' 
#' ## then back-calculate the original LMS curves
#' lms <- LMSfit(x=ages, y=v, sex=1, centiles=pnorm(zs)*100, plot=FALSE)
#' 
#' @export LMSfit
LMSfit <- function(x, y, sex, data = parent.frame(), centiles = c(3,10,25,50,75,90,97), df = c(6,10,8), L1 = FALSE, plot=TRUE, ...) {
#	x is a vector of ages, separately by sex
#	y is a matrix, nrows = length(x), ncols = length(centiles)
#	sex is a factor coded 1 for male, 2 for female, same length as x
#	fits a cubic smoothing spline to the empirical L, M and S values
#	df contains the cubic smoothing spline edf for L, M and S
#	L1 = True forces L to be 1
#	... allows for par vars
#	returns: LMS as data frame (L, M and S smoothed values)
#			 ey as expected centiles
#			 fit as mean, SD, min and max of z-scores (mean back-transformed to centile)

	z <- qnorm(centiles/100)
	ez <- ey <- y
	SD <- matrix(nrow = 4, ncol = length(z))
	rownames(SD) <- c("cmean", "zSD", "zmin", "zmax")
	colnames(SD) <- colnames(ey) <- z2cent(z)
	Lopt <- Mopt <- Sopt <- L <- M <- S <- rep(1, length(x))
	for (i in 1:length(x)) {
		c1 <- cor(z, as.numeric(y[i,]))
		c2 <- cor(z, as.numeric(log(y[i,])))
		c3 <- cor(z, as.numeric(-1/y[i,]))
		if (!L1) Lopt[i] <- (c3 - c1) / 2 / (c1 - 2 * c2 + c3)
		res <- lm(I(as.numeric(y[i,])) ^ Lopt[i] ~ z)
		Mopt[i] <- res$coef[1] ^ (1/Lopt[i])
		Sopt[i] <- res$coef[2] / res$coef[1] / Lopt[i]
	}
	mf <- c("Male", "Female")
	LMStitle <- c("L", "M", "S")
	for (sx in 1:2) {
		if (!any(sex == sx)) next
		LMSt <- cbind(Lopt, Mopt, Sopt)[sex == sx,]
		xt <- x[sex == sx]
		for (i in 1:3) {
			sms <- smooth.spline(xt, LMSt[, i], df = df[i])
			if (i == 1 && !L1) L[sex == sx] <- sms$y else
			  if (i == 2) M[sex == sx] <- sms$y else
		    	if (i == 3) S[sex == sx] <- sms$y
      if (plot) {
  			plot(xt, LMSt[, i], xlab = "age", ylab = LMStitle[i], ...)
  			lines(sms, ...)
  			title(main = paste(mf[sx], "sex   df =", df[i]))
			# x1 <- seq(min(xt), max(xt), length = 101)
			# der1 <- predict(sms, x1, der = 1)
			# plot(der1, type = "l", xlab = "age", ylab = "Velocity", ...)
			# plot(xt, sms$yin - sms$y, xlab = "age", ylab = "Residuals",
			# type = "b", ...)
      }
		}
		for (j in 1:length(z)) {
			ey[, j] <- M * (1 + L * S * z[j]) ^ (1/L)
			ez[, j] <- ((y[, j] / M) ^ L - 1) / L / S
			SD[1, j] <- pnorm(mean(ez[, j])) * 100
			SD[2, j] <- sd(ez[, j])
			SD[3, j] <- min(ez[, j])
			SD[4, j] <- max(ez[, j])
		}
	}
	return(list(LMS = data.frame(sex, x, L, M, S, Lopt, Mopt, Sopt),
	ey = ey, fit = SD))
}
