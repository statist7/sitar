LMSfit <- function(x, y, sex=1, LMSdf = c(6,10,8), centile = c(3,10,25,50,75,90,97), L1 = FALSE, col = 3, ...) {
#	x is a vector of ages, separately by sex
#	y is a matrix, nrows = length(x), ncols = length(centile)
#	sex is a factor coded 1 for male, 2 for female, same length as x (default 1)
#	fits a cubic smoothing spline to the empirical L, M and S values
#	LMSdf contains the cubic smoothing spline edf for L, M and S
#	L1 = True forces L to be 1
#	col is colour in plot
#	... allows for par vars
#	returns: LMS as data frame (L, M and S smoothed values)
#			 ey as expected centiles
#			 fit as mean, SD, min and max of z-scores (mean back-transformed to centile)

	z <- qnorm(centile/100)
	ez <- ey <- y
	SD <- matrix(nrow = 4, ncol = length(z))
	rownames(SD) <- c("mean", "zSD", "zmin", "zmax")
	colnames(SD) <- colnames(ey) <- centile
	Lopt <- Mopt <- Sopt <- L <- M <- S <- rep(1, length(x))
	for (i in 1:length(x)) {
		c1 <- cor(z, y[i,])
		c2 <- cor(z, log(y[i,]))
		c3 <- cor(z, -1/y[i,])
		if (!L1) Lopt[i] <- (c3 - c1) / 2 / (c1 - 2 * c2 + c3)
		res <- lm(y[i,] ^ Lopt[i] ~ z)
		Mopt[i] <- res$coef[1] ^ (1/Lopt[i])
		Sopt[i] <- res$coef[2] / res$coef[1] / Lopt[i]
	}
	mf <- c("Male", "Female")
	LMStitle <- c("L", "M", "S")
	for (sx in 1:2) {
		if (!any(sex == sx)) next
		LMSt <- cbind(Lopt, Mopt, Sopt)[sex == sx,]
		xt <- x[sex == sx]
		x1 <- seq(min(xt), max(xt), length = 101)
		for (i in 1:3) {
			sms <- smooth.spline(xt, LMSt[,i], df = LMSdf[i])
			if (i == 1 && !L1) L[sex == sx] <- sms$y else
			if (i == 2) M[sex == sx] <- sms$y else 
			if (i == 3) S[sex == sx] <- sms$y
			der1 <- predict(sms, x1, der = 1)
			plot(xt, LMSt[,i], xlab = "age", ylab = LMStitle[i], ...)
			lines(sms, ...)
			title(main = paste(mf[sx], "sex   df =", LMSdf[i]))
			# plot(der1, type = "l", xlab = "age", ylab = "Velocity", col = col)
			# plot(xt, sms$yin - sms$y, xlab = "age", ylab = "Residuals", 
			# type = "b", col = col)
		}
		for (j in 1:length(z)) {
			ey[,j] <- M * (1 + L * S * z[j]) ^ (1/L)
			ez[,j] <- ((y[, j] / M) ^ L - 1) / L / S
			SD[1,j] <- pnorm(mean(ez[, j])) * 100
			SD[2,j] <- sd(ez[, j])
			SD[3,j] <- min(ez[, j])
			SD[4,j] <- max(ez[, j])
		}
	}
	return(list(LMS = data.frame(sex, x, L, M, S, Lopt, Mopt, Sopt), 
	ey = ey, fit = SD))
}