	LMS2z <- function(x, y, sex, data = parent.frame(), measure, ref, toz=TRUE) {
#	converts measurement y to/from z-score adjusted for x & sex
#		using LMS reference 'ref' for 'measure'
#	x		age
#	y		measurement (or z-score if toz FALSE)
#	sex		sex variable (male=1, female=2)
#	data	source of x, y and sex
#	measure label for measurement, one of:
#		'ht' 'wt' 'bmi' 'head' 'sitht' 'leglen' 'waist' 'bfat'
#	ref		name of reference, one of: 'uk90' 'who06'
#	toz		if TRUE returns measurement converted to z-score using ref
#			if FALSE returns z-score converted to measurement using ref
	mcall <- as.list(match.call()[c('x', 'y', 'sex')])
	df <- as.data.frame(lapply(mcall, eval, envir = data, enclos = parent.frame()))
	ref <- get(ref)
	lms <- paste(c('L', 'M', 'S'), measure, sep='.')
	v <- matrix(nrow=nrow(df), ncol=3)
	for (i in 1:3) {
		for (ix in 1:2) {
			sexvar <- as.numeric(df[, 3]) == ix
			sexref <- as.numeric(ref$sex) == ix
			if (sum(sexvar) > 0) v[sexvar, i] <- spline(ref$years[sexref],
				ref[sexref, lms[i]], method='natural', xout=df[sexvar, 1])$y
		}
	}
	if (toz) zLMS(df[, 2], v[, 1], v[, 2], v[, 3])
		else cLMS(df[, 2], v[, 1], v[, 2], v[, 3])
}

	zLMS <- function(x, L = 1, M, S) {
	L0 <- L + 1e-7 * (L == 0)
	( (x / M) ^ L0 - 1) / L0 / S
}

	cLMS <- function(z, L = 1, M, S) {
	L0 <- L + 1e-7 * (L == 0)
	M * (1 + L0 * S * z) ^ (1 / L0)
}

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
