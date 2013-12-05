	lms2z <- function(x, y, sex, data=NULL, measure, ref, toz=TRUE) {
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
	mcall <- match.call()[-1]
	vars <- c('x', 'y', 'sex')
	if (!is.null(data)) df <- as.data.frame(lapply(as.list(mcall[vars]), eval, envir = data, enclos = parent.frame()))
	else df <- as.data.frame(cbind(x, y, sex))
	lms <- c('L', 'M', 'S')
	lmsv <- paste(lms, measure, sep='.')
	ref <- get(ref)
	v <- df[, 1]
	for (i in 1:3) {
		for (ix in 1:2) {
			sexvar <- as.numeric(df[, 3]) == ix
			sexref <- as.numeric(ref$sex) == ix
			if (sum(sexvar) > 0) v[sexvar] <- spline(ref$years[sexref], 
				ref[sexref, lmsv[i]], method='natural', xout=df[sexvar, 1])$y
		}
		assign(lms[[i]], v)
	}
	if (toz) zLMS(df[, 2], L, M, S)
		else cLMS(df[, 2], L, M, S)
}

	zLMS <- function(x, L, M, S, data=NULL) {
		with(data, {
			L0 <- L + 1e-7 * (L == 0)
			( (x / M) ^ L0 - 1) / L0 / S
		} )
	}
	
	cLMS <- function(z, L, M, S, data=NULL) {
		with(data, {
			L0 <- L + 1e-7 * (L == 0)
			c <- M * (1 + L0 * S * z) ^ (1 / L0)
			# c <- M * (1 + L0 * S %o% z) ^ (1 / L0)
			if (length(z) == 1) as.numeric(c)
				else c
		} )
	}
	
	z2cent <- function(z)
#	z is z-score
#	returns corresponding centile as label
{
	np <- ifelse(abs(z) < 2.33, 0, 1)
	ct <- round(pnorm(z) * 100, np)
	mod10 <- ifelse(np == 1, 0, floor(ct %% 10))
	th <- ifelse(mod10 == 0 | mod10 > 4, 4, mod10)
	th <- paste(ct, c('st','nd','rd','th')[th], sep='')
	th[th == '0th'] <- paste('SDS', round(z[th == '0th'], 1), sep='')
	th[th == '100th'] <- paste('SDS', round(z[th == '100th'], 1), sep='+')
	th
}

#	convert file of age, sex and weight to z-scores
data <- data.frame(age=1:20, sex=rep(1:2, times=10), wt=rnorm(20, mean=40, sd=5))
lms2z(age, wt, sex, data, measure = 'wt', ref = 'uk90')

zs <- -4:4*2/3 # z-scores for height centiles
ages <- 0:12/4 # 3-month ages
v <- as.data.frame(lapply(as.list(zs), function(z) {
	v <- lms2z(ages, z, sex = 1, measure = 'ht', ref = 'who06', toz=FALSE)
}))
names(v) <- z2cent(zs)
rownames(v) <- ages
round(v, 2)