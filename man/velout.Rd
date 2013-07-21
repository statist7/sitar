\name{velout}
\alias{velout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{velout}
\description{
Quickly identify putative outliers in a large number of growth curves.}
\usage{
velout(x, y, id, data, lag = 1, velpower = 0.5, limit = 5, linearise = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An age variable.}
  \item{y}{A outcome variable, typically weight or height.}
  \item{id}{A factor identifying each subject.}
  \item{data}{A data frame containing x, y and id.}
  \item{lag}{The lag between measurements for defining growth velocity.}
  \item{velpower}{A value, typically between 0 and 1, defining the power of delta x to use when calculating velocity as delta(y)/delta(x)^velpower. The default of 0.5 is midway between velocity and increment.}
  \item{limit}{The number of standard deviations beyond which a velocity is deemed to be an outlier.}
  \item{linearise}{If TRUE y is converted to a residual about the median curve of y versus x.}
}
\details{
The algorithm works by viewing serial measurements in each growth curve as triplets (A-B-C) and comparing the velocities between them. Velocity is calculated as 
diff(y, lag = lag) / diff(x, lag = lag) ^ velpower
Missing values for x or y are ignored. If any of the AB, BC or AC velocities are abnormal (more than \code{limit} SDs in absolute value from the median for the dataset) the code for B is non-zero.
}
\value{
Returns a data frame with columns: id, x, y (from the call), code (as described below), vel1, vel2 and vel3 (corresponding to the velocities AB, BC and AC above). The 'data' attribute contains the name of 'data'.

Code is a factor taking values between 0 and 8, with 0 normal (see table below). Values 1-6 depend on the pattern of abnormal velocities, while 7 and 8 indicate a duplicate age (7 for the first in an individual and 8 for later ones). Edge outliers, i.e. first or last for an individual, have just one velocity. Code 4 indicates a conventional outlier, with both AB and BC abnormal and AC normal. Code 6 is an edge outlier. Other codes are not necessarily outliers, e.g. codes 1 or 3 may be adjacent to a code 4. Use \code{codeplot} to look at individual curves, and \code{zapvelout} to delete outliers. 

code AB+BC AC	interpretation
  0	 0	0		no outlier
  0	 0	NA		no outlier
  1	 0	1		rare pattern
  2	 1	0		complicated - look at curve
  3	 1	1		adjacent to simple outlier
  4	 2	0		single outlier
  5	 2	1		double outlier
  6	 1	NA		edge outlier
  7	 -	-		first duplicate age
  8	 -	-		later duplicate age
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{
\code{\link{codeplot}}, \code{\link{zapvelout}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, id, data, lag = 1, velpower = 0.5, limit = 5, 
    linearise = FALSE) 
{
    if (missing(x) | missing(y) | missing(id) | missing(data)) 
        stop("x, y, id and data must all be specified")
    data.name <- deparse(substitute(data))
    if (data.name \%in\% search()) 
        stop("data already attached - needs detaching")
    nrow <- dim(data)[1]
    xyid <- c(deparse(substitute(x)), deparse(substitute(y)), 
        deparse(substitute(id)))
    dc <- data[, xyid]
    dc <- na.omit(cbind(dc, count = 1:nrow))
    dc <- dc[order(dc[, 3], dc[, 1]), ]
	if (linearise) {
#	fit spline curve to convert y to residual adjusted for x
		require(quantreg)
		spline.lm <- rq(dc[,2] ~ bs(dc[,1], df=5))
		dc[,2] <- residual(spline.lm)
    }
    dt1 <- diff(dc[, 1], lag = lag)
    vel1 <- diff(dc[, 2], lag = lag)/dt1^velpower
    dt2 <- diff(dc[, 1], lag = lag * 2)
    vel2 <- diff(dc[, 2], lag = lag * 2)/dt2^velpower
    dt1 <- dt1 == 0
    idlev <- as.numeric(dc[, 3])
    dt1[diff(idlev, lag = lag) != 0] <- FALSE
    vel1[diff(idlev, lag = lag) != 0] <- NA
    vel2[diff(idlev, lag = lag * 2) != 0] <- NA
    vel1 <- trunc(vel1/mad(vel1, na.rm = TRUE))
    vel2 <- trunc(vel2/mad(vel2, na.rm = TRUE))
    vel3 <- c(rep(NA, lag), vel2, rep(NA, lag))
    vel2 <- c(vel1, rep(NA, lag))
    vel1 <- c(rep(NA, lag), vel1)
    code <- (as.numeric(abs(vel1) >= limit) + as.numeric(abs(vel2) >= 
        limit)) * 2 + as.numeric(abs(vel3) >= limit)
    dt2 <- c(dt1, rep(FALSE, lag))
    dt1 <- c(rep(FALSE, lag), dt1)
    code[dt2 | dt1] <- 8
    code[dt2 & !dt1] <- 7
    t <- is.na(vel3)
    code[t] <- (as.numeric(!is.na(vel1[t]) & abs(vel1[t]) >= 
        limit) + as.numeric(!is.na(vel2[t]) & abs(vel2[t]) >= 
        limit)) * 6
    dc <- cbind(dc[, c(3, 1, 2, 4)], code, vel1, vel2, vel3)
    mat <- as.data.frame(matrix(nrow = nrow, ncol = dim(dc)[2], 
        dimnames = list(row.names(data), dimnames(dc)[[2]])))
    attr(mat, "data") <- data.name
    mat[dc$count, ] <- dc
    if (is.factor(dc[, 1])) {
        mat[, 1] <- as.factor(mat[, 1])
        levels(mat[, 1]) <- levels(dc[, 1])
    }
    mat$count <- NULL
    mat$code <- factor(mat$code)
    cat("code frequencies\n")
    print(summary(mat$code))
    invisible(mat)
  }
}
