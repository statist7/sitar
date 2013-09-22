\name{codeplot}
\alias{codeplot}
\alias{zapvelout}
\title{codeplot}
\description{
Handles output from \code{velout} function to display growth curves with outlying points, either plotting or zapping the outliers.}
\usage{
codeplot(outliers, icode = 4, ..., print = TRUE)
zapvelout(outliers, icode)
}
\arguments{
  \item{outliers}{ Data frame returned from velout.}
  \item{icode}{ The code number(s) defining the subset of curves to be displayed or zapped (between 1 and 6).}
  \item{\dots}{ Optional plot parameters.}
  \item{print}{ Option to print as well as plot information on each curve.}
}
\details{
The function \code{velout} identifies putative outliers for \code{y} in \code{data}, \code{codeplot} plots them, and \code{zapvelout} sets missing those confirmed as outliers.  Codes range from 0 (normal) to 8, where 4 and 6 are conventional outliers (see \code{\link{velout}}.
}
\value{
\code{codeplot} returns summary information on each curve with an outlier of the relevant code, and optionally plots the curve. \code{zapvelout} sets to NA values of \code{y} whose code is contained in \code{icode}, and returns the modified data frame.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{
\code{\link{velout}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (outliers, icode = 4, ..., print = TRUE) 
{
    cat("click or ESC in plot window to progress through cases", 
        "\nESC in console then ESC in plot window to quit\n")
    id <- names(outliers)[1]
    x <- names(outliers)[2]
    y <- names(outliers)[3]
    outliers$res <- outliers$code \%in\% icode
    ids <- unique(outliers[outliers$res, 1])
    nids <- length(ids)
    dat <- outliers[outliers[, 1] \%in\% ids, ]
    dat <- within(dat, {
        vel1[is.na(vel1)] <- 999
        vel2[is.na(vel2)] <- 999
        vel3[is.na(vel3)] <- 999
    })
    dat <- na.omit(dat)
    dat <- within(dat, {
        vel1[vel1 == 999] <- NA
        vel2[vel2 == 999] <- NA
        vel3[vel3 == 999] <- NA
    })
    el <- new.env()
    assign("kount", 0, envir = el)
    by(dat, factor(dat[, 1]), function(z) {
        kount <- get("kount", envir = el) + 1
        assign("kount", kount, envir = el)
        if (print) 
            cat("\ncase", kount, "of", nids, "-", id, as.character(z[1, 
                1]), "with", dim(z)[[1]], "points\n")
        inc <- z[z$res, ]
        if (print) 
            print(inc[, 1:7])
        do.call...("plot", list(z[, 2], z[, 3], type = "b", xlab = x, 
            ylab = y, pch = 46), ...)
        do.call...("points", list(inc[, 2], inc[, 3], pch = 8), 
            ...)
        title(paste(id, z[1, 1], "code", paste(unique(inc$code), 
            collapse = " ")))
        locator(1)
    })
    invisible()
  }
}
