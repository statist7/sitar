\name{codeplot}
\alias{codeplot}
\alias{zapvelout}
\title{Plot and zap velocity outliers in growth curves}
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
The function \code{velout} identifies putative outliers for \code{y} in \code{data}, \code{codeplot} plots them, and \code{zapvelout} sets missing those confirmed as outliers.  Codes range from 0 (normal) to 8, where 4 and 6 are conventional outliers (see \code{\link{velout}}).
}
\value{
\code{codeplot} returns summary information on each curve with an outlier of the relevant code, and optionally plots the curve. \code{zapvelout} sets to NA values of \code{y} whose code is contained in \code{icode}, and returns the modified data frame.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{
\code{\link{velout}}
}
\examples{
## identify outliers
outliers <- velout(age, height, id, heights, limit=2)

## plot outliers with code 4 or 6
codeplot(outliers, icode=c(4,6))

## set the 8 outliers missing
newheights <- zapvelout(outliers, icode=6)
}
