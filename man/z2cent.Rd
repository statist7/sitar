\name{z2cent}
\alias{z2cent}
\title{Express z-scores as centile character strings for plotting}
\description{
Converts z-scores, typically defining centiles in a growth chart, to character strings that can be used to label the centile curves.
}
\usage{
z2cent(z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{a scalar or vector of z-scores.}
}
\value{
A character string is returned, the same length as z. Z-scores in the range +/- 3.3 are converted to centiles with one or two significant figures (lower tail) or to their complement (upper tail). For z-scores such that \code{qnorm(z)} > 3.3 in absolute value, the character consists of "SDS" appended to the z-score rounded to one decimal place.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{link{cLMS}}, \code{link{zLMS}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z) 
{
    np <- ifelse(abs(z) < 2.33, 0, 1)
    ct <- round(pnorm(z) * 100, np)
    mod10 <- ifelse(np == 1, 0, floor(ct\%\%10))
    th <- ifelse(mod10 == 0 | mod10 > 4, 4, mod10)
    th <- paste(ct, c('st','nd','rd','th')[th], sep='')
    th[th == '0th'] <- paste('SDS', round(z[th == '0th'], 1), sep='')
    th[th == '100th'] <- paste('SDS', round(z[th == '100th'], 1), sep='+')
    th
  }
}
