\name{z2cent}
\alias{z2cent}
\title{Express z-scores as centile character strings for plotting}
\description{
Converts z-scores, typically defining centiles in a growth chart, to character strings that can be used to label the centile curves.
}
\usage{z2cent(z)}
\arguments{
  \item{z}{a scalar or vector of z-scores.}
}
\value{
A character string is returned, the same length as z. Z-scores in the range +/- 3.3 are converted to centiles with one or two significant figures (lower tail) or to their complement (upper tail). For z-scores exceeding 3.3 in absolute value the character consists of "SDS" appended to the z-score rounded to one decimal place.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{cLMS}}}
\examples{
z2cent(-4:4)
}
\keyword{character}