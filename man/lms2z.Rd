\name{lms2z}
\alias{lms2z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{lms2z}
\description{
A function to convert between measurements and z-scores using a growth reference previously fitted by the LMS method.}
\usage{
lms2z(x, y, sex, data = NULL, measure, ref, toz = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of ages}
  \item{y}{
vector of either measurements or z-scores, depending on the value of \code{toz}.}
  \item{sex}{
two-level factor with level 1 male and level 2 female}
  \item{data}{
dataframe containing \code{x}, \code{y} and \code{sex}}
  \item{measure}{
character indicating the name of the measurement, depending on the choice of \code{ref}.}
  \item{ref}{
character indicating the name of the growth reference, available as a data object.}
  \item{toz}{
logical set to TRUE for conversion from measurement to z-score, and FALSE for the reverse.}
}
\details{
Vectors of L, M and S corresponding to \code{x} and \code{sex} are created and passed to either \code{\link{cLMS}} or \code{\link{zLMS}}.}
\value{A vector of the same length as \code{y} containing the transformed values.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{cLMS}}, \code{\link{zLMS}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, sex, data = NULL, measure, ref, toz = TRUE) 
{
    attach(data)
    on.exit(detach(data))
    lms <- paste(c("L", "M", "S"), measure, sep = ".")
    lmsout <- matrix(nrow = length(x), ncol = 3)
    ref <- get(ref)
    for (i in 1:3) {
        for (ix in 1:2) {
            if (sum(as.numeric(sex) == ix, na.rm = TRUE) > 0) 
                lmsout[as.numeric(sex) == ix, i] <- spline(ref$years[ref$sex == 
                  ix], ref[ref$sex == ix, lms[i]], method = "natural", 
                  xout = x[as.numeric(sex) == ix])$y
        }
    }
    if (toz) 
        zLMS(y, lmsout[, 1], lmsout[, 2], lmsout[, 3])
    else cLMS(y, lmsout[, 1], lmsout[, 2], lmsout[, 3])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
