\name{lms2z}
\alias{lms2z}
\title{lms2z}
\description{
A function to convert between measurements and z-scores using a growth reference
previously fitted by the LMS method.}
\usage{
lms2z(x, y, sex, data = NULL, measure, ref, toz = TRUE)
}
\arguments{
  \item{x}{
vector of ages}
  \item{y}{
vector of either measurements or z-scores, depending on the value of \code{toz}.}
  \item{sex}{
two-level factor with level 1 male and level 2 female.}
  \item{data}{
dataframe containing \code{x}, \code{y} and \code{sex}.}
  \item{measure}{
character indicating the name of the measurement, depending on the choice of
\code{ref} (see e.g. references \code{uk90}, \code{who06} and \code{ukwhopt}).}
  \item{ref}{
character indicating the name of the growth reference, available as a \code{data} object.}
  \item{toz}{
logical set to TRUE for conversion from measurement to z-score, or FALSE for the reverse.}
}
\details{
Vectors of L, M and S corresponding to \code{x} and \code{sex} are extracted and
passed to either \code{\link{cLMS}} or \code{\link{zLMS}}, depending on \code{toz}.}
\value{A vector of the same length as \code{y} containing the transformed values.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{z2cent}}}
\examples{
data(who06)
## calculate z-score for length 60 cm in boys at age 0:12 months
lms2z(x = 0:12/12, y = 60, sex = 1, data = NULL, measure = 'ht', ref = 'who06')
}
\keyword{arith}
