\name{lms2z}
\alias{lms2z}
\title{Convert to/from measurement from/to z-score with growth reference}
\description{
A function to convert between measurements and z-scores using a growth reference
previously fitted by the LMS method.}
\usage{
lms2z(x, y, sex, data = parent.frame(), measure, ref, toz = TRUE)
}
\arguments{
  \item{x}{
vector of ages}
  \item{y}{
vector of either measurements or z-scores, depending on the value of \code{toz}.}
  \item{sex}{
two-level factor with level 1 male and level 2 female.}
  \item{data}{
optional data frame containing \code{x}, \code{y} and \code{sex}.}
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
\seealso{\code{\link{z2cent}}. The LMS method can be fitted to data using the package
  \code{gamlss} with the \code{BCCG} family, where nu (originally lambda), mu and sigma
  correspond to L, M and S respectively.}
\examples{
## convert girls' heights data to UK 90 z-scores
data(heights)
data(uk90)
lms2z(age, height, sex = 2, heights, measure = 'ht', ref = 'uk90')

## construct table of boys weight centiles for WHO standard on 3-month age grid
data(who06)
zs <- -4:4*2/3 # z-scores for centiles
ages <- 0:12/4 # 3-month ages
v <- as.data.frame(lapply(as.list(zs), function(z) {
  v <- lms2z(ages, z, sex = 1, measure = 'wt', ref = 'who06', toz = FALSE)
}))
names(v) <- z2cent(zs)
rownames(v) <- ages
round(v, 2)
}
\keyword{arith}
