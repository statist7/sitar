\name{LMSfit}
\alias{LMSfit}
\title{Estimate LMS curves from tabulated growth reference centiles}
\description{
A function to summarise an existing set of growth reference centiles as the L, M and S curves of the LMS method,
allowing measurements to be converted to z-scores.}
\usage{
LMSfit(x, y, sex = 1, data = parent.frame(), centiles = c(3, 10, 25, 50, 75, 90, 97),
  LMSdf = c(6, 10, 8), L1 = FALSE, ...)}
\arguments{
  \item{x}{
vector of tabulated ages.}
  \item{y}{
matrix of corresponding measurement centiles, e.g. of height or weight, with \code{nrows = length(x)} and \code{ncols = length(centiles)}.}
  \item{sex}{
two-level factor where level 1 corresponds to male and level 2 to female, default 1.}
  \item{data}{
optional data frame containing \code{x}, \code{y} and optionally \code{sex}.}
  \item{centiles}{
vector of centiles corresponding to the columns of \code{y}, default c(3, 10, 25, 50, 75, 90, 97).}
  \item{LMSdf}{
length-3 vector with the cubic smoothing spline estimated degrees of freedom (edf) for the L, M and S curves, default (c(6, 10, 8)).}
  \item{L1}{
optional logical specifying whether the L curve should be constrained to 1, i.e. a Normal distribution, default FALSE.}
  \item{\dots}{
optional graphical parameters for the plots.}
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
