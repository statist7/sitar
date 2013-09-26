\name{velout}
\alias{velout}
\title{velout}
\description{
Quickly identifies putative outliers in a large number of growth curves.}
\usage{
velout(x, y, id, data, lag = 1, velpower = 0.5, limit = 5, linearise = FALSE)
}
\arguments{
  \item{x}{age vector.}
  \item{y}{outcome vector, typically weight or height.}
  \item{id}{factor identifying each subject.}
  \item{data}{data frame containing x, y and id.}
  \item{lag}{lag between measurements for defining growth velocity.}
  \item{velpower}{a value, typically between 0 and 1, defining the power of delta x to use when calculating velocity as delta(y)/delta(x)^velpower. The default of 0.5 is midway between velocity and increment.}
  \item{limit}{the number of standard deviations beyond which a velocity is deemed to be an outlier.}
  \item{linearise}{if TRUE y is converted to a residual about the median curve of y versus x.}
}
\details{
The algorithm works by viewing serial measurements in each growth curve as triplets (A-B-C) and comparing the velocities between them. Velocity is calculated as 

diff(y, lag = lag) / diff(x, lag = lag) ^ velpower

 Missing values for x or y are ignored. If any of the AB, BC or AC velocities are abnormal (more than \code{limit} SDs in absolute value from the median for the dataset) the code for B is non-zero.
}
\value{
Returns a data frame with columns: id, x, y (from the call), code (as described below), vel1, vel2 and vel3 (corresponding to the velocities AB, BC and AC above). The 'data' attribute contains the name of 'data'.

Code is a factor taking values between 0 and 8, with 0 normal (see table below). Values 1-6 depend on the pattern of abnormal velocities, while 7 and 8 indicate a duplicate age (7 for the first in an individual and 8 for later ones). Edge outliers, i.e. first or last for an individual, have just one velocity. Code 4 indicates a conventional outlier, with both AB and BC abnormal and AC normal. Code 6 is an edge outlier. Other codes are not necessarily outliers, e.g. codes 1 or 3 may be adjacent to a code 4. Use \code{codeplot} to look at individual curves, and \code{zapvelout} to delete outliers. 

\tabular{cccl}{
code	AB+BC	AC	interpretation
0	0	0		no outlier
0	0	NA		no outlier
1	0	1		rare pattern
2	1	0		complicated - look at curve
3	1	1		adjacent to simple outlier
4	2	0		single outlier
5	2	1		double outlier
6	1	NA		edge outlier
7	-	-		first duplicate age
8	-	-		later duplicate age
}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{codeplot}}, \code{\link{zapvelout}}}
\examples{
outliers <- velout(age, height, id, heights, limit=3)
}
