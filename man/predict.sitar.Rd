\name{predict.sitar}
\alias{predict.sitar}
\title{Predict SITAR model}
\description{
Predict method for \code{sitar} objects, based on \code{predict.lme}.}
\usage{\method{predict}{sitar}(object, newdata, level, ...)}
\arguments{
  \item{object}{
an object inheriting from class \code{sitar}.}
  \item{newdata}{
an optional data frame to be used for obtaining the predictions. It requires a column named x for \code{x}, and a column named id for \code{id} if \code{level = 1}, while \code{a}, \code{b}, \code{c} and any other covariates can also be included. Coefficients for the regression spline are added internally.}
  \item{\dots}{
other optional arguments, including \code{asList}, \code{na.action} and \code{naPattern}.}
}
\value{A vector of the predictions.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
data(heights)
## fit model
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

## predictions at level 0
predict(m1, newdata=data.frame(age=9:16), level=0)

## predictions at level 1 for subject 5
predict(m1, newdata=data.frame(age=9:16, id=5), level=1)
}
