\name{predict.sitar}
\alias{predict.sitar}
\title{Predict SITAR model}
\description{
Predict method for \code{sitar} objects, based on \code{predict.lme}.}
\usage{\method{predict}{sitar}(object, newdata, level=1, ...)}
\arguments{
  \item{object}{
an object inheriting from class \code{sitar}.}
  \item{newdata}{
an optional data frame to be used for obtaining the predictions. It requires named columns for \code{x}, and for \code{id} if \code{level = 1}. Any covariates in \code{a.formula}, \code{b.formula} or \code{c.formula} can also be included. By default their values are set to the mean, so when \code{level = 0} the prediction represents the mean curve. Note that factors are coded as instrumental variables, so their names need the level appending.}
  \item{level}{
an optional integer vector giving the level of grouping to be used in obtaining the predictions, level zero corresponding to the population predictions. Defaults to level 1.}
  \item{\dots}{
other optional arguments, including \code{asList}, \code{na.action} and \code{naPattern}.}
}
\value{A vector of the predictions, or a list of vectors if \code{asList = TRUE} or \code{length(level) > 1}.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
data(heights)
## fit model
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

## predictions at level 0
predict(m1, newdata=data.frame(x=9:16), level=0)

## predictions at level 1 for subject 5
predict(m1, newdata=data.frame(x=9:16, id=5), level=1)
}
