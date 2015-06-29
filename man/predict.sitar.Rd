\name{predict.sitar}
\alias{predict.sitar}
\title{Predict SITAR model}
\description{
Predict method for \code{sitar} objects, based on \code{predict.lme}.}
\usage{\method{predict}{sitar}(object, newdata, level=1, ..., deriv=0, abc=NULL, xfun=I, yfun=I)}
\arguments{
  \item{object}{
an object inheriting from class \code{sitar}.}
  \item{newdata}{
an optional data frame to be used for obtaining the predictions. It requires named columns for \code{x}, and for \code{id} if \code{level = 1}, matching the names in \code{object}. Any covariates in \code{a.formula}, \code{b.formula} or \code{c.formula} can also be included. By default their values are set to the mean, so when \code{level = 0} the prediction represents the mean curve. Note that factors are coded as instrumental variables for each level, corresponding to the fixed effect coefficients in the model, so their names need the level appending.}
  \item{level}{
an optional integer giving the level of grouping to be used in obtaining the predictions, level 0 corresponding to the population predictions. Defaults to level 1.}
  \item{\dots}{
other optional arguments, including \code{na.action} and \code{naPattern}.}
  \item{deriv}{an optional integer specifying predictions corresponding to either the fitted curve or its first derivative (i.e. velocity), where \code{deriv = 0} (default) specifies the distance curve and \code{deriv = 1} the velocity curve.}
  \item{abc}{an optional named vector or single-row data frame containing values of a subset of \code{a}, \code{b} and \code{c}, default \code{NULL}. If it is set, \code{level} is ignored. It gives predictions for a subject with the specified values of \code{a}, \code{b} and \code{c}, where missing values are set to 0.}
  \item{xfun}{an optional function to apply to \code{x} to convert it back to the original scale, e.g. if x = log(age) then xfun = exp.}
  \item{yfun}{an optional function to apply to \code{y} to convert it back to the original scale, e.g. if y = sqrt(height) then yfun = function(z) z^2.}
}
\details{When \code{deriv = 1} the returned velocity is in units of \code{yfun(y)} per \code{xfun(x)}. So if \code{x} and/or \code{y} are transformed, velocity in the raw units can be obtained by specifying \code{xfun} and/or \code{yfun} to back-transform them appropriately.}
\value{A vector of the predictions, or a list of vectors if \code{asList = TRUE} or \code{length(level) > 1}.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
data(heights)
## fit model
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

## predictions at level 0
predict(m1, newdata=data.frame(age=9:16), level=0)

## predictions at level 1 for subject 5
predict(m1, newdata=data.frame(age=9:16, id=5), level=1)

## velocity predictions for subjects with early and late puberty
vel1 <- predict(m1, deriv=1, abc=data.frame(b=-1))
mplot(age, vel1, id, heights, col=id)
vel1 <- predict(m1, deriv=1, abc=data.frame(b=1))
mplot(age, vel1, id, heights, col=id, add=TRUE)
}
