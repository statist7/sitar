\name{sitar}
\alias{sitar}
\alias{update.sitar}
\title{Fit SITAR growth curve model}
\description{SITAR is a method of growth curve analysis, based on \pkg{nlme}, that summarises a set of growth
 curves with a mean growth curve as a regression spline, plus a set of up to three fixed and
 random effects (a, b and c) defining how individual growth curves differ from the mean curve.}
\usage{sitar(x, y, id, data, df, knots, fixed = random, random = "a+b+c", a.formula = ~1,
 b.formula = ~1, c.formula = ~1, bounds = 0.04, start, bstart = "mean", xoffset = "mean",
 returndata = FALSE, verbose = FALSE, correlation = NULL, 
 weights = NULL, subset = NULL, method = "ML", na.action = na.fail, 
 control = nlmeControl(returnObject = TRUE), newform = TRUE)
\method{update}{sitar}(object, ..., evaluate = TRUE)
}
\arguments{
  \item{x}{vector of ages.}
  \item{y}{vector of measurements.}
  \item{id}{factor of subject identifiers.}
  \item{data}{data frame containing variables \code{x}, \code{y} and \code{id}.}
  \item{df}{degrees of freedom for cubic regression spline.}
  \item{knots}{vector of values for knots (default \code{df} quantiles of \code{x} distribution).}
  \item{fixed}{character string defining which of a, b and c are fixed effects (default \code{random}).}
  \item{random}{character string defining which of a, b and c are random effects (default \code{"a+b+c"}).}
  \item{a.formula}{formula for fixed effect a (default \code{~ 1}).}
  \item{b.formula}{formula for fixed effect b (default \code{~ 1}).}
  \item{c.formula}{formula for fixed effect c (default \code{~ 1}).}
  \item{bounds}{span of \code{x} for regression spline, or fractional extension of range (default 0.04).}
  \item{start}{optional numeric vector of initial estimates for the fixed effects, 
  or list of initial estimates for the fixed and random effects (see \code{\link{nlme}}).}
  \item{bstart}{optional starting value for fixed effect b (either "mean" (default), "apv" or value).}
  \item{xoffset}{optional value of offset for \code{x} (either "mean" (default), "apv" or value).}
  \item{returndata}{logical which if TRUE causes the model matrix to be returned, or if FALSE (default)
  the fitted model. Setting returndata TRUE is useful in conjunction with \code{subset} and 
  \code{\link{subsample}} for simulation purposes.}
  \item{verbose}{optional logical value to print information on the evolution of the iterative algorithm (see \code{\link{nlme}}).}
  \item{correlation}{optional \code{corStruct} object describing the within-group correlation structure (see \code{\link{nlme}}).}
  \item{weights}{optional \code{varFunc} object or one-sided formula describing the within-group heteroscedasticity structure (see \code{\link{nlme}}).}
  \item{subset}{optional expression indicating the subset of the rows of data that should be used in the fit (see \code{\link{nlme}}).}
  \item{method}{character string, either "REML" or "ML" (default) (see \code{\link{nlme}}).}
  \item{na.action}{function for when the data contain NAs (see \code{\link{nlme}}).}
  \item{control}{list of control values for the estimation algorithm (see \code{\link{nlme}}).}
  \item{newform}{logical to distinguish between old and new forms of \code{sitar} object - default TRUE indicates new form.}
  \item{object}{object of class \code{sitar}.}
  \item{\dots}{further parameters for \code{update} consisting of any of the above \code{sitar} parameters.}
  \item{evaluate}{logical to control evaluation.  If TRUE (default) the expanded 
  \code{update} call is passed to \code{sitar} for evaluation, while if FALSE the expanded call itself is returned.}
}
\details{
The value of \code{xoffset} is ignored unless \code{newform} is FALSE. \code{bstart} (or
 \code{xoffset} if \code{!newform}) allow the origin of b to be varied, which affects its 
 random effect variance. \code{update} updates the model by taking the \code{object} call, 
 adding any new parameters and replacing changed ones. Where possible the fixed and random
 effects of the model being updated are passed via the \code{start} argument.
}
\value{
An object of class \code{sitar} representing the nonlinear mixed-effects model fit. Generic functions such 
 as \code{print}, \code{plot} and \code{summary} have methods to show the results of the fit. The 
 functions \code{resid}, \code{coef}, \code{fitted}, \code{fixed.effects}, and \code{random.effects} can 
 be used to extract some of its components. 
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
data(heights)
## fit simple model
m1 <- sitar(x=age, y=height, id=id, data=heights, df=5)

##	alternatively try sqrt transform for height and increase df
m2 <- update(m1, x=age, y=sqrt(height), df=6)
}
\keyword{ package }
\keyword{ nonlinear }
\keyword{ regression }
\keyword{ models }
