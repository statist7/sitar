\name{summary.sitar}
\alias{summary.sitar}
\title{Create summary of SITAR model}
\description{
A \code{summary} method for \code{sitar} objects based on \code{\link{summary.lme}}.}
\usage{
\method{summary}{sitar}(object, adjustSigma = TRUE, verbose = FALSE, ...)
}
\arguments{
  \item{object}{an object inheriting from class \code{sitar}.}
  \item{adjustSigma}{an optional logical value (see \code{\link{summary.lme}}.}
  \item{verbose}{an optional logical value used to control the amount of output in the \code{print.summary.sitar} method.}
  \item{\dots}{some methods for this generic require additional arguments. None are used in this method.}
}
\value{
an object inheriting from class \code{summary.sitar} with all components included in \code{object} (see \code{\link{lmeObject}} for a full description of the components) plus the components detailed for \code{\link{summary.lme}} and the following components:
\item{x.adj}{vector of length \code{x} in \code{object} with \code{x} values adjusted for subject-specific random effects b and c.}
\item{y.adj}{vector of length \code{y} in \code{object} with \code{y} values adjusted for subject-specific random effects a.}
\item{apv}{length 2 vector giving respectively age at peak velocity and peak velocity based on the fitted distance curve (using transformed \code{x} and \code{y} where specified).}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
