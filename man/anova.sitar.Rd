\name{anova.sitar}
\alias{anova.sitar}
\title{Compare Likelihoods of Fitted SITAR Objects}
\description{
anova method for \code{sitar} objects, based on \code{anova.lme}.}
\usage{\method{anova}{sitar}(object, ..., test, type, adjustSigma, Terms, L, verbose)}
\arguments{
  \item{object}{an object inheriting from class \code{sitar}.}
  \item{\dots}{other optional fitted model objects.}
  \item{test}{an optional logical value controlling whether likelihood ratio tests should be used.}
  \item{type}{an optional character string specifying the type of sum of squares to be used.}
  \item{adjustSigma}{see \code{\link{anova.lme}}.}
  \item{Terms}{see \code{\link{anova.lme}}.}
  \item{L}{see \code{\link{anova.lme}}.}
  \item{verbose}{an optional logical value.}
}
\value{a data frame inheriting from class "anova.lme".}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
