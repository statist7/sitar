\name{heights}
\alias{heights}
\docType{data}
\title{Serial heights measured in 12 girls}
\description{
Heights of 12 girls from the Chard Growth Study measured twice a year between 8 and
16 years of age.
}
\usage{data(heights)}
\format{
  A data frame with 124 observations on the following 4 variables:
  \describe{
    \item{\code{id}}{factor of subject ids (levels 1:12).}
    \item{\code{age}}{vector of ages (years).}
    \item{\code{height}}{vector of heights (cm).}
    \item{\code{men}}{vector of ages at menarche (years), where negative values are right censored.}
  }
}
\examples{
  require(graphics)
  data(heights)
  coplot(height ~ age | id, data = heights, panel=panel.smooth,
    show.given=FALSE, xlab='age (years)', ylab='height (cm)', pch=19)
}
\keyword{datasets}
