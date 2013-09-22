\name{heights}\alias{heights}
\docType{data}\title{Serial heights measured in 12 girls}\description{Heights of 12 girls from the Edinburgh Growth Study measured twice a year between 8 and 16 years of age.}\usage{data(heights)}\format{
  A data frame with 124 observations on the following 3 variables.
  \describe{
    \item{\code{id}}{a factor of subject ids (levels 1:12)}
    \item{\code{age}}{a numeric vector of age (years)}
    \item{\code{height}}{a numeric vector of height (cm)}
  }
}
\examples{
data(heights)
plot(heights)
}
\keyword{datasets}