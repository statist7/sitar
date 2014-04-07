\name{subsample}
\alias{subsample}
\title{Sample from SITAR dataset}
\description{
A function to sample from a SITAR dataset for experimental design purposes. Two different 
sampling schemes are offered, based on the values of \code{id} and \code{x}.}
\usage{
subsample(x, id, data, prob = 1, xlim = NULL)
}
\arguments{
  \item{x}{
vector of age.}
  \item{id}{
factor of subject identifiers.}
  \item{data}{
dataframe containing \code{x} and \code{id}.}
  \item{prob}{
scalar defining sampling probability. See Details.}
  \item{xlim}{
length 2 vector defining range of \code{x} to be selected. See Details.}
}
\details{
With the first sampling scheme \code{xlim} is set to \code{NULL} (default), and rows of
 \code{data} are sampled with probability \code{prob} without replacement.
  With the second sampling scheme \code{xlim} is set to a range within \code{range(x)}.
 Subjects \code{id} are then sampled with probability \code{prob} without replacement,
 and all their rows where \code{x} is within \code{xlim} are selected. 
 The second scheme is useful for testing the power of the model to predict later growth 
 when data only up to a certain age are available. Setting \code{xlim} to \code{range(x)} 
 allows data to be sampled by subject. The returned value can be used as the \code{subset} 
 argument in \code{sitar} or \code{update.sitar}.}
\value{
Returns a logical the length of \code{x} where \code{TRUE} indicates a sampled value.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{sitar}}}
\examples{
## draw 50% random sample
s50 <- subsample(age, id, heights, prob=0.5)

## truncate age range to 7-12 for 50% of subjects
t50 <- subsample(age, id, heights, prob=0.5, xlim=c(7, 12))
}
