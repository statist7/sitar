\name{varexp}
\alias{varexp}
\title{variance explained by SITAR model(s)}
\description{
varexp calculates the variance explained by one or more SITAR models, based on the residual variance about the cross-sectional spline curve compared to the residual variance about the SITAR spline curve.
}
\usage{
varexp(..., pattern)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{one or more SITAR models}
  \item{pattern}{a character string defining a search term for names of SITAR models}
}
\details{
pattern is passed to \code{ls(pattern=pattern)}, and any filenames so identified are concatenated with any explicitly given filenames. The variance explained is given by
        pc[i] <- 100 * (1 - (sigma2/sigma1)^2)
where sigma1 is the cross-sectional RSD and sigma2 the SITAR RSD.
}
\value{
A vector of percentages, in decreasing order, each named by the corresponding model.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{devadj}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (...) 
{
    ARG <- list(...)
    pc <- numeric(length(ARG))
    names(pc) <- match.call(expand.dots = FALSE)$...
    for (i in 1:length(ARG)) {
        model <- ARG[[i]]
        if (!"ns" \%in\% names(model)) {
            print(paste(names(pc)[i], "ns not found", sep = "$"))
            next
        }
        pc[i] <- 100 * (1 - (model$sigma/summary(model$ns)$sigma)^2)
    }
    round(pc[rev(order(pc))], 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
