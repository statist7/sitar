\name{lms2z}
\alias{lms2z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{lms2z}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lms2z(x, y, sex, data = NULL, measure, ref, toz = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{sex}{
%%     ~~Describe \code{sex} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{measure}{
%%     ~~Describe \code{measure} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
}
  \item{toz}{
%%     ~~Describe \code{toz} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, sex, data = NULL, measure, ref, toz = TRUE) 
{
    attach(data)
    on.exit(detach(data))
    lms <- paste(c("L", "M", "S"), measure, sep = ".")
    lmsout <- matrix(nrow = length(x), ncol = 3)
    ref <- get(ref)
    for (i in 1:3) {
        for (ix in 1:2) {
            if (sum(as.numeric(sex) == ix, na.rm = TRUE) > 0) 
                lmsout[as.numeric(sex) == ix, i] <- spline(ref$years[ref$sex == 
                  ix], ref[ref$sex == ix, lms[i]], method = "natural", 
                  xout = x[as.numeric(sex) == ix])$y
        }
    }
    if (toz) 
        zLMS(y, lmsout[, 1], lmsout[, 2], lmsout[, 3])
    else cLMS(y, lmsout[, 1], lmsout[, 2], lmsout[, 3])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
