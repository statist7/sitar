\name{makess}
\alias{makess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{makess}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makess(x, y, xfun, yfun)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{xfun}{
%%     ~~Describe \code{xfun} here~~
}
  \item{yfun}{
%%     ~~Describe \code{yfun} here~~
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
function (x, y, xfun, yfun) 
{
    if (!missing(xfun)) 
        x <- xfun(x)
    if (!missing(yfun)) 
        y <- yfun(y)
    o <- order(x)
    af <- unique(cbind(x[o], y[o]))
    ss <- smooth.spline(af)
    ss1 <- predict(ss, deriv = 1)
    ss2 <- predict(ss, deriv = 2)
    apv <- vector(length = 2)
    ts <- as.ts(ss2$y)
    ts <- ts * lag(ts) < 0
    if (any(ts, na.rm = TRUE)) {
        ts <- c(NA, ts)
        iapv <- which.max(ss1$y * ts)
        o <- (iapv - 1):(iapv + 1)
        mc <- lm(ss1$y[o] ~ poly(ss$x[o], 2, raw = TRUE))$coef
        apv[1] <- -mc[2]/2/mc[3]
        apv[2] <- mc[1] - mc[2]^2/4/mc[3]
        names(apv) <- c("apv", "pv")
    }
    else {
        iapv <- which.max(ss1$y)
        apv[1] <- ss$x[iapv]
        apv[2] <- ss1$y[iapv]
        names(apv) <- c("amv", "mv")
    }
    ss$apv <- apv
    ss
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
