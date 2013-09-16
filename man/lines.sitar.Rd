\name{lines.sitar}
\alias{lines.sitar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{lines.sitar}
\description{
<<<<<<< HEAD
The lines function for sitar to draw growth curves in various ways.
=======
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
>>>>>>> 62df7c832c5e30769135f5c9dd6ffbba872bb60d
}
\usage{
lines.sitar(model, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
<<<<<<< HEAD
  \item{model}{sitar model object}
  \item{\dots}{
Further graphical parameters (see par) may also be supplied as arguments, particularly, line type, lty, line width, lwd, and color, col.
=======
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
>>>>>>> 62df7c832c5e30769135f5c9dd6ffbba872bb60d
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
<<<<<<< HEAD
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
=======
\author{
%%  ~~who you are~~
}
>>>>>>> 62df7c832c5e30769135f5c9dd6ffbba872bb60d
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
function (model, ...) 
{
    call <- match.call()[-1]
    if (!"newplot" \%in\% names(call)) 
        call <- c(as.list(call), list(newplot = FALSE))
    do.call("plot.sitar", as.list(call))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
