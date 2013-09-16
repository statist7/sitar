\name{funcall}
\alias{funcall}
<<<<<<< HEAD
\title{funcall}
\description{
Applies an expression to vector v, optionally inverting the expression first. For example if the expression is log, funcall returns log(v) if inverse is FALSE, and exp(v) if inverse is TRUE.}
=======
%- Also NEED an '\alias' for EACH other topic documented here.
\title{funcall}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
>>>>>>> 62df7c832c5e30769135f5c9dd6ffbba872bb60d
\usage{
funcall(v, vcall, inverse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
<<<<<<< HEAD
  \item{v}{vector}
  \item{vcall}{expression}
  \item{inverse}{logical}
}
\details{
Inverse covers functions log, exp, sqrt, ^, *, /, +, -.}
\value{
Returns a vector of length v
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
=======
  \item{v}{
%%     ~~Describe \code{v} here~~
}
  \item{vcall}{
%%     ~~Describe \code{vcall} here~~
}
  \item{inverse}{
%%     ~~Describe \code{inverse} here~~
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
>>>>>>> 62df7c832c5e30769135f5c9dd6ffbba872bb60d

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (v, vcall, inverse = FALSE) 
{
    if (length(vcall) == 1) 
        vcall <- substitute(v)
    else {
        if (inverse) {
            fun <- vcall[[1]]
            if (fun == "log") 
                vcall[[1]] <- as.symbol("exp")
            else if (fun == "exp") 
                vcall[[1]] <- as.symbol("log")
            else if (fun == "sqrt") {
                vcall[[1]] <- as.symbol("^")
                vcall[[3]] <- 2
            }
            else if (fun == "^") 
                vcall[[3]] <- 1/vcall[[3]]
            else if (fun == "*") 
                vcall[[1]] <- as.symbol("/")
            else if (fun == "/") 
                vcall[[1]] <- as.symbol("*")
            else if (fun == "+") 
                vcall[[1]] <- as.symbol("-")
            else if (fun == "-") 
                vcall[[1]] <- as.symbol("+")
            else stop("funcall: unknown function")
        }
        vcall[[2]] <- substitute(v)
    }
    eval(vcall)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
