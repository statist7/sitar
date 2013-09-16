\name{funcall}
\alias{funcall}
\title{funcall}
\description{
Applies an expression to vector v, optionally inverting the expression first. For example if the expression is log, funcall returns log(v) if inverse is FALSE, and exp(v) if inverse is TRUE.}
\usage{
funcall(v, vcall, inverse = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
