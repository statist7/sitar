\name{nk2df}
\alias{nk2df}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{nk2df}
\description{
A routine to update old-style \code{sitar} objects containing the \code{nk} argument to change it to \code{df}. This can also be done using \code{update}.}
\usage{
nk2df(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
an object of class \code{sitar}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The object \code{obj} is returned, with the parameter \code{nk} in \code{call.sitar} changed to \code{df}.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
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
function (obj) 
{
    if (class(obj)[1] != "sitar") 
        stop(paste("'", deparse(substitute(obj)), "' not a sitar object", 
            sep = ""))
    nc <- match("nk", names(obj$call.sitar))
    if (is.na(nc)) 
        stop("'nk' arg not found")
    names(obj$call.sitar)[nc] <- "df"
    invisible(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
