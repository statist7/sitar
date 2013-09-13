\name{zapvelout}
\alias{zapvelout}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{zapvelout}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zapvelout(outliers, icode)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outliers}{
%%     ~~Describe \code{outliers} here~~
}
  \item{icode}{
%%     ~~Describe \code{icode} here~~
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
function (outliers, icode) 
{
    data <- get(attr(outliers, "data"))
    if (dim(data)[1] != dim(outliers)[1]) 
        stop("numbers of rows in data and outliers differ")
    y <- names(outliers)[3]
    zap <- outliers$code \%in\% icode
    cat(sum(zap), y, "values set missing\n")
    data[outliers$code \%in\% icode, y] <- NA
    invisible(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
