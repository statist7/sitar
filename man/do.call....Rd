\name{do.call...}
\alias{do.call...}
\title{Extended version of do.call to handle repeated arguments}
\description{
do.call fails if arguments are duplicated in the call. do.call... includes ... as an argument and checks for duplicate arguments. Those in ... overwrite earlier arguments.}
\usage{
do.call...(what, args, ..., quote = FALSE, envir = parent.frame())
}
\arguments{
  \item{what}{
either a function or a non-empty character string naming the function to be called.}
  \item{args}{
a list of arguments to the function call. The names attribute of args gives the argument names.}
  \item{\dots}{
a further list of arguments, some of which may also appear in args.}
  \item{quote}{
a logical value indicating whether to quote the arguments.}
  \item{envir}{
an environment within which to evaluate the call. This will be most useful if what is a character string and the arguments are symbols or quoted expressions.}
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
function (what, args, ..., quote = FALSE, envir = parent.frame()) 
{
    ddd <- list(...)
    if (!is.null(ddd)) {
        for (i in names(ddd)) {
            if (i \%in\% names(args)) 
                args[i] <- ddd[i]
            else {
                args <- c(args, i = ddd[i])
                names(args)[length(args)] <- i
            }
        }
    }
    do.call(what, args, quote = quote, envir = envir)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
