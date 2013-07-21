\name{recalib}
\alias{recalib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{recalib}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
recalib(xc, yc, id = NULL, data, xcnew = NULL, ycnew = NULL, model, from, to)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xc}{
%%     ~~Describe \code{xc} here~~
}
  \item{yc}{
%%     ~~Describe \code{yc} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{xcnew}{
%%     ~~Describe \code{xcnew} here~~
}
  \item{ycnew}{
%%     ~~Describe \code{ycnew} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{from}{
%%     ~~Describe \code{from} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
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
function (xc, yc, id = NULL, data, xcnew = NULL, ycnew = NULL, 
    model, from, to) 
{
    xcall <- model$call.sitar$x
    ycall <- model$call.sitar$y
    mux <- model$mux
    coef <- random.effects(model)
    dcoef <- coef[from, ] - coef[to, ]
    include <- with(data, {
        if (is.null(id)) 
            rep(TRUE, dim(data)[[1]])
        else as.character(id) \%in\% from
    })
    nxc <- length(xc)
    if (is.null(xcnew)) 
        xcnew <- paste("xnew", 1:nxc, sep = "")
    for (i in 1:nxc) {
        data[, xcnew[i]] <- data[, xc[i]]
        x <- c(data[include, xcnew[i]])
        x <- funcall(x, xcall)
        x <- (x - mux - c(coef[from, "b"])) * c(exp(dcoef$c)) + 
            c(coef[to, "b"]) + mux
        data[include, xcnew[i]] <- funcall(x, xcall, TRUE)
    }
    nyc <- length(yc)
    if (is.null(ycnew)) 
        ycnew <- paste("ynew", 1:nyc, sep = "")
    for (i in 1:nyc) {
        data[, ycnew[i]] <- data[, yc[i]]
        y <- c(data[include, ycnew[i]])
        y <- funcall(y, ycall)
        y <- y - dcoef$a
        data[include, ycnew[i]] <- funcall(y, ycall, TRUE)
    }
    invisible(data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
