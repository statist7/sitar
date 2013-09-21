\name{nk2df}
\alias{nk2df}
\title{nk2df}
\description{
A routine to update old-style \code{sitar} objects containing the \code{nk} argument to change it to \code{df}. This can also be done using \code{update}.}
\usage{nk2df(obj)}
\arguments{
  \item{obj}{
an object of class \code{sitar}.}
}
\value{
The object \code{obj} is returned, with the parameter \code{nk} in \code{call.sitar} changed to \code{df}.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
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
