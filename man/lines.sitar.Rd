\name{lines.sitar}
\alias{lines.sitar}
\title{lines.sitar}
\description{
The \code{lines} function for \code{sitar} to draw growth curves in various ways.
}
\usage{lines.sitar(model, ...)}
\arguments{
  \item{model}{sitar model object}
  \item{\dots}{
Further graphical parameters (see par) may also be supplied as arguments, particularly, line type, lty, line width, lwd, and color, col.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{plot.sitar}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, ...) 
{
    call <- match.call()[-1]
    if (!"add" \%in\% names(call)) 
        call <- c(as.list(call), list(add = TRUE))
    do.call("plot.sitar", as.list(call))
  }
}
