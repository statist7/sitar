\name{revise}
\alias{revise}
\title{revise}
\description{Change old-style \code{sitar} obect to new-style. Outdated, and better to use \code{update}.}
\usage{
revise(obj)
}
\arguments{
  \item{obj}{
an old-style \code{sitar} object.}
}
\value{
Returns an new-style \code{sitar} object.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{
\code{\link{nk2df}}, \code{\link{update}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obj) 
{
    objname <- deparse(substitute(obj))
    if (grepl("[", objname, fixed = TRUE)) 
        stop("object name must not be indexed")
    if (length(class(obj)) == 1 && class(obj) == "by") {
        new <- vector("list", length(obj))
        class(new) <- class(obj)
        for (i in 1:length(obj)) {
            if (is.null(obj[[i]]$sitar)) 
                stop("not an old-style sitar object")
            new[[i]] <- obj[[i]]$sitar
            if (!is.null(obj[[i]]$ns)) 
                new[[i]]$ns <- obj[[i]]$ns
            if (!is.null(obj[[i]]$call)) 
                new[[i]]$call.sitar <- obj[[i]]$call
            if (!is.null(obj[[i]]$xoffset)) 
                new[[i]]$xoffset <- obj[[i]]$mux
            if (!"sitar" \%in\% class(new[[i]])[1]) 
                class(new[[i]]) <- c("sitar", class(new[[i]]))
            nc <- match("nk", names(new[[i]]$call.sitar))
            if (is.na(nc)) 
                stop("'nk' arg not found")
            names(new[[i]]$call.sitar)[nc] <- "df"
        }
    }
    else {
        if (is.null(obj$sitar)) 
            stop("not an old-style sitar object")
        new <- obj$sitar
        if (!is.null(obj$ns)) 
            new$ns <- obj$ns
        if (!is.null(obj$call)) 
            new$call.sitar <- obj$call
        if (!is.null(obj$xoffset)) 
            new$xoffset <- obj$mux
        if (!"sitar" \%in\% class(new)[1]) 
            class(new) <- c("sitar", class(new))
        nc <- match("nk", names(new$call.sitar))
        if (is.na(nc)) 
            stop("'nk' arg not found")
        names(new$call.sitar)[nc] <- "df"
    }
    new
  }
}
