\name{recalib}
\alias{recalib}
\title{recalib}
\description{
A function to recalibrate x,y data using SITAR random effects}
\usage{
recalib(xc, yc, id = NULL, data, xcnew = NULL, ycnew = NULL, model, from, to)
}
\arguments{
  \item{xc}{character vector defining column name(s) of \code{x} data to be recalibrated.}
  \item{yc}{character vector defining column name(s) of \code{y} data to be recalibrated.}
  \item{id}{factor defining \code{from} and \code{to} rows. If \code{NULL} then recalibrate all rows.}
  \item{data}{dataframe containing \code{xc}, \code{yc} and \code{id}.}
  \item{xcnew}{column names for replacement columns \code{xc}. If default \code{NULL} then use names xcnew1... .}
  \item{ycnew}{column names for replacement columns \code{yc}. If default \code{NULL} then use names ycnew1... .}
  \item{model}{\code{sitar} model defining the random effects to be used for recalibration.}
  \item{from}{level of \code{id} defining existing data (must be a single row in \code{coef{model}}).}
  \item{to}{level of \code{id} defining data to be recalibrated (must be a single row in \code{coef{model}}).}
}
\details{
\code{recalib} recalibrates the values of \code{xc} and \code{yc} based on \code{model}. 
\code{xc} values are changed to \code{(xc - c(coef[from, 'b'])) * exp(coef[from, 'c'] - coef[to, 'c']) + coef[to, 'b']}. 
\code{yc} values are changed to \code{yc - coef[from, 'a'] + coef[to, 'a']}.
}
\value{
Returns the dataframe \code{data} with the \code{from} rows of \code{xc} and \code{yc} recalibrated.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
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
