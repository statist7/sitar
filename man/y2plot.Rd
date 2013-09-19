\name{y2plot}
\alias{y2plot}
\title{y2plot}
\description{
Function to plot two y variables, y1 and y2, against a single x variable, with the y1 and y2 axes on the left and right of the plot respectively.}
\usage{
y2plot(x, y1, y2 = NULL, labels, y2par = NULL, add = FALSE, xy, xlegend = "topleft", inset = 0.05, ...)
}
\arguments{
  \item{x}{vector of ages.}
  \item{y1}{vector of measurements for plotting on left y axis.}
  \item{y2}{optional vector of measurements for plotting on right y axis.}
  \item{labels}{character vector containing labels for x, y1 and y2.}
  \item{y2par}{named list of par arguments to format the y2 axis.}
  \item{add}{logical flag to specify if a new plot (with axes etc) is to be drawn (FALSE) or an existing plot is to be added to (TRUE).}
  \item{xy}{obsolete.}
  \item{xlegend}{position for legend.}
  \item{inset}{inset for legend.}
  \item{\dots}{optional par arguments.}
}
\details{
y2plot draws up to two superimposed plots, one with the y axis on the left and the other on the right, with suitable adjustment for \code{par('mar')} and including a legend. The format for y1 is controlled by par arguments, and that for y2 by the list y2par. }
\value{
Returns the list
  \item{usr}{par('usr') for y1 axis}
  \item{usr2}{par('usr') for y2 axis}
In addition the variable par.usr2, equal to usr2, is created in globalenv(), which obviates the need for xy.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{plot.sitar}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y1, y2 = NULL, labels, y2par = NULL, newplot = TRUE, 
    xy, xlegend = "topleft", inset = 0.05, ...) 
{
    if (missing(labels)) 
        labels <- c(deparse(substitute(x)), deparse(substitute(y1)), 
            deparse(substitute(y2)))
    if (is.null(y2par$ylab)) 
        y2par$ylab <- labels[3]
    ypar <- list(...)
    opar <- par(no.readonly = TRUE)
    lty <- 1:2
    lwd <- c(1, 1)
    col <- rep(par("col"), 2)
    for (i in c("lty", "lwd", "col")) {
        j <- get(i)
        if (i \%in\% names(ypar)) {
            j[1] <- ypar[i]
            assign(i, unlist(j))
        }
        else ypar[i] <- j[1]
    }
    if (!missing(y2)) {
        for (i in c("lty", "lwd", "col")) {
            j <- get(i)
            if (i \%in\% names(y2par)) {
                j[2] <- y2par[i]
                assign(i, unlist(j))
            }
            else y2par[i] <- j[2]
        }
    }
    if (newplot) {
        if (!is.null(ypar$mar)) {
            mar <- ypar$mar
        }
        else {
            mar <- c(5, 4, 4, 2) + 0.1
            if (!missing(y2)) 
                mar[4] <- 4.1
        }
        par(mar = mar)
        do.call("plot", c(list(x = x, y = y1, type = "l"), ypar))
        xy$usr <- par("usr")
        if (!missing(y2)) {
            par(new = TRUE)
            if (!is.null(ypar$xlim) && is.null(y2par$xlim)) 
                y2par$xlim <- ypar$xlim
            do.call("plot", c(list(x = x, y = y2, ann = FALSE, 
                bty = "n", xaxt = "n", yaxt = "n", type = "l"), 
                y2par))
            xy$usr2 <- par("usr")
            assign("par.usr2", par("usr"), globalenv())
            if (par("mar")[4] >= 2) 
                axis(4)
            y2par$col <- NULL
            if (par("mar")[4] >= 3) 
                do.call("mtext", c(list(text = labels[3], side = 4, 
                  line = 3, cex = par()$cex), y2par))
            if (!is.null(xlegend) && !is.null(inset)) 
                legend(xlegend, legend = labels[2:3], bty = "o", 
                  lty = lty, lwd = lwd, col = col, inset = inset)
            par(usr = xy$usr)
        }
    }
    else {
        do.call("lines", c(list(x, y1), ypar))
        if (!missing(y2)) {
            if (exists("xy$usr2")) {
                par(usr = xy$usr2)
                assign("par.usr2", par("usr"), globalenv())
            }
            else if (exists("par.usr2")) {
                par(usr = par.usr2)
                xy$usr2 <- par("usr")
            }
            else stop("Error in y2plot: second y axis requires previous call to set it\n", 
                call. = FALSE)
            do.call("lines", c(list(x = x, y = y2), y2par))
            par(usr = opar$usr)
        }
    }
    invisible(xy)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
