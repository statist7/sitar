\name{plot.sitar}
\alias{plot.sitar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.sitar}
\description{
A version of \code{plot} for objects of class \code{sitar}, providing various flavours of plot of the fitted growth curves. }
\usage{
plot.sitar(model, opt = "dv", labels, apv = FALSE, xfun, yfun, subset = NULL, abc = c(a = 0, b = 0, c = 0), add = FALSE, xy = NULL, nlme = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
an object of class \code{sitar}.}
  \item{opt}{
character string containing a subset of letters corresponding to the options: 'd' for fitted Distance curve, 'v' for fitted Velocity curve, 'f' for individual Fitted growth curves, 'u' for Unadjusted individual growth curves, and 'a' for Adjusted individual growth curves. Options 'dvf' give spline curves, while 'ua' give curves made up as line segments. If both 'd' and 'v' are specified, the axis for the velocity curve appears on the right side of the plot (y2), and legend identifying the distance and velocity curves is provided.}
  \item{labels}{
optional character vector containing plot labels for \code{x}, \code{y} and \code{y} velocity from \code{model}. The first two elements can alternatively be provided via \code{\link{par}} parameters \code{xlab} and \code{ylab}, and the third element via \code{y2par} (see Details). Default labels are the names of \code{x} and \code{y}, and "\code{y} velocity".}
  \item{apv}{
optional logical specifying whether or not to calculate the age at peak velocity from the velocity curve. If TRUE, age at peak velocity is calculated as the age when the second derivative of the fitted curve changes sign (if necessary after applying \code{xfun} and \code{yfun}). Age at peak velocity is marked in the plot with a vertical dotted line, and its value, along with peak velocity, is returned.}
  \item{xfun}{
optional function to be applied to \code{x} prior to plotting.}
  \item{yfun}{
optional function to be applied to \code{y} prior to plotting.}
  \item{subset}{
optional logical vector of length \code{x} defining a subset of data rows to be plotted.}
  \item{abc}{
a character string matching an /code{id} level in \code{ranef(model)}, or else a list of up to three values labelled a, b and c, e.g. abc=list(a=1, b=0.1, c=-0.1), or where the values of a, b and c are random effects defining an individual growth curve from the average curve.}
  \item{add}{
optional logical defining whether the plot is pre-existing (TRUE) or new (FALSE). Alternatively use \code{lines}.}
  \item{xy}{
defunct.}
  \item{nlme}{
optional logical which set TRUE plots the model as an \code{nlme} object, using \code{\link{plot.nlme}} arguments.}
  \item{\dots}{
Further graphical parameters (see \code{par}) may also be supplied as arguments, e.g. axis labels \code{xlab} and \code{ylab}, line type, lty, line width, lwd, and color, col. For the velocity (y2) plot \code{y2par} can be used (see Details).}
\details{
For option 'dv' (the default) the velocity curve plot (with right axis) can be annotated with \code{par} parameters given as a named list called \code{y2par}.}
\value{
Returns a list of three objects:
\item{ss}{the \code{smooth.spline} object corresponding to the fitted distance curve.}
\item{usr}{the value of \code{par('usr')}.}
\item{usr2}{the value of \code{par('usr')} for the velocity curve.}
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
\code{\link{lines.sitar}}, \code{\link{mplot}}, \code{\link{plot.clean}}, \code{\link{y2plot}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (model, opt = "dv", labels, apv = FALSE, xfun, yfun, 
    subset = NULL, abc = c(a = 0, b = 0, c = 0), newplot = TRUE, 
    xy = NULL, nlme = FALSE, ...) 
{
    if (nlme) {
        call <- match.call()[-1]
        names(call)[[1]] <- "x"
        do.call("plot.lme", as.list(call))
    }
    else {
        call <- model$call.sitar
        cnames <- names(call)
        mnames <- c("x", "y", "id", "data")
        cnames <- cnames[match(mnames, cnames, 0)]
        mcall <- call[cnames]
        data <- eval(mcall$data)
        attach(data)
        on.exit(detach(data))
        x <- eval(mcall$x)
        y <- eval(mcall$y)
        id <- eval(mcall$id)
        nf <- length(fitted(model))
        if (nf != length(y)) 
            stop(paste("model (length=", nf, ") incompatible with data (rows=", 
                length(y), ")", sep = ""))
        if (is.null(subset)) 
            subset <- rep(TRUE, nf)
        ARG <- list(...)
        if (!"xlab" \%in\% names(ARG)) {
            xl <- ifelse(missing(labels), deparse(mcall$x), labels[1])
            ARG <- c(ARG, list(xlab = xl))
        }
        else xl <- ARG$xlab
        if (!"ylab" \%in\% names(ARG)) {
            yl <- ifelse(missing(labels), deparse(mcall$y), labels[2])
            ARG <- c(ARG, list(ylab = yl))
        }
        else yl <- ARG$ylab
        if (grepl("d", opt) || grepl("v", opt) || apv) {
            xt <- x[subset]
            yt <- fitted(model, level = 0)[subset]
            if (!missing(abc)) {
                if (!is.null(names(abc))) {
                  random <- dimnames(ranef(model))[[2]]
                  for (l in letters[1:3]) if (is.na(abc[l]) || 
                    !l \%in\% random) 
                    abc[l] <- 0
                }
                else if (length(abc) == 1) {
                  idabc <- dimnames(ranef(model))[[1]] \%in\% abc
                  if (sum(idabc) == 0) 
                    stop(paste("id", abc, "not found"))
                  abc <- ranef(model)[idabc, ]
                }
                else stop("abc should be either single id level or up to three named random effect values")
                xoffset <- model$xoffset
                if (is.null(xoffset)) 
                  xoffset <- mean(x)
                xt <- (xt - xoffset + abc[["b"]])/exp(abc[["c"]]) + 
                  xoffset
                yt <- yt + abc[["a"]]
            }
            xy$ss <- ss <- makess(xt, yt, xfun = xfun, yfun = yfun)
            ss1 <- predict(ss, deriv = 1)
            ss2 <- predict(ss, deriv = 2)
            if (missing(labels)) 
                labels <- c(xl, yl, paste(yl, "velocity"))
            if (grepl("d", opt) && grepl("v", opt)) {
                xy <- do.call("y2plot", c(list(x = ss$x, y1 = ss$y, 
                  y2 = ss1$y, labels = labels, newplot = newplot, 
                  xy = xy), ARG))
                newplot <- FALSE
            }
            else if (grepl("d", opt)) {
                xy <- do.call("y2plot", c(list(x = quote(ss$x), 
                  y1 = quote(ss$y), newplot = newplot, xy = xy), 
                  ARG))
                newplot <- FALSE
            }
            else if (grepl("v", opt)) {
                ARG$ylab <- labels[3]
                xy <- do.call("y2plot", c(list(x = quote(ss$x), 
                  y1 = quote(ss1$y), newplot = newplot, xy = xy), 
                  ARG))
                newplot <- FALSE
            }
        }
        if (grepl("u", opt)) {
            if (!missing(xfun)) 
                x <- xfun(x)
            if (!missing(yfun)) 
                y <- yfun(y)
            do.call("mplot", c(list(x = quote(x), y = quote(y), 
                id = quote(id), subset = subset, newplot = newplot), 
                ARG))
            newplot <- FALSE
        }
        if (grepl("f", opt)) {
            y <- fitted(model, level = 1)
            if (!missing(xfun)) 
                x <- xfun(x)
            if (!missing(yfun)) 
                y <- yfun(y)
            do.call("mplot", c(list(x = quote(x), y = quote(y), 
                id = quote(id), subset = subset, newplot = newplot), 
                ARG))
            newplot <- FALSE
        }
        if (grepl("a", opt)) {
            fred <- summary(model)
            x <- fred$x.adj
            y <- fred$y.adj
            if (!missing(xfun)) 
                x <- xfun(x)
            if (!missing(yfun)) 
                y <- yfun(y)
            do.call("mplot", c(list(x = quote(x), y = quote(y), 
                id = quote(id), subset = subset, newplot = newplot), 
                ARG))
            newplot <- FALSE
        }
        if (apv) {
            xy$apv <- ss$apv
            if (!is.na(opt)) 
                print(signif(xy$apv, 4))
            if (!newplot) {
                if (is.null(ARG$y2par$lty)) 
                  ARG$y2par$lty <- 3
                do.call("abline", c(list(v = xy$apv[1]), ARG$y2par))
            }
        }
        invisible(xy)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
