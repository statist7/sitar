\name{plot.sitar}
\alias{plot.sitar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot.sitar}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.sitar(model, opt = "dv", labels, apv = FALSE, xfun, yfun, subset = NULL, abc = c(a = 0, b = 0, c = 0), newplot = TRUE, xy = NULL, nlme = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{opt}{
%%     ~~Describe \code{opt} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{apv}{
%%     ~~Describe \code{apv} here~~
}
  \item{xfun}{
%%     ~~Describe \code{xfun} here~~
}
  \item{yfun}{
%%     ~~Describe \code{yfun} here~~
}
  \item{subset}{
%%     ~~Describe \code{subset} here~~
}
  \item{abc}{
%%     ~~Describe \code{abc} here~~
}
  \item{newplot}{
%%     ~~Describe \code{newplot} here~~
}
  \item{xy}{
%%     ~~Describe \code{xy} here~~
}
  \item{nlme}{
%%     ~~Describe \code{nlme} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
