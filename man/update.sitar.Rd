\name{update.sitar}
\alias{update.sitar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{update.sitar}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update.sitar(object, ..., evaluate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{evaluate}{
%%     ~~Describe \code{evaluate} here~~
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
function (object, ..., evaluate = TRUE) 
{
    call <- object$call.sitar
    if (is.null(call)) 
        stop("need an object with call.sitar component")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras) > 0) {
        all.pars <- c(as.list(call), formals(sitar))[-1]
        for (a in letters[1:3]) {
            pos.e <- pmatch(paste(a, "form", sep = "."), names(extras), 
                0)
            if (pos.e) {
                pos.p <- pmatch(paste(a, "formula", sep = "."), 
                  names(all.pars), 0)
                if (pos.p) 
                  extras[[pos.e]] <- update.formula(all.pars[[pos.p]], 
                    extras[[pos.e]])
            }
        }
        existing <- pmatch(names(extras), names(call))
        if (sum(existing, na.rm = TRUE)) 
            for (a in 1:length(existing)) call[existing[a]] <- extras[a]
        existing <- !is.na(existing)
        if (any(!existing)) 
            call <- as.call(c(as.list(call), extras[!existing]))
    }
    if (sum(pmatch(names(extras), c("x", "y", "id", "data", "fixed", 
        "random", "a.formula", "b.formula", "c.formula", "start", 
        "bstart", "xoffset")), na.rm = TRUE) == 0) {
        start. <- list(fixed = fixef(object), random = ranef(object))
        data <- eval(call$data)
        id <- factor(eval(call$id, data))
        if (nlevels(id) != dim(ranef(object))[1]) {
            idcheck <- rownames(ranef(object)) \%in\% levels(id)
            start.$random <- ranef(object)[idcheck, ]
            cat(dim(ranef(object))[1] - sum(idcheck), "subjects omitted\n")
            newid <- !levels(id) \%in\% rownames(ranef(object))
            if (sum(newid) > 0) {
                newre <- matrix(0, nrow = sum(newid), ncol = dim(ranef(object))[2], 
                  dimnames = list(levels(id)[newid], dimnames(ranef(object))[[2]]))
                start.$random <- rbind(start.$random, newre)
                cat(sum(newid), "subjects added\n")
            }
        }
        if (sum(pmatch(names(extras), c("df", "knots", "bounds")), 
            na.rm = TRUE) > 0) {
            x <- eval(call$x, data)
            xoffset <- object$xoffset
            if (is.null(xoffset)) 
                xoffset <- mean(x)
            knots <- attr(object$ns$model$ns, "knots") + xoffset
            df <- object$ns$rank - 1
            bounds <- attr(object$ns$model$ns, "Boundary.knots") + 
                xoffset
            if (length(fixef(object)) > df + 1) 
                fixed.extra <- (df + 2):length(fixef(object))
            else fixed.extra <- NULL
            if (!is.null(extras$knots)) {
                knots <- extras$knots
                df <- length(knots) + 1
            }
            else if (!is.null(extras$df)) {
                df <- extras$df
                knots <- quantile(x, (1:(df - 1))/df)
            }
            if (!is.null(extras$bounds)) {
                bounds <- extras$bounds
                if (length(bounds) == 1) 
                  bounds <- range(x) + abs(bounds) * c(-1, 1) * 
                    diff(range(x))
            }
            b <- fixef(object)["b"]
            b <- ifelse(!is.na(b), b, 0)
            c <- fixef(object)["c"]
            c <- ifelse(!is.na(c), c, 0)
            xc <- (x - xoffset - b) * exp(c)
            knots <- knots - xoffset
            bounds <- bounds - xoffset
            spline.lm <- lm(fitted(object, level = 0) ~ ns(xc, 
                knots = knots, Bound = bounds))
            start.$fixed <- c(coef(spline.lm)[c(2:(df + 1), 1)], 
                fixef(object)[fixed.extra])
        }
        assign("start.", start., parent.frame())
        if (!"start" \%in\% names(call)) 
            call <- as.call(c(as.list(call), start = quote(start.)))
    }
    else call$start <- NULL
    if (evaluate) 
        eval(call, parent.frame())
    else call
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
