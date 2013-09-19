\name{print.summary.sitar}
\alias{print.summary.sitar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{print.summary.sitar}
\description{
A \code{print.summary} method for \code{sitar} objects.}
\usage{
print.summary.sitar(x, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object inheriting class \code{nlme}.}
  \item{verbose}{
a logical to control the amount of output.}
  \item{\dots}{
to specify extra arguments.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A formatted summary of the object.
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
function (x, verbose = FALSE, ...) 
{
    dd <- x$dims
    verbose <- verbose || attr(x, "verbose")
    cat("SITAR nonlinear mixed-effects model fit by ")
    cat(ifelse(x$method == "REML", "REML\n", "maximum likelihood\n"))
    cat("  Call:", deparse(x$call.sitar), fill = TRUE)
    print(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = c(x$logLik), 
        row.names = " "))
    if (verbose) 
        cat("Convergence at iteration:", x$numIter, "\n")
    cat("\n")
    print(summary(x$modelStruct), sigma = x$sigma, reEstimates = x$coef$random, 
        verbose = verbose)
    cat("Fixed effects: ")
    fixF <- x$call$fixed
    if (inherits(fixF, "formula") || is.call(fixF)) {
        cat(deparse(x$call$fixed), "\n")
    }
    else {
        cat(deparse(lapply(fixF, function(el) as.name(deparse(el)))), 
            "\n")
    }
    xtTab <- as.data.frame(x$tTable)
    wchPval <- match("p-value", names(xtTab))
    for (i in names(xtTab)[-wchPval]) {
        xtTab[, i] <- format(zapsmall(xtTab[, i]))
    }
    xtTab[, wchPval] <- format(round(xtTab[, wchPval], 4))
    if (any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 
        0))) {
        levels(xtTab[, wchPval])[wchLv] <- "<.0001"
    }
    row.names(xtTab) <- dimnames(x$tTable)[[1]]
    print(xtTab)
    if (nrow(x$tTable) > 1) {
        corr <- x$corFixed
        class(corr) <- "correlation"
        print(corr, title = " Correlation:", ...)
    }
    cat("\nStandardized Within-Group Residuals:\n")
    print(x$residuals)
    cat("\nNumber of Observations:", x$dims[["N"]])
    cat("\nNumber of Groups: ")
    Ngrps <- dd$ngrps[1:dd$Q]
    if ((lNgrps <- length(Ngrps)) == 1) {
        cat(Ngrps, "\n")
    }
    else {
        sNgrps <- 1:lNgrps
        aux <- rep(names(Ngrps), sNgrps)
        aux <- split(aux, array(rep(sNgrps, lNgrps), c(lNgrps, 
            lNgrps))[!lower.tri(diag(lNgrps))])
        names(Ngrps) <- unlist(lapply(aux, paste, collapse = " \%in\% "))
        cat("\n")
        print(rev(Ngrps))
    }
    invisible(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
