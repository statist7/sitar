\name{sitar}
\alias{sitar}
\alias{update.sitar}
% \alias{summary.sitar}
% \alias{print.sitar}
% \alias{plot.sitar}
% \alias{print.summary.sitar}
% \alias{lines.sitar}
\title{SITAR growth curve analysis}
\description{SITAR is a method of growth curve analysis, based on \code{nlme}, that takes a set of growth curves and estimates
 a mean growth curve as a regression spline, plus a set of up to three fixed and
 random effects (a, b and c) defining how individual growth curves differ from the mean
 curve. Effect a (or alpha) measures size, and is a random intercept relative to the
 spline curve intercept. Effect b (or beta) measures tempo, the timing of the growth
 process, and reflects a shift on the x scale relative to the mean. Effect c (or gamma)
 is velocity, and indicates how the x scale is stretched or shrunk reflecting the rate
 at which 'time' passes for individuals. The aim is for individual curves, adjusted for
 a, b and c, to be superimposed on the mean curve. SITAR means SuperImposition by Translation And Rotation.}
 \usage{sitar(x, y, id, data, df, knots, fixed = random, random = "a+b+c", a.formula = ~1,
 b.formula = ~1, c.formula = ~1, bounds = 0.04, start, bstart = "mean", xoffset = "mean",
 file = ".sitar.temp.R", returndata=FALSE, verbose = FALSE, correlation = NULL, weights = NULL, subset = NULL, method = "ML", na.action = na.fail, control = nlmeControl(returnObject = TRUE), newform = TRUE)
 
\method{update}{sitar}(object, ..., evaluate = TRUE)
}
\arguments{
  \item{x}{vector of ages.}
  \item{y}{vector of measurements.}
  \item{id}{factor of subject identifiers.}
  \item{data}{optional data frame containing variables \code{x}, \code{y} and \code{id}.}
  \item{df}{degrees of freedom for regression spline.}
  \item{knots}{vector of values for knots (default df quantiles of \code{x} distribution).}
  \item{fixed}{character string defining which of a, b and c are fixed effects (default random).}
  \item{random}{character string defining which of a, b and c are random effects (default "a+b+c").}
  \item{a.formula}{formula for fixed effect a (default ~ 1).}
  \item{b.formula}{formula for fixed effect b (default ~ 1).}
  \item{c.formula}{formula for fixed effect c (default ~ 1).}
  \item{bounds}{span of \code{x} for regression spline, or extension (default range(\code{x}) +/- 0.04).}
  \item{start}{optional numeric vector, or list of initial estimates for the fixed effects and random effects.}
  \item{bstart}{optional starting value for fixed effect b (either "mean" (default), "apv" or value).}
  \item{xoffset}{optional value of offset for \code{x} (either "mean" (default), "apv" or value).}
  \item{file}{filename for intermediate code (default .sitar.temp.R).}
  \item{returndata}{logical which if TRUE causes the working matrix to be returned, or if FALSE (default) returns the fitted model. Setting returndata TRUE is useful in conjunction with \code{subset} and \code{subsample} for simulation purposes.}
  \item{verbose}{optional logical value to print information on the evolution of the iterative algorithm (default FALSE).}
  \item{correlation}{optional \code{corStruct} object describing the within-group correlation structure (see NLME documentation).}
  \item{weights}{optional \code{varFunc} object or one-sided formula describing the within-group heteroscedasticity structure (see \code{\link{nlme}}).}
  \item{subset}{optional expression indicating the subset of the rows of data that should be used in the fit (see \code{\link{nlme}}).}
  \item{method}{character string, either "REML" or "ML" (default) (see \code{\link{nlme}}).}
  \item{na.action}{function for when the data contain NAs (see \code{\link{nlme}}).}
  \item{control}{list of control values for the estimation algorithm (see \code{\link{nlme}}).}
  \item{newform}{logical to distinguish between old and new forms of \code{sitar} object - default TRUE indicates new form.}
  \item{object}{object of class \code{sitar}.}
  \item{\dots}{further parameters for \code{update} consisting of any of the parameters for \code{sitar}.}
  \item{evaluate}{logical to control evaluation.  If TRUE (default) the expanded \code{update} call is passed to \code{sitar} for evaluation, while if FALSE the expanded call itself is returned.}
}
\value{
An object of class \code{sitar}, based on \code{nlme}, representing the nonlinear mixed-effects model fit. Generic functions such as \code{print}, \code{plot} and \code{summary} have methods to show the results of the fit. The functions \code{resid}, \code{coef}, \code{fitted}, \code{fixed.effects}, and \code{random.effects} can be used to extract some of its components. \code{update} updates the model.
}
\references{
The idea of SITAR growth curve analysis arose from the paper by Beath (2007) and was first described in Cole et al (2010). The other references describe applications of SITAR to a variety of data forms.
Beath KJ. Infant growth modelling using a shape invariant model with random efffects. Stat Med 2007;26:2547-64.
Cole TJ, Cortina Borja M, Sandhu J, et al. Nonlinear growth generates age changes in the moments of the frequency distribution: the example of height in puberty. Biostatistics 2008;9:159-71.
Cole TJ, Donaldson MD, Ben-Shlomo Y. SITAR--a useful instrument for growth curve analysis. Int J Epidemiol 2010;39:1558-66.
Gault EJ, Perry RJ, Cole TJ, et al. Effect of oxandrolone and timing of pubertal induction on final height in Turner's syndrome: randomised, double blind, placebo controlled trial. BMJ 2011;342:d1980.
Johnson L, Llewellyn CH, van Jaarsveld CHM, et al. Genetic and environmental influences on infant growth: prospective analysis of the Gemini twin birth cohort. PLoS ONE 2011;6:e19918.
Prentice A, Dibba B, Sawo Y, et al. The effect of prepubertal calcium carbonate supplementation on the age of peak height velocity in Gambian adolescents. Am J Clin Nutr 2012;96:1042-50.
Dean MC, Cole TJ. Human life history evolution explains dissociation between the timing of tooth eruption and peak rates of root growth. PLoS ONE 2013;8:e54534.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y, id, data, df, knots, fixed = random, random = "a+b+c", 
    a.formula = ~1, b.formula = ~1, c.formula = ~1, bounds = 0.04, 
    start, bstart = 0, xoffset = "mean", file = ".sitar.temp.R", 
    verbose = FALSE, correlation = NULL, weights = NULL, method = "ML", 
    na.action = na.fail, control = nlmeControl(returnObject = TRUE)) 
{
    if (!missing(data)) {
        if (!deparse(substitute(data), width = 99) \%in\% search()) {
            on.exit(detach(data))
            attach(data)
        }
    }
    if (missing(df) & missing(knots)) 
        stop("either df or knots must be specified")
    if (!missing(df) & !missing(knots)) 
        cat("both df and knots specified - df redefined from knots\n")
    if (missing(knots)) 
        knots <- quantile(na.action(x), (1:(df - 1))/df)
    else df <- length(knots) + 1
    if (length(bounds) == 1) 
        bounds <- range(x) + abs(bounds) * c(-1, 1) * diff(range(x))
    if (length(bounds) != 2) 
        stop("bounds should be length 1 or 2")
    if (xoffset == "mean") 
        xoffset <- mean(x)
    else if (xoffset == "apv") {
        spline.lm <- lm(na.action(y) ~ ns(x, knots = knots, Bound = bounds))
        xoffset <- makess(x, fitted(spline.lm))$apv[1]
    }
    else if (!is.numeric(xoffset)) 
        stop("xoffset must be \"mean\" or \"apv\" or numeric")
    xc <- x - xoffset
    knots <- knots - xoffset
    bounds <- bounds - xoffset
    spline.lm <- lm(y ~ ns(xc, knots = knots, Bound = bounds))
    if (missing(start)) 
        start <- coef(spline.lm)[c(2:(df + 1), 1)]
    fix <- fixed
    if (!grepl("a", fix)) 
        fix <- paste("a", fix, sep = "+")
    pars <- fixed <- ss <- paste("s", 1:df, sep = "")
    names(model) <- model <- letters[1:3]
    mm.formula <- as.formula("~ 1")
    for (l in model) {
        if (!grepl(l, fix) && !grepl(l, random)) {
            model[l] <- NA
            next
        }
        pars <- c(pars, l)
        formula <- get(paste(l, "formula", sep = "."))
        if (formula == as.formula("~ -1") || formula == as.formula("~ 1-1") || 
            !grepl(l, fix)) 
            next
        if (formula == as.formula("~ 1")) {
            fixed <- c(fixed, l)
            if (!is.list(start)) {
                if (l != "a") 
                  start <- c(start, 0)
            }
        }
        else {
            if (formula != mm.formula) {
                mm <- model.matrix(formula, data)
                if (dim(mm)[[1]] < length(y)) 
                  stop("Missing values in data")
                mm.formula <- formula
                colnames(mm) <- gsub(" ", "_", colnames(mm))
                colnames(mm) <- gsub(":", ".", colnames(mm), 
                  fixed = TRUE)
                mm.intercept <- colnames(mm)[1] == "(Intercept)"
                cmm <- apply(mm, 2, function(x) sd(x, na.rm = TRUE) > 
                  0)
                mm.names <- colnames(mm)[cmm]
                mm <- as.matrix(mm[, mm.names])
                mm <- scale(mm, scale = FALSE)
                colnames(mm) <- mm.names
            }
            if (mm.intercept) {
                fixed <- c(fixed, l)
                if (!is.list(start)) {
                  if (l != "a") 
                    start <- c(start, 0)
                }
            }
            for (i in 1:dim(mm)[[2]]) {
                var <- colnames(mm)[i]
                if (var == "(Intercept)") 
                  next
                rc <- paste(l, var, sep = ".")
                pars <- c(pars, rc)
                fixed <- c(fixed, rc)
                if (!is.list(start)) 
                  start <- c(start, 0)
                model[l] <- paste(model[l], "+", rc, "*", var, 
                  sep = "")
                if (!var \%in\% pars) {
                  pars <- c(pars, var)
                  assign(var, mm[, i])
                }
            }
        }
    }
    pars <- paste(pars, collapse = ",")
    fixed <- paste(fixed, collapse = "+")
    sscomma <- paste(ss, collapse = ",")
    nsd <- paste(model["a"], "+")
    nsf <- paste("(xc", ifelse(!is.na(model["b"]), paste("- (", 
        model["b"], "))"), ")"))
    if (!is.na(model["c"])) 
        nsf <- paste(nsf, "* exp(", model["c"], ")")
    if (verbose) {
        cat("df", df, "xoffset", xoffset, "\nknots\n", knots, 
            "\nbounds\n", bounds)
        if (is.list(start)) {
            cat("\nstarting values\n  fixed effects\n", start$fixed)
            if (!is.null(start$random)) {
                cat("\n  random effects\n")
                print(start$random)
            }
        }
        else cat("\nstarting values\n", start, "\n")
    }
    code <- c("\tfitnlme <- function(xc,$pars) {", "#\t\t** created by sitar - edits will be lost **", 
        "\t\tsplinecoefs <- as.matrix(cbind($sscomma))", "\t\tas.vector( $nsd", 
        "\t\t(splinecoefs * as.matrix(ns($nsf,", "\t\t\tknots=knots, Boundary.knots=bounds))) \%*\% ", 
        "\t\t\tmatrix(rep(1,df), ncol=1))", "\t}", "\tassign('fitnlme', fitnlme, globalenv())", 
        "\tnlme.out <- nlme(y ~ fitnlme(xc, $pars),", "\tfixed = $fixed ~ 1,", 
        "\trandom = $random ~ 1 | id,", "\tstart = start,", "#\tgroups, subset, (unimplemented nlme arguments)", 
        "\tcorrelation = correlation, weights = weights, method = method,", 
        "\tna.action = na.action, control = control, verbose = verbose)")
    for (i in c("random", "pars", "fixed", "sscomma", "nsd", 
        "nsf", "file")) {
        code <- gsub(paste("$", i, sep = ""), get(i), code, fixed = TRUE)
    }
    write(code, file = file, sep = "\n")
    if (verbose) {
        cat("file", file, "\n")
        cat(code, sep = "\n")
    }
    elapsed <- system.time(source(file, local = TRUE))
    nlme.out$call.sitar <- match.call()
    nlme.out$xoffset <- xoffset
    nlme.out$ns <- spline.lm
    nlme.out$elapsed <- elapsed
    if (class(nlme.out)[1] == "nlme") 
        class(nlme.out) <- c("sitar", class(nlme.out))
    nlme.out
  }
}
\keyword{ package }
\keyword{ nonlinear }
\keyword{ regression }
\keyword{ models }
