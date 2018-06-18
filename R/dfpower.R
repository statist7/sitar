
#' Tabulate BIC of SITAR models by degrees of freedom and xy power transformations
#'
#' \code{dfpower} fits a series of SITAR models tabulated by specified degrees of freedom
#' and power transformations of x and y, returning a table of function values (e.g. BIC)
#' applied to each model.
#'
#' The function provides a convenient way to optimise the model's degrees of freedom
#' and explore transformations of x and y, based by default on adjusted BIC.
#' The function value is returned with changed sign if there is a warning, or as
#' NA if there is an error. The run-time can be shortened by reducing \code{maxIter},
#' as the models that converge do so in relatively few iterations, and much of
#' the run-time is spent on models that fail to converge.
#'
#'\code{FUN} can be any function returning a single numerical value.
#'
#' The returned table can be rearranged using \code{\link{aperm}}.
#'
#' @param model fitted \pkg{sitar} model to be updated.
#' @param df vector of degrees of freedom to be fitted (defaults to df in \code{model}).
#' @param xpowers vector of powers to apply to x (defaults to x power in \code{model}).
#' @param ypowers vector of powers to apply to y (defaults to y power in \code{model}).
#' @param FUN function to be tabulated (default BICadj, or e.g. AICadj or varexp).
#' @param maxIter maximum number of iterations per fit.
#' @param verbose logical controlling monitoring.
#' @return Table or vector of returned values.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' m1 <- sitar(log(age), height, id, heights, 4)
#' dfpower(m1, df=4:5, xpowers=0:1, maxIter=4)
#' @export dfpower
dfpower <- function(model, df, xpowers, ypowers, FUN=BICadj,
                    maxIter=nlmeControl()$maxIter, verbose=FALSE) {

  vars <- function(expr) {
    (av <- all.vars(expr))[grep('^pi$', av, invert=TRUE)]
  }

  tidy <- function(lab, p) {
    f <- matrix(c(-1, 0, 0.5, 1, Inf,
                  '1/', 'log(', 'sqrt(', '', '',
                  '', ')', ')', '', '^'), ncol=3)
    . <- match(p, as.numeric(f[, 1]), nrow(f))
    res <- paste0(f[., 2], lab, f[., 3])
    res[. == nrow(f)] <- paste0(res, p)[. == nrow(f)]
    res
  }

  trapWarn <- function(expr) {
    warnings <- NULL
    wHandler <- function(w) {
      warnings <<- c(warnings, list(w))
      invokeRestart("muffleWarning")
    }
    val <- withCallingHandlers(expr, warning = wHandler)
    list(value = val, warnings = warnings)
  }

  runWarn <- function(expr)
    # returns value and (list of any) warnings, else error message
    tryCatch(trapWarn(expr), warning=function(w) {},
             error=function(e) e)

  as.lang <- function(x) parse(text=x)[[1]]

  lab <- deparse(substitute(model))
  cs <- model$call.sitar
  xlab0 <- deparse(cs$x)
  xlab <- vars(cs$x)
  ylab0 <- deparse(cs$y)
  ylab <- vars(cs$y)
  df0 <- cs$df
  if (missing(df)) df <- df0
  df <- unique(df[df > 1])
  if (missing(xpowers)) xpowers <- getL(cs$x)
  if (missing(ypowers)) ypowers <- getL(cs$y)
  mat <- array(dim=c(length(df), length(xpowers), length(ypowers)),
               dimnames=list(df, tidy(xlab, xpowers), tidy(ylab, ypowers)))
  control <- nlmeControl(maxIter=maxIter, returnObject=TRUE)
  expr <- quote(runWarn(update(model, df=df, y=y, x=x, control=control)))
  if (verbose)
    cat('df', 'iter', 'y ~ x', deparse(substitute(FUN)), '\n')
  for (idf in seq_along(df)) {
    for (iy in seq_along(ypowers)) {
      yp <- dimnames(mat)[[3]][iy]
      for (ix in seq_along(xpowers)) {
        xp <- dimnames(mat)[[2]][ix]
        if (df[idf] == df0 && yp == ylab0 && xp == xlab0) {
          . <- list()
          .$value <- model
        }
        else {
          . <- eval(do.call('substitute', list(expr, list(df=df[idf], y=as.lang(yp), x=as.lang(xp)))))
        }
        if (!inherits(., 'error')) {
          B <- FUN(.$value)
          if (length(.$warnings))
            B <- -B
          mat[idf, ix, iy] <- B
          if (verbose)
            cat(df[idf], .$value$numIter, yp, '~', xp, B, '\n')
        }
      }
    }
  }
  drop(mat)
}
