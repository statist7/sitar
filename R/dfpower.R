
#' Tabulate BIC of SITAR models by degrees of freedom and xy power transformations
#'
#' \code{dfpower} fits a series of SITAR models for specified degrees of freedom
#' and power transformations of x and y, returning a table of adjusted BIC.
#'
#' The function provides a convenient way to optimise the model's degrees of freedom
#' and explore transformations of x and y.
#' Adjusted BIC is obtained using \code{\link{BICadj}}, and
#' is set negative for models failing to converge; the
#' run-time can be shortened by reducing \code{maxIter} appropriately. For models
#' failing to fit it returns NA.
#'
#'\code{FUN} can be any function returning a single numerical value.
#'
#' The returned table can be rearranged using \code{\link{aperm}}.
#'
#' @param model fitted \pkg{sitar} model to be updated.
#' @param df vector of degrees of freedom to be fitted (defaults to df in \code{model}).
#' @param xpowers vector of powers to apply to x (defaults to x power in \code{model}).
#' @param ypowers vector of powers to apply to y (defaults to y power in \code{model}).
#' @param FUN function to be tabulated (e.g. BICadj or AICadj).
#' @param maxIter maximum number of iterations per fit.
#' @param verbose logical controlling monitoring.
#' @return 3-way table of returned values by df, x and y powers.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' m1 <- sitar(log(age), height, id, heights, 4)
#' dfpower(m1, 4:5)
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
    expr <- quote(try(update(model, df=df, y=y, x=x, control=control), silent=TRUE))
    if (verbose)
      cat('df', 'iter', 'y ~ x', deparse(substitute(FUN)), '\n')
    for (idf in 1:length(df)) {
      for (iy in 1:length(ypowers)) {
        yp <- dimnames(mat)[[3]][iy]
        for (ix in 1:length(xpowers)) {
          xp <- dimnames(mat)[[2]][ix]
          if (df[idf] == df0 && yp == ylab0 && xp == xlab0)
            . <- model
          else {
            . <- eval(do.call('substitute', list(expr, list(df=df[idf], y=as.lang(yp), x=as.lang(xp)))))
          }
          if (!inherits(., 'try-error')) {
            B <- FUN(.)
            if (.$numIter >= maxIter)
              B <- -B
            mat[idf, ix, iy] <- B
            if (verbose)
              cat(df[idf], .$numIter, yp, '~', xp, B, '\n')
          }
        }
      }
    }
    mat
  }
