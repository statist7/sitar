#' Tabulate BIC of SITAR models by degrees of freedom and xy power transformations
#'
#' \code{dfpower} fits a series of SITAR models for specified degrees of freedom
#' and power transformations of x and y, returning a table of adjusted BIC.
#'
#' The function provides a convenient way to optimise the model's degrees of freedom
#' and explore transformations of x and y.
#' Adjusted BIC is obtained using \code{\link{BICadj}}, and
#' is set negative for models failing to converge; the
#' run-time can be shortened by reducing \code{maxIter} appropriately.
#' The returned table can be rearranged using \code{\link{aperm}}.
#'
#' @param model fitted \pkg{sitar} model to be updated.
#' @param df vector of degrees of freedom to be fitted.
#' @param xpowers vector of powers to apply to x.
#' @param ypowers vector of powers to apply to y.
#' @param maxIter maxium number of iterations.
#' @param verbose logical controlling monitoring.
#' @return 3-way table of BIC by df, x and y powers.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' m1 <- sitar(log(age), height, id, heights, 4)
#' dfpower(m1, 4:6, 1:0)
#' @export dfpower
  dfpower <- function(model, df, xpowers=1, ypowers=xpowers,
                      maxIter=nlmeControl()$maxIter, verbose=FALSE) {
    vars <- function(expr) {
      (av <- all.vars(expr, unique=FALSE))[grep('^pi$', av, invert=TRUE)]
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
    lab <- deparse(substitute(model))
    cs <- model$call.sitar
    xlab0 <- deparse(cs$x)
    xlab <- vars(cs$x)
    ylab0 <- deparse(cs$y)
    ylab <- vars(cs$y)
    df0 <- cs$df
    df <- df[df > 1]
    mat <- array(dim=c(length(df), length(xpowers), length(ypowers)),
                 dimnames=list(df, tidy(xlab, xpowers), tidy(ylab, ypowers)))
    control <- paste0(', control=nlmeControl(maxIter=', maxIter, ', returnObject=TRUE)')
    if (verbose) cat('df', 'iter', 'y ~ x', 'BIC\n')
    for (idf in 1:length(df)) {
      for (iy in 1:length(ypowers)) {
        yp <- dimnames(mat)[[3]][iy]
        for (ix in 1:length(xpowers)) {
          xp <- dimnames(mat)[[2]][ix]
          . <- 'obj'
          if (df[idf] == df0 && yp == ylab0 && xp == xlab0) . <- lab
          else {
            eval(parse(text=paste0("assign('", ., "', update(", lab, ", df=",
                                   df[idf], ", y=", yp, ", x=", xp, control, "))")))
          }
          . <- get(.)
          B <- BICadj(.)
          if (.$numIter >= maxIter) B <- -B
          mat[idf, ix, iy] <- B
          if (verbose) cat(df[idf], .$numIter, yp, '~', xp, B, '\n')
        }
      }
    }
    mat
  }
