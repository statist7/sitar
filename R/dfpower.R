
#' Tabulate BIC of SITAR models by degrees of freedom, fixed effects and xy power transformations
#'
#' \code{dfpower} fits a series of \pkg{sitar} models tabulated by combinations of
#' a) specified degrees of freedom for the spline curve,
#' b) specified fixed effects a, b, c, d,
#' c) specified power transformations of x, and
#' d) specified power transformations of y,
#' returning a four-way array of function values (e.g. BIC) applied to each model.
#' The function provides a convenient way to optimise the model.
#'
#' \code{xpowers} and \code{ypowers} treat power 0 as \code{log}. The formula for
#' \code{x} in \code{object} must be of the form \code{x^power} or \code{fun(x)}, e.g.
#' \code{x}, \code{x^0.5} or \code{log(x)}. More complex formulae e.g. \code{log(x + 1)}
#' will fail. In this case fit the model with the variable \code{x1 = x + 1} instead.
#'
#' \code{FUN} can be any function returning a single numerical value, e.g.
#' \code{BICadj}, \code{BIC}, \code{AIC}, \code{varexp} or \code{sigma}.
#'
#' Other fixed effects in \code{object} for covariates in \code{a.formula}, \code{b.formula},
#' \code{c.formula} or \code{d.formula} are propagated through all the models.
#' This also applies to the \code{control} argument if set in \code{object}.
#'
#' The run-time can be shortened by reducing \code{maxIter},
#' as models often converge quickly or not at all.
#'
#' @param object fitted \pkg{sitar} model to be updated.
#' @param df vector of integer spline degrees of freedom to be fitted (defaults to \code{df} in \code{object}).
#' @param fixed character vector of fixed effects to be included
#' (defaults to \code{fixed} in \code{object}, typically 'a + b + c').
#' @param xpowers vector of powers to apply to x (defaults to the power of x in \code{object}).
#' @param ypowers vector of powers to apply to y (defaults to the power of y in \code{object}).
#' @param FUN function to be tabulated (default \code{BICadj}).
#' @param maxIter maximum number of iterations per fit (default \code{nlmeControl()$maxIter}).
#' @param drop logical which if TRUE (default) drops redundant dimensions and labels from the returned array.
#' @param verbose logical controlling monitoring, which gives \code{numIter} for each model.
#'
#' @return Four-way array of returned values, ranked with the largest dimensions first,
#' and by default with single-level dimensions dropped.
#'
#' Values are returned with changed sign if the model fit generates a warning, or as
#' NA if there is an error.
#'
#' @seealso \code{\link{aperm}} transposes the returned array;
#' \code{\link{addmargins}} adds margins.
#'
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' m1 <- sitar(x = age, y = height, id = id, data = heights, df = 4)
#' \dontshow{
#' dfpower(m1, df = 4:5, fixed = 'a+c', xpowers = 0, maxIter = 4)
#' }
#' \donttest{
#' dfpower(m1, df = 4:6, fixed = c('a', 'a+b', 'a+c', 'a+b+c'),
#'   xpowers = 0:1, ypowers = 0:1, maxIter = 8)
#' }
#' @export dfpower
#' @md
dfpower <- function(object, df, fixed, xpowers, ypowers, FUN=BICadj,
                    maxIter=50, drop = TRUE, verbose=FALSE) {

  vars <- function(expr) {
    (av <- all.vars(expr))[grep('^pi$', av, invert=TRUE)]
  }

  power.tidy <- function(lab, p) {
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

  fixed.tidy <- function(f) {
    unlist(lapply(f, function(x)
    paste0(all.names(str2expression(x), functions = FALSE, unique = TRUE),
           collapse = '+')))
  }

  stopifnot('object not a sitar model' = inherits(object, 'sitar'))
  lab <- deparse(substitute(object))
  cs <- object$call.sitar
  df0 <- cs$df
  fixed0 <- names(fixef(object))
  fixed0 <- fixed0[grep('^[abcd]$', fixed0)]
  d_only <- identical(sort(fixed0), c('a', 'd'))
  d_plus <- d_only && length(fixed0) > 2L
  fixed0 <- paste(fixed0, collapse = '+')
  xlab0 <- deparse(cs$x)
  xlab <- vars(cs$x)
  ylab0 <- deparse(cs$y)
  ylab <- vars(cs$y)
  if (missing(df))
    df <- df0
  df <- as.integer(unique(df[df > 0]))
  if (d_only)
    df <- 0L
  if (missing(fixed))
    fixed <- fixed0
  fixed <- fixed.tidy(fixed)
  if (missing(xpowers))
    xpowers <- getL(cs$x)
  xpowers <- xpowers[!is.na(xpowers)]
  if (missing(ypowers))
    ypowers <- getL(cs$y)
  ypowers <- ypowers[!is.na(ypowers)]
  if (is.null(cs$control))
    control <- nlmeControl(maxIter=maxIter, returnObject=TRUE)
  else {
    control <- as.list(cs$control)
    control$maxIter <- maxIter
    control$returnObject <- TRUE
    control <- as.call(control)
  }
  stopifnot('df not non-negative integer(s)' = all(df >= 0L) & length(df) > 0L,
            'fixed not character(s)' = is.character(fixed) & length(fixed) > 0L,
            'no valid xpowers' = length(xpowers) > 0L,
            'no valid ypowers' = length(ypowers) > 0L,
            'verbose not logical' = is.logical(verbose))
  mat <- array(dim=c(length(df), length(fixed), length(xpowers), length(ypowers)),
               dimnames=list(df, fixed, power.tidy(xlab, xpowers), power.tidy(ylab, ypowers)))
  expr <- quote(runWarn(update(object, fixed=fixed, df=df, y=y, x=x, control=control)))
  if (verbose)
    cat('df', 'iter', 'fixed', 'y ~ x', deparse(substitute(FUN)), '\n')
  for (idf in seq_along(df)) {
    dft <- df[idf]
    for (iy in seq_along(ypowers)) {
      yp <- dimnames(mat)[[4]][iy]
      for (ix in seq_along(xpowers)) {
        xp <- dimnames(mat)[[3]][ix]
        for (ifix in seq_along(fixed)) {
          fixt <- fixed[ifix]
          . <- eval(do.call('substitute', list(
            expr, list(df=dft, fixed=fixt, y=as.lang(yp), x=as.lang(xp)))))
          if (!inherits(., 'error')) {
            B <- FUN(.$value)
            if (length(.$warnings) > 0L) {
              signB <- sign(B)
              B <- -B
            }
            mat[idf, ifix, ix, iy] <- B
            if (verbose)
              cat(dft, .$value$numIter, fixt, yp, '~', xp, B, '\n')
          }
        }
      }
    }
  }
  if (exists('signB'))
    cat(paste('\n**', ifelse(signB > 0, 'negative', 'positive'),
              'values indicate the model fitted with warnings **\n\n'))
  mat <- aperm(mat, rev(order(dim(mat))))
  if (drop)
    mat <- drop(mat)
  mat
}
