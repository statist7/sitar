#' Invert an expression defining a data transformation
#'
#' Given a transformed variable and the expression used to transform it, \code{ifun} creates
#' a function containing the inverse expression that will back-transform the variable.
#'
#' \code{ifun} returns the inverting function such that
#' \code{ifun(expr)(eval(expr)) = varname}, where
#' \code{expr} can include any of the invertible functions in the 'Math' and 'Ops'
#' groups.
#'
#' To illustrate its use, consider variants of the \code{sitar} model
#' \code{height ~ age} where \code{age} and/or \code{height} are transformed,
#' e.g. \code{height ~ log(age)} or \code{log(height) ~ sqrt(age)}. Each model
#' is of the form \code{y ~ x} but the units of \code{x} and \code{y} vary.
#'
#' The models are compared by plotting the fitted curves in their original units,
#' by first applying suitable functions to back-transform \code{x} and \code{y}.
#' For example with \code{log(age)}, where \code{expr = quote(log(age))},
#' the function \code{ifun = function(x) exp(x)} back-transforms
#' \code{eval(expr)} to give \code{age}. See the first example.
#'
#' \code{ifun} generalises this process for increasingly complex \code{expr}, as
#' the next two examples show.
#'
#' The final example shows \code{ifun} in action with \code{\link{plot.sitar}},
#' which uses \code{ifun} as the default function for arguments \code{xfun} and
#' \code{yfun} - they are used to back-transform \code{x} and \code{y} using the
#' values of \code{expr} for \code{x} and \code{y} extracted from the model's
#' \code{sitar} call.
#'
#' Structuring \code{expr} suitably ensures it can be inverted - it should
#' contain a single mention of a single variable (\code{varname} here),
#' and possibly functions such as \eqn{f(.)}, \eqn{g(.)}, \eqn{h(.)} etc
#' such that \code{expr} = \eqn{f(g(h((varname))))}. The number of such functions
#' is in principle unlimited. \code{ifun} returns \code{function(x)}
#' \eqn{h^{-1}(g^{-1}(f^{-1}((x))))}{h^-1(g^-1(f^-1((x))))},
#' which ensures that
#' \code{expr} is invertible so long as the individual functions are invertible.
#'
#' @param expr a single-variable call or quoted expression to be inverted.
#' The variable's name in \code{expr} is referred to here as \code{varname}.
#' @param verbose a logical controlling printing of the intermediate functions
#' \eqn{f(.)}, \eqn{g(.)}, \eqn{h(.)} etc (see 'Details').
#' @return The required inverting function, with single argument \code{x}. Its
#' \code{"varname"} attribute contains \code{varname} as a character string.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{plot.sitar}}
#' @examples
#' ## define varname variable
#' age <- 1:9
#'
#' ## simple case - age transformed to log(age)
#' (expr <- quote(log(age)))
#' ## transformed age
#' eval(expr)
#' ## inverting function, with "varname" attribute set to "age"
#' ifun(expr)
#' ## inverted transformed age identical to age
#' all.equal(age, ifun(expr)(eval(expr)))
#'
#' ## more complex case - age transformed to log age since conception
#' expr <- quote(log(age + 0.75))
#' ## inverting function
#' ifun(expr)
#' ## inverted transformed age identical to age
#' all.equal(age, ifun(expr)(eval(expr)))
#'
#' ## ludicrously complex case involving exp, log10, ^, pi and trigonometry
#' (expr <- quote((exp(sin(pi * log10(age + 0.75)/2) - 1)^4)))
#' ## inverting function
#' ifun(expr, verbose=TRUE)
#' ## identical to original
#' all.equal(age, ifun(expr)(eval(expr)))
#'
#' ## example of plot.sitar back-transforming transformed x and y in sitar models
#' ## fit sitar models
#' m1 <- sitar(x=age, y=height^2, id=id, data=heights, df=6)
#' m2 <- update(m1, x=log(age+0.75), y=height)
#'
#' ## default plot options for xfun & yfun back-transform x & y to original scales
#' ## xfun=ifun(x$call.sitar$x)
#' ## yfun=ifun(x$call.sitar$y)
#' ## compare mean curves for the three models with x & y on the original scales
#' plot(m1, 'd', las=1)
#' lines(m2, 'd', col=2)
#' @export ifun
ifun <- function(expr, verbose=FALSE) {

# returns number of names in expression ignoring pi
  vars <- function(expr) {
    (av <- all.vars(expr, unique=FALSE))[grep('^pi$', av, invert=TRUE)]
  }

# returns inverse function
  recur <- function(fun, funinv=quote(x), verbose=verbose) {
    fun <- as.expression(fun)[[1]]
#   if verbose print intermediate functions
    if (verbose) {
      print(fun)
      print(funinv)
      cat('---\n')
    }
#   if bracketed drop brackets
    while (length(fun) > 1 && fun[[1]] == as.name('('))
      fun <- fun[[2]]
    if (!is.name(fun)) {
#     element of expression containing x (either 2 or 3, ignoring leading symbol)
      x1 <- which(vapply(fun, function(f) length(vars(f)) == 1, TRUE))[-1]
#     if [cospi, sinpi, tanpi] drop 'pi' and multiply x by pi
      if (grepl('pi', fname <- as.name(fun[[1]]))) {
        fun[[1]] <- as.name(sub('pi', '', fname))
        f <- quote(x * pi)
        f[[2]] <- fun[[2]]
        fun[[2]] <- f
      }
#     expressions matching leading symbol with same length and x position
      nf <- which(vapply(fns, function(f)
        f[[1]] == fun[[1]] && length(f) == length(fun) && f[[x1]] == 'x', TRUE))
      if (length(nf) == 0)
        stop (paste('unrecognised name:', deparse(fun[[1]])))
#     if multiple matches check if numeric args are equal
      if (length(nf) > 1 && length(fun) == 3) {
#     compare elements not containing x
        nft <- which(vapply(fns[nf], function(f) f[[5 - x1]] == fun[[5 - x1]], TRUE))
        if (length(nft)) nf <- nf[nft]
      }
#     if more than one match use the first
      nf <- nf[[1]]
#     use complement of pair as inverse function
      fn2 <- fns[[nf - 1 + 2 * (nf %% 2)]]
#     element containing x in inverse function
      x2 <- which(as.list(fn2) == 'x')
#     if length 3 copy n
      if (length(fn2) == 3) {
#     function returns value for n
        f <- function(n) {}
        body(f) <- fn2[[5 - x2]]
        fn2[[5 - x2]] <- f(eval(fun[[5 - x1]]))
      }
#     update function and inverse function and repeat as necessary
      fun <- fun[[x1]]
      fn2[[x2]] <- funinv
      funinv <- fn2
      if (!is.name(fun)) {
        results <- recur(fun, funinv, verbose=verbose)
        fun <- results$fun
        funinv <- results$funinv
      }
    }
    return(list(funinv=funinv, fun=fun))
  }

# inverse function pairs
  fns <- quote(c(
    x+n, x-n,
    x*n, x/n,
    x^n, x^(1/n),
    sqrt(x), x^2,
    exp(x), log(x),
    expm1(x), log1p(x),
    n^x, log(x,n),
    log10(x), 10^x,
    log2(x), 2^x,
    n+x, x-n,
    n-x, n-x,
    n*x, x/n,
    n/x, n/x,
    +x, +x,
    -x, -x,
    I(x), I(x),
    cos(x), acos(x),
    sin(x), asin(x),
    tan(x), atan(x),
    cosh(x), acosh(x),
    sinh(x), asinh(x),
    tanh(x), atanh(x)
    ))
  fns[[1]] <- NULL

  varname <- vars(expr)
  if (length(varname) != 1)
    stop('expression should contain just one instance of one name')
  fn <- function(x) {}
  body(fn) <- with(fns, recur(expr, verbose=verbose))$funinv
  attr(fn, 'varname') <- varname
  fn
}
