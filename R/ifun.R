#' Invert an expression defining a data transformation
#'
#' Enables a transformed variable to be back-transformed to its original scale,
#' e.g. for plotting purposes, by creating a function to invert the
#' transforming expression.
#'
#' \code{ifun} handles all the invertible functions in the 'Math' and 'Ops'
#' groups.
#'
#' To explain how \code{ifun} works, consider as an example variants of the
#' \code{sitar} model \code{height ~ age} where \code{age} and/or \code{height}
#' may be transformed, e.g. \code{height ~ log(age)} or \code{log(height) ~
#' sqrt(age)}. Thus each model is of the form \code{y ~ x} but the units of
#' \code{x} and \code{y} vary.
#'
#' It is useful to be able to compare the models by plotting the fitted curves
#' in the original units. This is done by applying suitable functions to invert
#' the transformed variables prior to plotting. For example the
#' \code{function(x) exp(x)} back-transforms \code{log(age)} and
#' \code{log(height)}. \code{ifun} automates this process by creating a
#' function to back-transform the general expression \code{fun(x)} to \code{x}.
#'
#' Structuring \code{fun} suitably ensures it can be inverted. It should
#' contain a single mention of a single variable (named \code{varname} here),
#' and possibly functions such as \eqn{f(.)}, \eqn{g(.)}, \eqn{h(.)} etc,
#' each of which involves (a function of) (a single mention of) \code{varname}
#' and up to one numerical expression, such that \code{fun} =
#' \eqn{f(g(h((varname))))}. The number of such functions is in principle
#' unlimited.
#'
#' \code{ifun} then returns
#' \eqn{h^{-1}(g^{-1}(f^{-1}((eval(fun)))))}{h^-1(g^-1(f^-1((eval(fun)))))}
#' as a function,
#' with any numerical expressions evaluated, which means that \code{fun} is
#' invertible so long as the individual functions are invertible. Note though
#' that the function being invertible does not guarantee that variables can
#' always be back-transformed, as for example the function may introduce
#' \code{NaN}s.
#'
#' Note that \code{\link{plot.sitar}} uses \code{ifun} as the default for
#' arguments \code{xfun} and \code{yfun}, where the corresponding values of
#' \code{fun} are extracted from the model's \code{sitar} call.
#'
#' @param fun a language object defining the expression to be inverted, best
#' \code{quote}d to avoid evaluation.
#' @return A list of length two: \item{fn}{the inverse function with argument
#' \code{x} which applied to \code{eval(fun)} returns \code{varname}.}
#' \item{varname}{the name of the variable in \code{fun} (given as a character
#' string).}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{plot.sitar}}
#' @examples
#'
#' ## define age
#' age <- 1:9
#'
#' ## simple case - transform age to log(age)
#' (fun <- quote(log(age)))
#'
#' (transformed.age <- eval(fun))
#' (inverting.function <- ifun(fun)$fn)
#' (inverted.transformed.age <- inverting.function(transformed.age))
#'
#' ## is inverted transformed age identical to age?
#' all.equal(age, inverted.transformed.age)
#'
#'
#' ## more complex case - transform age to log age since conception
#' fun <- quote(log(age + 0.75))
#'
#' (transformed.age <- eval(fun))
#' (inverting.function <- ifun(fun)$fn)
#' (inverted.transformed.age <- inverting.function(transformed.age))
#'
#' ## identical to original?
#' all.equal(age, inverted.transformed.age)
#'
#'
#' ## ludicrously complex case involving exp, log10, ^, pi and trigonometry
#' (fun <- quote((exp(sin(pi * log10(age + 0.75)/2) - 1)^4)))
#'
#' (transformed.age <- eval(fun))
#' (inverting.function <- ifun(fun)$fn)
#' (inverted.transformed.age <- inverting.function(transformed.age))
#'
#' ## identical to original?
#' all.equal(age, inverted.transformed.age)
#'
#'
#' ## example of plot.sitar back-transforming transformed x and y in sitar models
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=6)
#' m2 <- update(m1, y=height^2)
#' m3 <- update(m1, x=log(age+0.75))
#'
#' ## default plot settings back-transform x and y to original scales
#' plot(m1, 'd')
#' lines(m2, 'd', col=2)
#' lines(m3, 'd', col=3)
#'
#' @export ifun
ifun <- function(fun) {

# returns number of names in function ignoring pi
  nvars <- function(fun) {
    length((av <- all.vars(fun, unique=FALSE))[grep('^pi$', av, invert=TRUE)])
  }

  if (nvars(fun) != 1)
    stop('expression should contain just one instance of one name')
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

# returns inverse function
  recur <- function(fun, funinv=quote(x)) {
    fun <- as.expression(fun)[[1]]
#   if bracketed drop brackets
    while (length(fun) > 1 && fun[[1]] == as.name('('))
      fun <- fun[[2]]
    if ((ne <- length(fun)) > 1) {
#     element of expression containing varname (either 2 or 3, ignoring leading symbol)
      x1 <- which(vapply(fun, function(f) nvars(f) == 1, TRUE))[-1]
#     element not containing varname (only relevant when ne == 3)
      n1 <- 5 - x1
#     if [cospi, sinpi, tanpi] drop 'pi' and multiply x by pi
      if (grepl('pi', fname <- as.name(fun[[1]]))) {
        fun[[1]] <- as.name(sub('pi', '', fname))
        f <- quote(x * pi)
        f[[2]] <- fun[[2]]
        fun[[2]] <- f
      }
#     expressions matching leading symbol with same length and x arg position
      nf <- which(vapply(fns, function(f)
        f[[1]] == fun[[1]] && length(f) == length(fun) && f[[x1]] == 'x', TRUE))
      if (length(nf) == 0)
        stop (paste('unrecognised name:', deparse(fun[[1]])))
#     if multiple matches check if numeric args are equal
      if (length(nf) > 1 && ne == 3) {
        nft <- which(vapply(fns[nf], function(f) f[[n1]] == fun[[n1]], TRUE))
        if (length(nft)) nf <- nf[nft]
      }
#     if more than one match use the first
      nf <- nf[[1]]
#     use complement of pair as inverse function
      fn2 <- fns[[nf - 1 + 2 * (nf %% 2)]]
#     identify position of x arg in inverse function
      x2 <- which(as.list(fn2) == 'x')
#     if length 3 copy n arg
      if (length(fn2) == 3) {
        n2 <- 5 - x2
  #       function returns value for n arg
        f <- function(n) {}
        body(f) <- fn2[[n2]]
        fn2[[n2]] <- f(eval(fun[[n1]]))
      }
#     copy x from current inverse function
      fn2[[x2]] <- funinv
#     update function and inverse function and repeat as necessary
      fun <- fun[[x1]]
      if (is.name(fun))
        funinv <- fn2
      else {
        funinv <- (results <- recur(fun, fn2))$fn
        fun <- results$varname
      }
    }
    return(list(fn=funinv, varname=fun))
  }

  results <- with(fns, recur(fun))
  fn <- function(x) {}
  body(fn) <- results$fn
  return(list(fn=fn, varname=deparse(results$varname)))
}
