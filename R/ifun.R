ifun <- function(fun) {
#	returns inverse of function fun
###########################
  recur <- function(fun, funinv=quote(x)) {
    fun <- as.expression(fun)[[1]]
#   if bracketed drop brackets
    if (length(fun) > 1 && fun[[1]] == as.name('(')) fun <- fun[[2]]
    ne <- length(fun)
    if (ne > 1) {
#     expression element containing language (ignoring leading symbol)
      x1 <- which(vapply(fun, is.language, TRUE))[[2]]
      if (length(x1) != 1) stop('there should be just one name in the expression')
#     element containing numeric (may be absent)
      n1 <- which(vapply(fun, is.numeric, TRUE))
#     expressions of same length matching leading symbol and x arg position
      nf <- which(vapply(fns, function(f) length(f) == length(fun) && f[[1]] == fun[[1]] && f[[x1]] == quote(x), TRUE))
      if (length(nf) == 0) stop (paste('unrecognised symbol:', deparse(fun[[1]])))
      if (length(nf) > 1) {
#       multiple matches - check if numeric args are equal
        nft <- which(vapply(fns[nf], function(f) length(n1) > 0 && f[[n1]] == fun[[n1]], TRUE))
        if (length(nft)) nf <- nf[nft]
      }
#     if more than one match use the first
      nf <- nf[[1]]
#     use complement of pair as inverse function
      fn2 <- if (nf %% 2) fns[[nf + 1]] else fns[[nf - 1]]
#     identify position of x arg in inverse function
      x2 <- which(as.list(fn2) == quote(x))
#     if length 3 then position of n arg is complement of 5
      if (length(fn2) == 3) {
        n2 <- 5 - x2
#       function returns value for n
        f <- function(n) {}
        body(f) <- fn2[[n2]]
        fn2[[n2]] <- f(fun[[n1]])
      }
#     copy x from current inverse function
      fn2[[x2]] <- funinv
#     update function and inverse function and repeat as necessary
      fun <- fun[[x1]]
      if (is.symbol(fun)) funinv <- fn2 else {
        funinv <- (results <- recur(fun, fn2))$fn
        fun <- results$varname
      }
    }
    return(list(fn=funinv, varname=fun))
  }
###########################
  fns <- quote(c(x+n, x-n, x*n, x/n, x^n, x^(1/n), sqrt(x), x^2, exp(x), log(x),
                 expm1(x), log1p(x), n^x, log(x,n), log10(x), 10^x, log2(x),
                 2^x, n+x, x-n, n-x, n-x, n*x, x/n, n/x, n/x, +x, +x, -x, -x,
                 cos(x), acos(x), sin(x), asin(x), tan(x), atan(x),
                 cosh(x), acosh(x), sinh(x), asinh(x), tanh(x), atanh(x)))
  fns[[1]] <- NULL
  fn <- function(x) {}
  results <- with(fns, recur(fun))
  body(fn) <- results$fn
  return(list(fn=fn, varname=deparse(results$varname)))
}
