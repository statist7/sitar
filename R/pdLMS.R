#' Plot frequency distributions(s) for given L, M and S values in LMS method
#'
#' The LMS method defines frequency distributions in terms of L, M and S parameters.
#' \code{pdLMS} plots one or more LMS distributions and optionally returns specified
#' centiles on each distribution.
#'
#' L, M and S should all be the same length, recycled if necessary.
#'
#' @param L vector of Box-Cox transformation (lambda) values, L in the LMS
#' method (default 1 corresponding to the Normal distribution).
#' @param M vector of medians (mu), M in the LMS method (default 1).
#' @param S vector of coefficients of variation (sigma), S in the LMS method
#' (default 0.2).
#' @param zcent optional vector of z-scores for conversion to the measurement
#' scale under each distribution.
#' @param zlim scalar defining z-score limits underlying x-axis (default 3.5).
#' @param N number of points per distribution curve (default 1000).
#' @param plot logical for plotting (default TRUE).
#' @param \dots Further graphical parameters (see \code{\link{par}}) may also
#' be supplied as arguments, particularly colour \code{col}, line type \code{lty},
#' line width \code{lwd} and character \code{pch}.
#' @return An invisible list with the following components:
#' \item{x}{vector of x values for plotting.}
#' \item{density}{matrix of densities for each distribution.}
#' \item{centile}{matrix of measurement centiles corresponding to \code{zcent}
#' under each distribution.}
#' The distributions can be plotted with \code{matplot(x, density, type='l')}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{z2cent}}, \code{\link{LMS2z}}, \code{\link{cLMS}}
#' @keywords arith
#' @examples
#'
#' ## plot normal distribution
#' pdLMS()
#' ## compare variety of distributions
#' ## with centiles corresponding to +3 z-scores
#' pdLMS(L=-2:3, M=2:3, S=1:3/10, zcent=3, lty=1)
#'
#' @importFrom graphics matplot matpoints
#' @export pdLMS
pdLMS <- function(L = 1, M = 1, S = 0.2, zcent = NULL, zlim = 3.5,
                  N = 1000, plot = TRUE, ...) {
  L[L == 0] <- 1e-7
  LMS <- data.frame(L, M, S)
  LSz <- L * S * abs(zlim)
  xr <- 0:1
  if (min(1 - LSz) > 0)
    xr[1] <- min(M * (1 - LSz) ^ (1 / L))
  x <- cLMS(abs(zlim), L, M, S)
  xr[2] <- if (length(which.max(x)) > 0)
    x[which.max(x)]
  else
    M * exp(S * abs(zlim))
  delta <- (xr[2] - xr[1]) / N
  x <- xr[1] + delta * (1:N - 0.5)
  density <- vapply(seq(nrow(LMS)), function(i) {
    with (LMS[i,], {
      z <- zLMS(x, L, M, S)
      z <- dnorm(z) * x ^ (L - 1)
      z / sum(z) / delta
    })
  }, x)
  if (!is.null(zcent)) {
    centile <- t(cLMS(as.matrix(zcent), L, M, S))
    dimnames(centile) <- list(zcent, seq(nrow(LMS)))
  }
  else
    centile <- NULL
  if (plot) {
    matplot(x, density, type = "l", ...)
    abline(h = 0, col = 8)
    if (!is.null(zcent))
      matpoints(centile, rep(0, nrow(centile)), ...)
  }
  invisible(list(x = x, density = density, centile = centile))
}
