#' Simulate datasets for SITAR models
#'
#' \code{simulate.sitar} Takes a SITAR model and generates datasets that on
#' average match its mean curve, while retaining the subject data structure.
#'
#' The values of \code{x} and \code{y} in \code{object} are adjusted to match
#' the model mean curve, with \code{y} the fitted values. Then new subject
#' random effects are drawn from the \code{object} covariance matrix, and used
#' to convert subject data from the mean curve to their own individual curves.
#' Residual error is then added back to the \code{y} values.
#'
#' @param object a SITAR model.
#' @param nsim number of datasets to simulate (default 1).
#' @param seed integer to initialize the random number generator (default NULL).
#' @param \dots optional additional arguments. None are used.
#' @return A list of \code{nsim} data frames (or if \code{nsim = 1} a data frame),
#' returned as tibbles, with the same variable names and dimensions as the data
#' used for \code{object}, and the \code{x} and \code{y} values updated.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#'
#' data(heights)
#' ## fit sitar model for height
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=4)
#'
#' ## plot unadjusted data as growth curves
#' plot(m1, opt='u')
#'
#' ## draw fitted height distance curves coloured by subject, using ggplot
#' \dontrun{
#' require(ggplot2)
#' ggplot(plot_D(m1), aes(.x, .y, colour=.id)) +
#' labs(x='age', y='height') +
#' geom_line(show.legend=FALSE)
#' }
#'
#' @importFrom tibble tibble
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm
#' @export
simulate.sitar <- function(object, nsim=1, seed=NULL, ...) {
  set.seed(seed)
  data <- tibble::as_tibble(getData(object))
  nxy <- unlist(lapply(object$call.sitar[2:3], deparse))
  nx <- which(names(data) %in% nxy[['x']])
  ny <- which(names(data) %in% nxy[['y']])
  N <- object$dims$N
  n <- object$dims$ngrps[1]
  gvc <- getVarCov(object)
  mu <- rep(0, dim(gvc)[1])
  rsd <- object$sigma
  id <- getGroups(object)
  old <- xyadj(object, y=fitted(object))

  . <- setNames(lapply(seq.int(nsim), function (i) {
    ref <- MASS::mvrnorm(N, mu, gvc)
    ref <- scale(ref, TRUE, FALSE)
    new <- xyadj(object, x=old$x, y=old$y, tomean=FALSE,
                  abc=ref[id, , drop=FALSE])
    data[, nx] <- new$x
    data[, ny] <- new$y + rnorm(N, sd=rsd)
    data
  }), seq.int(nsim))
  if (nsim == 1)
    .[[1]]
  else
    invisible(.)
}
