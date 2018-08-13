#' Simulate datasets for SITAR models
#'
#' \code{simulate.sitar} takes a SITAR model and simulates new datasets for
#' the model that on average produce the same mean curve, while retaining
#' the subject data structure.
#'
#' Each simulation generates a \code{tibble} with new values of \code{x} and
#' \code{y}, and other variables (e.g. \code{id}) unchanged. The values of
#' \code{y} are first set to the level 0 fitted values, and those for \code{x}
#' are adjusted correspondingly using the subject random effects. Thus
#' \code{plot(y ~ x)} gives the mean curve (option 'd' in \code{plot}). New
#' subject random effects are then drawn randomly from the \code{object}
#' covariance matrix, centred if the corresponding fixed effects are in
#' \code{object}, and used to convert \code{x} and \code{y} from the mean curve
#' to the individual subject curves (option 'D' in \code{plot}). Residual error
#' (SD = \code{object$sigma}) is then added to \code{y}.
#'
#' @param object SITAR model.
#' @param nsim number of datasets to simulate (default 1).
#' @param seed integer to initialize the random number generator (default NULL).
#' @param \dots optional additional arguments. None are used.
#' @return A data frame (or if \code{nsim > 1} a list of data frames), each
#' returned as a \code{tibble}, identical to the original data frame used
#' for \code{object} except with updated values for \code{x} and \code{y}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' ## fit sitar model for height
#' model <- sitar(x = age, y = height, id = id, data = heights, df = 4)
#'
#' ## simulate sitar datasets
#' sim <- simulate(model, 3)
#'
#' ## fit sitar models to data
#' models <- lapply(sim, function(zz) try(sitar(age, height, id, zz, 4)))
#'
#' ## delete any dud models
#' dud <- vapply(models, function(m) 'try-error' %in% class(m), FALSE)
#' sim[dud] <- models[dud] <- NULL
#'
#' # calculate mean/sd of age at peak velocity
#' . <- vapply(seq.int(sim), function(i) {
#'   zz <<- sim[[i]]
#'   plot(models[[i]], apv=TRUE)$apv
#' }, c(0, 0))
#' rm(zz)
#' plot(model, apv=TRUE)
#' apply(., 1, mean, na.rm=TRUE)
#' apply(., 1, sd, na.rm=TRUE)
#' @importFrom tibble tibble
#' @importFrom MASS mvrnorm
#' @importFrom stats simulate rnorm
#' @export
simulate.sitar <- function(object, nsim=1, seed=NULL, ...) {
  set.seed(seed)
  data <- tibble::as_tibble(getData(object))
  nxy <- unlist(lapply(object$call.sitar[2:3], deparse))
  nx <- which(names(data) %in% nxy[['x']])
  ny <- which(names(data) %in% nxy[['y']])
  N <- object$dims$N
  n <- object$dims$ngrps[1]
  center <- letters[1:3] %in% names(fixef(object))
  gvc <- getVarCov(object)
  mu <- rep(0, dim(gvc)[1])
  rsd <- object$sigma
  id <- getGroups(object)
  old <- xyadj(object, y=fitted(object))

  . <- setNames(lapply(seq.int(nsim), function (i) {
    ref <- MASS::mvrnorm(n, mu, gvc)
    ref[, center] <- scale(ref[, center], TRUE, FALSE)
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
