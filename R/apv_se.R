#' Bootstrap standard errors for SITAR peak velocity and age at peak velocity
#'
#' \code{apv_se} bootstraps a SITAR model to generate standard errors for
#' age at peak velocity (apv) and peak velocity (pv).
#'
#' If \code{plot} is TRUE, the mean velocity curve is plotted along with
#' each bootstrap sample's pv versus apv.
#'
#' @param object SITAR model.
#' @param nboot number of bootstrap samples (default 10).
#' @param seed integer to initialize the random number generator (default NULL).
#' @param plot logical to control plotting (default FALSE).
#' @param \dots optional plot arguments.
#' @return a 2x2 array giving the mean and se of apv and pv, with attribute "bs"
#' a tibble containing the bootstrap estimates of apv and pv, with NAs removed.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' ## fit sitar model for height
#' model <- sitar(x = age, y = height, id = id, data = heights, df = 4)
#'
#' ## bootstrap standard errors for age at peak velocity and peak velocity
#' output <- apv_se(model, nboot=3, seed=111, plot=TRUE)
#' @importFrom dplyr %>%
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr nest unnest
#' @importFrom rsample bootstraps
#' @importFrom purrr map
#' @importFrom rlang enquo
#' @export
apv_se <- function(object, nboot=10, seed=NULL, plot=FALSE, ...) {
  id <- object$call.sitar$id
  id <- enquo(id)

  set.seed(seed)
  bs <- getData(object) %>% # generate bootstrap splits
    nest(-!!id) %>%
    bootstraps(times=nboot)
  df <- map(bs$splits, ~as_tibble(.) %>% # generate bootstrap samples
              unnest(.id='.newid'))

  run_sitar <- quote(update(object, data=.df, id=.newid, start=fixef(object)))
  apv <- na.omit(as_tibble(t(vapply(df, function(z) { # fit sitar and extract apv and pv
    eval(parse(text=".df <<- z"))
    obj <- try(eval(run_sitar))
    if (any(class(obj) %in% 'sitar'))
      getPeakTrough(plot_v(obj))
    else
      c(NA, NA)
  }, vector('numeric', 2)))))

  vel <- plot_v(object)
  peak <- getPeakTrough(vel)
  names(apv) <- names(vel) <- names(peak) <- c('apv', 'pv')

  if (plot) {
    xy <- rbind(vel, apv)
    ARG <- setNames(lapply(xy, range), c('xlim', 'ylim')) # default plot
    dots <- list(...)
    ARG[names(dots)] <- dots
    do.call('plot', c(list(x=apv), ARG)) # plot apv and pv
    lines(vel, lwd=2) # add velocity curve
    abline(v=peak[1], h=peak[2], lty=3)
  }

  se <- vapply(apv, sd, 1)
  output <- rbind(peak, se)
  attr(output, 'bs') <- apv
  return(output)
}
