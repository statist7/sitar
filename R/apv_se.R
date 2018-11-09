#' Bootstrap standard errors for SITAR peak velocity and age at peak velocity
#'
#' \code{apv_se} bootstraps a SITAR model to generate standard errors for
#' age at peak velocity (apv) and peak velocity (pv).
#'
#' The mean velocity curve to be bootstrapped can be modified with
#' arguments \code{subset}, \code{abc}, \code{xfun}, \code{yfun} or \code{ns}.
#'
#' If \code{plot} is TRUE, the original velocity curve is plotted along with
#' each bootstrap sample's pv versus apv.
#'
#' @param object SITAR model.
#' @param nboot number of bootstrap samples (default 10).
#' @param seed integer to initialize the random number generator (default NULL).
#' @param plot logical to control plotting (default FALSE).
#' @param \dots optional arguments defining the velocity
#' curve to be bootstrapped, and the plot. See Details.
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

  dots <- eval(substitute(alist(...)))
  edots <- lapply(dots, eval, getData(object))
  vel <- do.call('plot', c(list(x=object, opt='v', returndata=TRUE), edots))
  peak <- getPeakTrough(vel)

  run_sitar <- quote(update(object, data=.df, id=.newid, start=fixef(object)))

  apv <- na.omit(as_tibble(t(vapply(df, function(z) { # fit sitar and extract apv and pv
    eval(parse(text=".df <<- z"))
    obj <- try(eval(run_sitar))
    if (any(class(obj) %in% 'sitar')) {
      edots <- lapply(dots, eval, getData(obj))
      vel <- do.call('plot', c(list(x=obj, opt='v', returndata=TRUE), edots))
      getPeakTrough(vel)
    } else
      c(NA, NA)
  }, vector('numeric', 2)))))

  names(apv) <- names(vel) <- names(peak) <- c('apv', 'pv')

  if (plot) {
    localplot <- function(..., subset, abc, xfun, yfun, ns) plot(...)
    xy <- rbind(vel, apv)
    ARG <- setNames(lapply(xy, range), c('xlim', 'ylim')) # default plot
    do.call('localplot', c(list(x=apv), ARG, edots)) # plot apv and pv
    lines(vel, lwd=2) # add velocity curve
    abline(v=peak[1], h=peak[2], lty=3)
  }

  se <- vapply(apv, sd, 1.0)
  output <- rbind(peak, se)
  attr(output, 'bs') <- apv
  return(output)
}
