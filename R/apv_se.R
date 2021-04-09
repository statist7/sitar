#' Bootstrap standard errors for SITAR peak velocity and age at peak velocity
#'
#' \code{apv_se} bootstraps a SITAR model to generate standard errors for
#' age at peak velocity (apv) and peak velocity (pv).
#'
#' If \code{plot} is TRUE, the original velocity curve is plotted along with
#' each bootstrap sample's pv versus apv.
#'
#' @param object SITAR model.
#' @param fun function to extract apv and pv from velocity curve (default getPeak),
#' alternative getTakeoff or getTrough.
#' @param nboot number of bootstrap samples (default 10).
#' @param seed integer to initialize the random number generator (default NULL).
#' @param plot logical to control plotting (default FALSE).
#' @param \dots optional arguments defining the velocity curve to be bootstrapped
#' (\code{plot.sitar} arguments \code{xfun}, \code{yfun}, \code{subset}, \code{ns}
#' or \code{abc}), and graphical \code{par} parameters.
#' @return a 2x2 array giving the mean and standard error of apv and pv, with
#' attribute "bs" a tibble containing the bootstrap estimates of apv and pv,
#' with NAs removed.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' ## fit sitar model for height
#' model <- sitar(x = age, y = height, id = id, data = heights, df = 4)
#'
#' ## bootstrap standard errors for age at peak velocity and peak velocity
#' output <- apv_se(model, nboot=3, seed=111, plot=TRUE)
#' @importFrom dplyr %>% everything
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr drop_na nest unnest
#' @importFrom rsample analysis bootstraps
#' @importFrom purrr map map_dfr
#' @importFrom rlang .data
#' @export
apv_se <- function(object,
                   fun = getPeak,
                   nboot = 10,
                   seed = NULL,
                   plot = FALSE,
                   ...) {
  # get data frame and object id
  df <- getData(object)
  id <- as.character(object$call.sitar$id)

  # get args
  dots <- eval(substitute(alist(...)))
  edots <- lapply(dots, eval, df)

  # function to process bootstrap samples
  fit_sitar_on_bootstrap <- function(split) {
    .x <- analysis(split) %>%
      unnest(cols = c(.data$data))
    eval(parse(text = ".x <<- .x"))
    obj <- try(update(object, data = .x, start = fixef(object)))
    if (any(class(obj) %in% 'sitar')) {
      vel <- do.call('plot_v', c(list(x = obj), edots))
      fun(vel)
    } else
      c(x = NA_real_, y = NA_real_)
  }

  set.seed(seed)

  df <- df %>%
    # generate bootstrap splits
    nest(data = -id) %>%
    bootstraps(times = nboot) %>%
    # generate bootstrap samples
    mutate(model = map(.data$splits, fit_sitar_on_bootstrap))

  # extract apv and pv
  apv <- map_dfr(df$model, ~ as_tibble(t(.x))) %>%
    drop_na

  vel <- do.call('plot_v', c(list(x = object), edots))
  peak <- fun(vel)
  names(apv) <- names(vel) <- names(peak) <- c('apv', 'pv')

  if (plot) {
    localplot <- function(..., subset, abc, xfun, yfun, ns)
      plot(...)
    xy <- rbind(vel, apv)
    ARG <-
      setNames(lapply(xy, range, na.rm = TRUE), c('xlim', 'ylim')) # default plot
    do.call('localplot', c(list(x = apv), ARG, edots)) # plot apv and pv
    lines(vel, lwd = 2) # add velocity curve
    abline(v = peak[1], h = peak[2], lty = 3)
  }

  se <- vapply(apv, sd, na.rm = TRUE, 1)
  output <- rbind(peak, se)
  attr(output, 'bs') <- apv
  return(output)
}
