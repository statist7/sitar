#' Extract individual estimates of peak velocity and age at peak velocity
#'
#' \code{getDV} takes a point on the mean fitted velocity curve, defined
#' typically as the age at peak velocity or an equivalent landmark age, and maps
#' it onto the corresponding point on individual velocity curves taking into
#' account their timing and intensity.
#'
#' SITAR is a shape-invariant model, so if there is a turning point on the mean
#' velocity curve, e.g. a peak in puberty, then there will be a corresponding
#' turning point on the velocity curves of all individuals, at an age depending
#' on their timing and intensity. This applies even to individuals who lack
#' measurements at that age, have stopped earlier or started later, so their
#' growth curve is incomplete and lacks the turning point. The returned variable
#' \code{missing} flags such individuals.
#'
#' Note that 'D' and 'V' in \code{getDV} correspond to the plot options for
#' individual Distance and Velocity curves.
#'
#' @param object a SITAR model.
#' @param x the age of interest, specified either as 'apv', the mean age at peak
#'   velocity, or 'ato', the mean age at take-off, or a numerical age.
#' @return A tibble with one row per individual and five columns, where id, age
#'   and distance have names corresponding to those used in the SITAR call:
#'   \item{id}{subject id.} \item{age}{subject's age corresponding to \code{x}.}
#'   \item{distance}{subject's distance at specified age.}
#'   \item{velocity}{subject's velocity at specified age.}
#'   \item{missing}{logical where TRUE means subject's specified age lies outside
#'   their measurement range.}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#'
#' data(heights)
#' library(dplyr)
#' library(ggplot2)
#' theme_set(theme_bw())
#'
#' # fit sitar model
#' model <- sitar(x = log(age), y = height, id = id, data = heights, df = 5)
#' (dv <- getDV(model))
#' (id_missing <- dv %>%
#'   filter(missing) %>%
#'   pull(id))
#'
#' # plot individual velocity curves
#' ggplot(plot_V(model), aes(age, height, group = id, colour = id)) +
#'     theme(legend.position = 'inside',
#'           legend.position.inside = c(0.9, 0.6)) +
#'     geom_line() +
#'
#' # add individual peak velocities
#' geom_point(aes(y = velocity), data = dv) +
#'
#' # highlight subjects 2, 7 and 12 with dashed lines
#' # despite incomplete curves their peak velocities are estimated
#' geom_line(data = . %>% filter(id %in% id_missing), linetype = 2, colour = 'white') +
#'             geom_point(aes(y = velocity), data = dv %>%
#'                          filter(id %in% id_missing), shape = 1, size = 3)
#'
#' @importFrom dplyr all_of between inner_join mutate nest_by pick pull
#'   rename_with select summarise
#' @importFrom tibble tibble
#' @export
getDV <- function(object, x = 'apv') {
  stopifnot('object must be a SITAR model' = inherits(object, 'sitar'),
            'x should be either "apv" or "ato" or a number' = x %in% c('apv', 'ato') | is.numeric(x),
            'x should not be missing' = !is.na(x))
  id <- rownames(ranef(object))
  idname <- all.vars(object$call.sitar$id)
  yname <- all.vars(object$call.sitar$y)
  xfun <- object$call.sitar$x
  xname <- all.vars(xfun)
  # mean age at peak/trough
  if (x == 'apv') {
    assign(xname, getPeak(plot_v(object, ns = 201))[['x']])
  } else if (x == 'ato') {
    assign(xname, getTrough(plot_v(object, ns = 201))[['x']])
  } else {
    assign(xname, x)
  }
  stopifnot('no turning point' = !is.na(get(xname)))
  # subject age at peak/trough
  df <- xyadj(object, x = eval(xfun), id = id, tomean = FALSE)[['x']] %>%
    ifun(xfun)(.) %>%
    tibble(id = id, x = .) %>%
    rename_with(~c(idname, xname))
  # subject size and velocity at peak/trough
  df <- inner_join(df,
                   predict(object, df, deriv = 0:1) %>%
                     rename_with(~c(idname, yname, 'velocity')),
                   by = idname)
  inner_join(df,
             getData(object) %>%
               select(all_of(c(idname, xname))) %>%
               nest_by(pick(1)) %>%
               summarise(xmin = min(pick(1)),
                         xmax = max(pick(1)),
                         .groups = 'drop'),
             by = idname) %>%
    mutate(missing = !between(pick(2) %>% pull(), .data$xmin, .data$xmax)) %>%
    select(-c(.data$xmin, .data$xmax))
}

