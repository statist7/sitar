#' Convert between IOTF, WHO and CDC prevalence rates for child thinness,
#' overweight and obesity
#'
#' Child thinness, overweight and obesity are defined as the child's body mass
#' index (BMI) lying beyond a pre-specified reference cutoff. Three references
#' are compared here: IOTF (International Obesity Task Force), WHO (World Health
#' Organization) and CDC (US Centers for Disease Control and Prevention), each
#' of which have their own cutoffs. \code{ob_convertr} takes age-sex-specific
#' prevalence rates of thinness, overweight and obesity based on one cutoff, and
#' converts them to rates based on a different cutoff, using a novel estimation
#' algorithm.
#'
#' The IOTF cutoffs correspond to the value of BMI (kg/m^2) at age 18: IOTF 35
#' (morbid obesity), IOTF 30 (obesity), IOTF 25 (overweight), IOTF 18.5 (grade 1
#' thinness), IOTF 17 (grade 2 thinness) and IOTF 16 (grade 3 thinness).
#'
#' The WHO cutoffs correspond to BMI z_scores. Age 5-19 years, WHO +2 (obesity),
#' WHO +1 (overweight) and WHO -2 (thinness). Age 0-5 years, WHO +3 (obesity),
#' WHO +2 (overweight) and WHO -2 (thinness).
#'
#' The CDC cutoffs correspond to BMI centiles: CDC 95 (obesity), CDC 85
#' (overweight) and CDC 5 (thinness).
#'
#' Note 1: the overweight category needs to be analysed as overweight plus
#' obesity. To predict overweight excluding obesity, first calculate predicted
#' overweight plus obesity then subtract predicted obesity.
#'
#' Note 2: the category labels are harmonised and not necessarily as originally
#' defined.
#'
#' The conversion algorithm exploits the fact that all three references are
#' based on the LMS method, which allows prevalence to be converted to a common
#' BMI centile and z-score scale.
#'
#' The algorithm is commutative, which means that converting a prevalence rate
#' from cutoff A to cutoff B and then from B to A returns the original value.
#'
#' @param prev vector of age-sex-specific percentage prevalence rates.
#' @param age vector of ages between 2 and 18 years corresponding to each rate.
#' @param sex vector of the sexes corresponding to each rate, coded as either
#'   'boys/girls' or 'male/female' or '1/2' (upper or lower case, and only the
#'   first character considered).
#' @param from the BMI cutoff (see Details) on which the prevalence is based.
#' @param to the BMI cutoff (see Details) on which to base the converted
#'   prevalence.
#' @param prev_true optional vector of known percentage prevalence rates
#'   corresponding to \code{to}, for validation purposes.
#' @param report character controlling the format of the returned data: 'vector'
#'   for the estimated prevalence rates, 'wider' for the working tibble in wide
#'   format, i.e. the \code{from} and \code{to} data side by side, or 'longer'
#'   for the tibble in long format, i.e. two rows per rate, one for \code{from}
#'   and one for \code{to}.
#' @param plot character controlling what if anything is plotted: 'no' for no
#'   plot, 'density' to display the BMI density distributions and cutoffs
#'   corresponding to \code{from} and \code{to}, or 'compare' to display the
#'   predicted prevalence rates plotted against the observed rates in
#'   \code{prev_true}.
#' @param data data frame containing \code{prev}, \code{age}, \code{sex} and
#'   \code{prev_true}.
#'
#' @return Either the converted prevalence rates or a plot visualizing the
#'   findings, depending on the \code{report} and \code{plot} settings.
#'
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' ## convert 10% IOTF overweight prevalence (cutoff IOTF 25) in 8-year-old boys
#' ## to the overweight prevalence based on WHO, i.e. cutoff WHO +1
#' ob_convertr(prev = 10, age = 8, sex = 'boys', from = 'IOTF 25', to = 'WHO +1')
#'
#' ## compare the BMI density functions and cutoffs for IOTF 25 and WHO +1
#' ob_convertr(prev = 10, age = 8, sex = 'boys', from = 'IOTF 25', to = 'WHO +1', plot = 'density')
#'
#' ## convert IOTF overweight prevalence to WHO overweight prevalence
#' ## and compare with true value - boys and girls age 7-17
#' ## note the need to first add obesity prevalence to overweight prevalence
#' data(deren)
#' deren <- within(deren, {
#'   IOTF25 = IOTF25 + IOTF30
#'   `WHO+1` = `WHO+1` + `WHO+2`})
#' ob_convertr(prev = IOTF25, age = Age, sex = Sex, from = 'IOTF 25', to = 'WHO +1',
#'    prev_true = `WHO+1`, data = deren, plot = 'compare')

#' @importFrom forcats fct_inorder fct_collapse fct_relabel
#' @importFrom ggplot2 ggplot xlab ylab geom_path geom_point geom_vline
#'   geom_abline scale_x_continuous scale_y_continuous aes
#' @importFrom tibble tibble column_to_rownames
#' @importFrom dplyr mutate rename filter transmute bind_rows bind_cols select
#'   across left_join pull contains starts_with ends_with if_else n rowwise
#'   group_by ungroup
#' @importFrom rlang .data enquo
#' @importFrom stats pnorm qnorm
#' @importFrom tidyr pivot_wider pivot_longer drop_na
#' @export ob_convertr
ob_convertr <- function(prev = 50, age, sex, from, to, prev_true = NA, report = c('vector', 'wider', 'longer'),
                     plot = c('no', 'density', 'compare'), data = parent.frame()) {

  # check sex contains only 1/M/B/TRUE or 2/F/G
  test_sex <- function(sex) {
    fsex <- toupper(substr(sex, 1, 1))
    fsex <- factor(fsex, levels = c(1:2, 'M', 'F', 'B', 'G', 'T'))
    levels(fsex) <- c(rep(1:2, 3), 1)
    fsex <- fct_collapse(fsex, boys = '1', girls = '2')
    droplevels(fsex)
  }

  # return sitar code for reference
  ref_sitar <- function(x) {
    x <- unique(x)
    stopifnot('cutoff not unique' = length(x) == 1L)
    x <- sub('^(.*) .*$', '\\1', x) # drop any cutoff
    f <- factor(x, levels = c('CDC', 'IOTF', 'WHO'))
    levels(f) <- c('cdc2000', 'iotf', 'who0607')
    as.character(f)
  }

  # create cutoffs matrix
  ref <- cutoff <- f <- LMS <- L <- M <- S <- NULL
  cutoffs <- tibble(ref = c('CDC', 'IOTF', 'WHO'),
                    cutoff = list(c(5, 85, '95'),
                                  c(16:17, 18.5, 25, 30, '35'),
                                  c(-2:-1, paste0('+', 1:3))),
                    f = list(function(cutoff, sex) qnorm(cutoff / 100),
                             function(cutoff, sex) LMS2z(18, cutoff, sex, 'bmi', 'iotf'),
                             function(cutoff, sex) cutoff),
                    sex = list(c('boys', 'girls'))) %>%
    unnest(cutoff) %>% unnest(sex) %>% rowwise %>%
    mutate(z = do.call(.data$f, list(as.numeric(.data$cutoff), sex)),
           cutoff = paste(.data$ref, .data$cutoff)) %>%
    select(-c(ref, f)) %>%
    pivot_wider(names_from = sex, values_from = 'z') %>%
    column_to_rownames('cutoff') %>%
    as.matrix

  # check from and to
  from <- unique(toupper(from))
  to <- unique(toupper(to))
  stopifnot('`from` should be length 1' = length(from) == 1L,
            '`to` should be length 1' = length(to) == 1L,
            '`to` is the same as `from`' = to != from,
            '`from` not recognised' = from %in% dimnames(cutoffs)[[1]],
            '`to` not recognised' = to %in% dimnames(cutoffs)[[1]])

  # create tibble of age sex and prev
  data <- if (identical(data, parent.frame())) {
    eval.parent(substitute(
      data.frame(age = age, sex = test_sex(sex), prev = prev, prev_true = prev_true)))
  } else {
    prev <- enquo(prev)
    age <- enquo(age)
    sex <- enquo(sex)
    prev_true <- enquo(prev_true)
    data %>%
      transmute(age = !!age, sex = test_sex(!!sex), prev = !!prev, prev_true = !!prev_true)
  }

  if (max(data$prev, na.rm = TRUE) <= 1)
    warning('\nis prevalence fractional? - should be percentage\n')

  # create meanz dz and LMS
  data <- data %>%
    mutate(from = !!from,
           to = !!to,
           age = round(age * 2) / 2,
           n = 1:n()) %>%
           select(age, n, everything()) %>%
    pivot_longer(c(from, to), values_to = 'cutoff', names_to = NULL) %>%
    mutate(z = diag(cutoffs[.data$cutoff, .data$sex]),
           cor = fct_relabel(fct_inorder(.data$cutoff), ~c(!!to, !!from))) %>%
    group_by(cutoff) %>%
    mutate(bmi = LMS2z(.data$age, .data$z, .data$sex, 'bmi', ref_sitar(cutoff), toz = FALSE),
           zrev = LMS2z(.data$age, .data$bmi, .data$sex, 'bmi', ref_sitar(cor))) %>%
    ungroup %>%
    mutate(cutoff = fct_relabel(fct_inorder(.data$cutoff), ~c('from', 'to')),
           zmean = (.data$z + .data$zrev) / 2,
           cor = NULL) %>%
    pivot_wider(names_from = .data$cutoff, values_from = .data$z:.data$zmean) %>%
    mutate(dz = .data$zmean_to - .data$zmean_from,
           prev_new = .data$dz + qnorm(prev / 100) * -sign(.data$z_from),
           prev_new = pnorm(.data$prev_new * -sign(.data$z_to)) * 100,
           from = !!from,
           to = !!to,
           n = NULL) %>%
    select(age, sex, starts_with('prev'), contains('z'), starts_with('bmi'), everything())

  # check for plot
  report <- match.arg(report)
  plot <- match.arg(plot)
  if (plot == 'density') {
    if (nrow(data) > 1L)
      cat('\ndensity plot is best with just one row of data\n')
    report <- 'longer'
  }
  else if (plot == 'compare') {
    stopifnot('compare plot needs prev_true set' = !any(is.na(data$prev_true)))
    report <- 'wider'
  }

  # format prevalence as vector or tibble
  data <- switch(report,
                 # vector
                 vector = data %>%
                   pull(.data$prev_new),
                 # wider
                 wider = data,
                 # longer
                 longer = data %>%
                   pivot_longer(ends_with(c('_from', '_to')),
                                names_to = c('.value', 'cutoff'), names_sep = '_') %>%
                   select(-c(from, to)) %>%
                   mutate(cutoff = fct_relabel(fct_inorder(.data$cutoff), ~c(!!from, !!to))))

  # return data or plot
  switch(plot,
         # return data
         no = return(data),
         # density plot
         density = {
           data %>%
             group_by(cutoff) %>%
             mutate(LMS = attr(LMS2z(.data$age, .data$z, .data$sex, 'bmi', ref_sitar(cutoff), toz = FALSE, LMStable = TRUE), 'LMStable')) %>%
             ungroup %>%
             rowwise %>%
             mutate(density = with(LMS, list(as_tibble(pdLMS(L, M, S, plot=F)[1:2])))) %>%
             unnest(density) %>%
             mutate(density = drop(density)) %>%
             ggplot(aes(.data$x, .data$density, colour = .data$cutoff)) +
             xlab(expression(body~mass~index~~(kg/m^2))) +
             geom_path() +
             geom_vline(aes(xintercept = .data$bmi, colour = .data$cutoff), data = data, linetype = 2)
         },
         # comparison plot
         compare = {
           # compare true and estimated prevalence (output with report = 'wider')
           ggplot(data, aes(prev_true, .data$prev_new, group = sex, colour = sex)) +
             xlab('observed prevalence (%)') + ylab('predicted prevalence (%)') +
             geom_point() +
             geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'gray') +
             scale_x_continuous(trans='log2') +
             scale_y_continuous(trans='log2')
         }
  )
}
