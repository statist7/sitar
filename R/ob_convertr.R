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
#' The IOTF cutoffs correspond to the value of BMI (kg/m^2) at age 18: IOTF35
#' (morbid obesity), IOTF30 (obesity), IOTF25 (overweight), IOTF18.5 (grade 1
#' thinness), IOTF17 (grade 2 thinness) and IOTF16 (grade 3 thinness).
#'
#' The WHO cutoffs correspond to BMI z_scores. Age 5-19 years, WHO+2 (obesity),
#' WHO+1 (overweight) and WHO-2 (thinness). Age 0-5 years, WHO+3 (obesity),
#' WHO+2 (overweight) and WHO-2 (thinness).
#'
#' The CDC cutoffs correspond to BMI centiles: CDC95 (obesity), CDC85
#' (overweight) and CDC5 (thinness).
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
#' ## convert 10% IOTF overweight prevalence (cutoff IOTF25) in 8-year-old boys
#' ## to the overweight prevalence based on WHO, i.e. cutoff WHO+1
#' ob_convertr(prev = 10, age = 8, sex = 'boys', from = 'IOTF25', to = 'WHO+1')
#'
#' ## compare the BMI density functions and cutoffs for IOTF25 and WHO+1
#' ob_convertr(prev = 10, age = 8, sex = 'boys', from = 'IOTF25', to = 'WHO+1', plot = 'density')
#'
#'#' ## convert IOTF overweight prevalence to WHO overweight prevalence
#' ## and compare with true value - boys and girls age 7-17
#' ## note the need to first add obesity prevalence to overweight prevalence
#' data(deren)
#' deren <- within(deren, {
#'   IOTF25 = IOTF25 + IOTF30
#'   `WHO+1` = `WHO+1` + `WHO+2`})
#' ob_convertr(prev = IOTF25, age = Age, sex = Sex, from = 'IOTF25', to = 'WHO+1',
#'    prev_true = `WHO+1`, data = deren, plot = 'compare')

#' @importFrom forcats fct_inorder fct_collapse
#' @importFrom ggplot2 ggplot xlab ylab geom_path geom_point geom_vline
#'   geom_abline scale_x_continuous scale_y_continuous aes
#' @importFrom tibble tibble
#' @importFrom dplyr mutate rename filter transmute bind_rows bind_cols select
#'   across left_join pull contains starts_with ends_with if_else n
#' @importFrom rlang .data enquo
#' @importFrom stats pnorm qnorm
#' @importFrom tidyr pivot_wider pivot_longer drop_na
#' @export ob_convertr
ob_convertr <- function(prev = 50, age, sex, from, to, prev_true = NA, report = c('vector', 'wider', 'longer'),
                     plot = c('no', 'density', 'compare'), data = parent.frame()) {

  cutoffs <- tibble(ref = c(rep('IOTF', 6), rep('WHO', 5), rep('CDC', 3)),
                    cutoff = c(16:17, 18.5, 25, 30, 35, "-2", "-1", "+1", "+2", "+3", centile <- c(5, 85, 95)),
                    boys = c(-2.565, -1.877, -1.014, 1.31, 2.288, 2.93, -2:-1, 1:3, qnorm(centile / 100)),
                    girls = c(-2.436, -1.789, -0.975, 1.244, 2.192, 2.822, -2:-1, 1:3, qnorm(centile / 100))) %>%
    mutate(cutoff = paste0(.data$ref, .data$cutoff))

  # check sex contains only 1/M/B/TRUE or 2/F/G
  test_sex <- function(sex) {
    fsex <- toupper(substr(sex, 1, 1))
    fsex <- factor(fsex, levels = c(1:2, 'M', 'F', 'B', 'G', 'T'))
    levels(fsex) <- c(rep(1:2, 3), 1)
    fsex <- fct_collapse(fsex, boys = '1', girls = '2')
    droplevels(fsex)
  }

  # create tibble of age sex and prev
  if (identical(data, parent.frame())) {
    data <- eval.parent(substitute(
      data.frame(age = age, sex = test_sex(sex), prev = prev, prev_true = prev_true)))
  } else {
    prev <- enquo(prev)
    age <- enquo(age)
    sex <- enquo(sex)
    prev_true <- enquo(prev_true)
    data <- data %>%
      transmute(age = !!age, sex = test_sex(!!sex), prev = !!prev, prev_true = !!prev_true)
  }

  if (max(data$prev) <= 1)
    warning('\nis prevalence fractional? - should be percentage\n')

  from <- unique(toupper(from))
  to <- unique(toupper(to))
  stopifnot('`from` not recognised' = from %in% cutoffs$cutoff,
            '`to` not recognised' = to %in% cutoffs$cutoff,
            '`from` should be length 1' = length(from) == 1L,
            '`to` should be length 1' = length(to) == 1L)

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

  # z-score cutoffs by sex and ref
  cutoffs <- cutoffs %>%
    filter(.data$cutoff == from) %>%
    mutate(co = 'from') %>%
    bind_rows(cutoffs %>%
                filter(.data$cutoff == to) %>%
                mutate(co = 'to')) %>%
    pivot_longer(c(.data$boys, .data$girls), names_to = 'sex', values_to = 'z')

  # LMS data by ref sex and age
  refdata <- sitar::iotf %>%
    mutate(ref = 'IOTF') %>%
    bind_rows(sitar::who0607 %>%
                select(c(.data$years, .data$sex, ends_with('bmi'))) %>%
                filter(dplyr::lead(.data$years) > .data$years) %>%
                mutate(ref = 'WHO'),
              sitar::cdc2000 %>%
                select(c(.data$years, .data$sex, ends_with('bmi'))) %>%
                drop_na %>%
                mutate(years = round(.data$years * 24) / 24,
                       across(everything(), ~{(.x + dplyr::lag(.x)) / 2}),
                       ref = 'CDC')) %>%
    rename(age = .data$years,
           L = .data$L.bmi,
           M = .data$M.bmi,
           S = .data$S.bmi) %>%
    filter(age >= 2, age <= 18, age == (round(age * 2) / 2)) %>%
    mutate(across(c(sex, .data$ref), factor),
           sex = test_sex(.data$sex))

  # combine prevalence, z-score cutoffs and LMS data
  data <- left_join(data, cutoffs, by = "sex") %>%
    # round age to half-years
    # add count to ensure unique rows
    mutate(age = round(age * 2) / 2,
           n = (1:n() + 1) %/% 2) %>%
    select(age, n, everything()) %>%
    left_join(refdata, by = c("age", "sex", "ref"))

  stopifnot('no data' = nrow(data %>% select(-prev_true) %>% drop_na) > 0L)

  data <- data %>%
    # calculate bmi from LMS
    mutate(bmi = .data$M * (1 + .data$L * .data$S * .data$z)^(1/.data$L)) %>%
    # swap bmi from and to
    pivot_wider(names_from = .data$co, values_from = c(.data$ref:.data$cutoff, .data$z:.data$bmi)) %>%
    rename(bmi = .data$bmi_to,
           bmi_to = .data$bmi_from) %>%
    rename(bmi_from = .data$bmi) %>%
    # calculate z-scores for swapped bmi
    pivot_longer(-c(age:prev_true), names_to = c('.value', 'co'), names_sep = '_') %>%
    mutate(zrev = ((.data$bmi / .data$M)^.data$L - 1) / (.data$L * .data$S)) %>%
    # swap bmi back again
    pivot_wider(names_from = .data$co, values_from = .data$ref:.data$zrev) %>%
    rename(bmi = .data$bmi_to,
           bmi_to = .data$bmi_from) %>%
    rename(bmi_from = .data$bmi) %>%
    # calculate delta z-score and adjust prevalence
    mutate(dz = (.data$z_to - .data$z_from + .data$zrev_from - .data$zrev_to) / 2,
           prev_new = .data$dz + qnorm(prev / 100) * -sign(.data$z_from),
           prev_new = pnorm(.data$prev_new * -sign(.data$z_to)) * 100,
           n = NULL) %>%
    select(c(.data$age:.data$prev, .data$prev_new, .data$prev_true, contains('z'),
             .data$bmi_from, .data$bmi_to, everything()))

  # format prevalence as vector or tibble
  data <- switch(report,
                 # vector
                 vector = data %>%
                   pull(.data$prev_new),
                 # wider
                 wider = data,
                 # longer
                 longer = data %>%
                   pivot_longer(ends_with(c('_from', 'to')),
                                names_to = c('.value', 'co'), names_sep = '_')
  )

  # return data or plot
  switch(plot,
         # return data
         no = return(data),
         # density plot
         density = {
           with(with(data, pdLMS(L, M, S, plot = FALSE)), {
             dimnames(density)[[2]] = data$cutoff
             tibble(x) %>%
               bind_cols(as.data.frame(density))
           }) %>%
             pivot_longer(-.data$x, names_to = 'cutoff', values_to = 'density') %>%
             mutate(cutoff = fct_inorder(.data$cutoff)) %>%
             ggplot(aes(.data$x, .data$density, group = .data$cutoff)) +
             xlab(expression(body~mass~index~~(kg/m^2))) +
             geom_path(aes(colour = .data$cutoff)) +
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
