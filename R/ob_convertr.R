#' Convert between IOTF, WHO and CDC prevalence rates for child underweight, overweight and obesity
#'
#' Child thinness, overweight and obesity are defined by the child's body mass index (BMI) exceeding
#' (above or below) a pre-specified reference cutoff. Three references are considered:
#' IOTF (International Obesity Task Force), WHO (World Health Organization) and CDC
#' (US Centers for Disease Control and Prevention). They each have separate cutoffs
#' for thinness, overweight and obesity. \code{ob_convertr} takes age-sex-specific
#' prevalence rates based on one reference cutoff and estimates what the rates would be
#' if based on a different cutoff, using a novel algorithm.
#'
#' Each reference has cutoffs corresponding to respectively underweight, overweight and obesity.
#' They are "IOTF16", "IOTF17", "IOTF18.5", "IOTF25", "IOTF30", "IOTF35",
#' "WHO-2", "WHO-1", "WHO+1", "WHO+2", "CDC5", "CDC85" and "CDC95".
#'
#' See ... for an explanation of how the conversion calculation is done.
#'
#' @param prev vector of percentage prevalence rates in groups of children.
#' @param age vector of mean ages in years corresponding to each rate.
#' @param sex vector of the sexes corresponding to each rate, coded as either 'boys/girls'
#' or 'male/female' or '1/2' (only the first character in each case counts).
#' @param from name of the BMI cutoff (see below) on which \code{prev} is based.
#' @param to name of the BMI cutoff to be used to convert the prevalence to.
#' @param prev_true optional vector of known prevalence rates corresponding to cutoff \code{to},
#' for comparison purposes.
#' @param report character controlling the format of the returned data: 'vector'
#' for the estimated prevalence rates, 'wider' for the working tibble in wide format,
#' i.e. the 'from' and 'to' data side by side, and 'longer' for the tibble in
#' long format, i.e. two rows, one for 'from' and one for 'to'.
#' @param plot character controlling what if anything is plotted: 'no' for no
#' plot, 'density' to compare the BMI distributions and cutoffs corresponding
#' to \code{from} and \code{to}, or 'compare' to compare the estimated prevalence
#' rates with the true rates in \code{prev_true}.
#' @param data data frame containing prev, age, sex and optionally prev_true.
#'
#' @return Either the estimated prevalence rates or a plot
#' visualizing the findings, depending on the \code{report} and \code{plot} settings.
#'
#' Note that \code{from} and \code{to} cannot be the same, but they can be based on the
#' same reference, e.g. IOTF25 and IOTF30. They must be from the same tail of the distribution
#' (i.e. both thinness or both overweight/obesity) except when \code{plot = density}.
#'
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' ## convert 10% IOTF overweight prevalence (cutoff IOTF25) in 8-year-old boys
#' ## to the overweight prevalence based on WHO, i.e. cutoff WHO+1
#' ob_convertr(prev = 10, age = 8, sex = 'boys', from = 'IOTF25', to = 'WHO+1')
#'
#' ## compare the BMI density functions and cutoffs for IOTF25 and WHO+1
#' ob_convertr(10, 8, 'boys', from = 'IOTF25', to = 'WHO+1', plot = 'density')
#' @importFrom utils globalVariables
#' @importFrom forcats fct_inorder fct_collapse
#' @importFrom ggplot2 ggplot xlab ylab geom_path geom_point geom_vline geom_abline scale_x_continuous scale_y_continuous aes
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate rename filter transmute bind_rows bind_cols select across left_join pull contains starts_with ends_with if_else lead
#' @importFrom rlang .data as_label enquo
#' @importFrom stats pnorm qnorm
#' @importFrom tidyr pivot_wider pivot_longer drop_na
#' @export ob_convertr
ob_convertr <- function(prev = 50, age, sex, from, to, prev_true = NA, report = c('vector', 'wider', 'longer'),
                     plot = c('no', 'density', 'compare'), data = parent.frame()) {


  cutoffs <- tibble(ref = c(rep('IOTF', 6), rep('WHO', 4), rep('CDC', 3)),
                    cutoff = c(16:17, 18.5, 25, 30, 35, "-2", "-1", "+1", "+2", centile <- c(5, 85, 95)),
                    boys = c(-2.565, -1.877, -1.014, 1.31, 2.288, 2.93, -2:-1, 1:2, qnorm(centile / 100)),
                    girls = c(-2.436, -1.789, -0.975, 1.244, 2.192, 2.822, -2:-1, 1:2, qnorm(centile / 100))) %>%
    mutate(cutoff = paste0(.data$ref, .data$cutoff))

  # globalVariables(c('!!', '.', 'iotf', 'who0607', 'cdc2000'))

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

  from <- toupper(from)
  to <- toupper(to)
  stopifnot('`from` not recognised' = from %in% cutoffs$cutoff,
            '`to` not recognised' = to %in% cutoffs$cutoff,
            '`from` and `to` the same' = from != to)

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
    filter(.data$cutoff == from | .data$cutoff == to) %>%
    mutate(co = if_else(.data$cutoff == from, 'from', 'to')) %>%
    pivot_longer(c(.data$boys, .data$girls), names_to = 'sex', values_to = 'z')

  stopifnot('cannot compare lower and upper tail cutoffs' = sign(min(cutoffs$z)) == sign(max(cutoffs$z)) |
              plot == 'density')

  # LMS data by ref sex and age
  refdata <- sitar::iotf %>%
    mutate(ref = 'IOTF') %>%
    bind_rows(sitar::who0607 %>%
                select(c(.data$years, .data$sex, ends_with('bmi'))) %>%
                filter(lead(.data$years) > .data$years) %>%
                mutate(ref = 'WHO'),
              sitar::cdc2000 %>%
                select(c(.data$years, .data$sex, ends_with('bmi'))) %>%
                drop_na %>%
                mutate(years = round(.data$years * 24) / 24,
                       across(everything(), ~{(.x + lag(.x)) / 2}),
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
    mutate(age = round(age * 2) / 2) %>%
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
    pivot_longer(-c(age, sex, prev, prev_true), names_to = c('.value', 'co'), names_sep = '_') %>%
    mutate(zrev = ((.data$bmi / .data$M)^.data$L - 1) / (.data$L * .data$S)) %>%
    # swap bmi back again
    pivot_wider(names_from = .data$co, values_from = .data$ref:.data$zrev) %>%
    rename(bmi = .data$bmi_to,
           bmi_to = .data$bmi_from) %>%
    rename(bmi_from = .data$bmi) %>%
    # calculate delta z-score and adjust prevalence
    mutate(dz = (.data$z_to - .data$z_from + .data$zrev_from - .data$zrev_to) / 2 * -sign(.data$z_from),
           prev_new = qnorm(prev / 100),
           prev_new = pnorm(.data$prev_new + .data$dz) * 100) %>%
    select(c(.data$age, .data$sex, .data$prev, .data$prev_new, .data$prev_true, contains('z'),
                    starts_with('bmi'), everything()))

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
             xlab('true prevalence (%)') + ylab('estimated prevalence (%)') +
             geom_point() +
             geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'gray') +
             scale_x_continuous(trans='log2') +
             scale_y_continuous(trans='log2')
         }
  )
}
