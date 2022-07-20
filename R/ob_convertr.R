#' Convert between IOTF, WHO and CDC prevalence rates for child thinness,
#' overweight and obesity
#'
#' Child thinness, overweight and obesity are defined as the child's body mass
#' index (BMI) lying beyond a pre-specified reference cutoff. Three references
#' are compared: IOTF (International Obesity Task Force), WHO (World Health
#' Organization) and CDC (US Centers for Disease Control and Prevention), each
#' of which have their own cutoffs. \code{ob_convertr} takes age-sex-specific
#' prevalence rates of thinness, overweight or obesity based on one of the cutoffs,
#' and converts them to the corresponding rates based on a different cutoff.
#' \code{ob_convertr2} uses paired prevalence rates of overweight and obesity on
#' one cutoff to estimate those based on another cutoff.
#'
#' The IOTF cutoffs correspond to the value of BMI
#' (\ifelse{html}{\out{kg/m<sup>2</sup>}}{\eqn{kg/m^2}}) at age 18: IOTF 35
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
#' Note: the overweight category needs to be analysed as overweight prevalence plus
#' obesity prevalence, i.e. the prevalence above the overweight cutoff. To predict
#' overweight prevalence excluding obesity prevalence, first calculate predicted
#' overweight prevalence including obesity then subtract predicted obesity prevalence.
#'
#' The algorithms for \code{ob_convertr} and \code{ob_convertr2} are distinguished
#' by the number of prevalence rates used for the prediction. For \code{ob_convertr}
#' (Cole & Lobstein, 2022) just one
#' rate is used -- in this case the algorithm is commutative, meaning that
#' converting a prevalence rate from cutoff A to cutoff B and then from B to A
#' returns the original value. \code{from} and \code{to} are
#' the names of the cutoffs, and \code{pfrom} and optionally \code{pto} are vectors
#' of percentage prevalence rates.
#'
#' \code{ob_convertr2} uses two known prevalence rates (Cole & Lobstein, 2023),
#' typically overweight and obesity based on one reference, returning the corresponding
#' rates based on another reference. It is more accurate than \code{ob_convertr} though
#' not exactly commutative. \code{from} and \code{to} are the names of the cutoffs as length-2
#' character strings, while \code{pfrom} and optionally \code{pto} are length-2
#' character strings giving the names of the corresponding vector prevalence rates.
#' For convenience the \code{from} or \code{to} names 'CDC', 'IOTF' or 'WHO'
#' expand to the corresponding pairs of cutoffs for overweight and obesity,
#' e.g. 'CDC' expands to c('CDC 85', 'CDC 95').
#'
#' Alternatively \code{ob_convertr2} can be used to interpolate or extrapolate
#' to one or more specified z-score cutoffs assuming the same reference for all cutoffs.
#' Here the values of \code{from} and \code{to} are numerical z-score cutoffs,
#' with at least two for \code{from}. See the final example.
#'
#' The algorithms require the prevalences of obesity and overweight net of obesity
#' to be non-zero, and if they are zero they are set to missing.
#'
#' @param age vector of ages between 2 and 18 years corresponding to prevalence rates \code{pfrom}.
#' @param sex vector of sexes corresponding to \code{pfrom}, coded as either
#'   'boys/girls' or 'male/female' or '1/2' (upper or lower case, based on the
#'   first character).
#' @param from name(s) of the BMI cutoff(s) on which the prevalence \code{pfrom}
#' is based (see Details).
#' @param to name(s) of the BMI cutoff(s) on which to base the
#' predicted prevalence (see Details).
#' @param pfrom vector of age-sex-specific percentage prevalence rates
#' based on \code{from} (\code{ob_convertr}) or the names of two or more such
#' prevalence rates (\code{ob_convertr2}).
#' @param pto vector (needed for \code{plot = "compare"}) of known percentage prevalence rates
#'   based on \code{to} (\code{ob_convertr}) or the names of two or more such
#'   prevalence rates (\code{ob_convertr2}).
#' @param report character controlling the format of the returned data: 'vector'
#'   for the estimated prevalence rates, 'wider' for the working tibble in wide
#'   format, i.e. the \code{from} and \code{to} data side by side, or 'longer'
#'   for the tibble in long format, i.e. two rows per rate, one for \code{from}
#'   and one for \code{to}. For \code{ob_convertr2} the three settings return
#'   progressively more information.
#' @param plot character controlling what if anything is plotted: 'no' for no
#'   plot, 'density' to display the BMI density distributions and cutoffs
#'   corresponding to \code{from} and \code{to}, or 'compare' to display the
#'   predicted prevalence rates plotted against the observed rates (\code{pto}).
#' @param data optional data frame containing \code{age}, \code{sex}, \code{pfrom} and
#'   \code{pto}.
#'
#' @return The predicted prevalence rates, optionally with a plot visualizing the
#'   findings, depending on the \code{report} and \code{plot} settings. Each
#'   predicted rate is given the name of the relevant cutoff followed by " pred".
#'
#'   With \code{report} set to "wider" or "longer", extra information
#'   is returned reflecting the internal workings of the algorithms. In particular
#'   \code{ob_convertr2} returns \code{b} the regression coefficient of z-score
#'   prevalence on z-score cutoff as described in Cole & Lobstein (2023).
#'
#'   If a \code{plot} is selected, the underlying data and plot are returned
#'   invisibly with names \code{data} and \code{plot}.
#'
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @references
#' Cole TJ, Lobstein T. Exploring an algorithm to harmonize International Obesity
#' Task Force and World Health Organization child overweight and obesity prevalence
#' rates. Pediatr Obes 2022;In press.
#'
#' Cole TJ, Lobstein T. An improved algorithm to harmonize child overweight and
#' obesity prevalence rates. Pediatr Obes 2023;In press.
#'
#' The CDC reference for children aged 2-20 years is: Must A, Dallal GE, Dietz
#' WH. Reference data for obesity: 85th and 95th percentiles of body mass index
#' (wt/ht2) and triceps skinfold thickness. American Journal of Clinical
#' Nutrition 1991; 53: 839-46. Available at:
#' \url{https://academic.oup.com/ajcn/article/53/4/839/4715058}
#'
#' The IOTF reference for children aged 2-18 years is: Cole TJ, Bellizzi MC,
#' Flegal KM, Dietz WH. Establishing a standard definition for child overweight
#' and obesity worldwide: international survey. BMJ 2000; 320: 1240-5. Available
#' at \doi{10.1136/bmj.320.7244.1240}
#'
#' The WHO reference for children aged 0-5 years is: WHO Child Growth Standards:
#' Length/height-for-age, weight-for-age, weight-for-length, weight-for-height
#' and body mass index-for-age: Methods and development. Geneva: World Health
#' Organization, 2006. Available at:
#' \url{https://www.who.int/toolkits/child-growth-standards/standards}
#'
#' The WHO reference for children aged 5-19 years is: de Onis M, Onyango AW,
#' Borghi E, Siyam A, Nishida C, Siekmann J. Development of a WHO growth
#' reference for school-aged children and adolescents. Bulletin of the World
#' Health Organization 2007; 85: 660-7. Available at:
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2636412/pdf/07-043497.pdf/}
#' @examples
#' ## convert 10% IOTF overweight prevalence (cutoff IOTF 25, including obesity)
#' ## in 8-year-old boys to overweight prevalence for cutoff WHO +1
#' ob_convertr(age = 8, sex = 'boys', from = 'IOTF 25', to = 'WHO +1', pfrom = 10)
#'
#' ## compare the BMI density functions and cutoffs for IOTF and WHO
#' ## in 8-year-old boys
#' ob_convertr2(age = 8, sex = 'boys', from = 'IOTF', to = 'WHO', plot = 'density')
#'
#' ## convert IOTF overweight prevalence to WHO overweight prevalence
#' ## and compare with true value - boys and girls aged 7-17 (22 groups)
#' ## note the need to first add obesity prevalence to overweight prevalence
#' data(deren)
#' deren <- within(deren, {
#'   CDC85 = CDC85 + CDC95
#'   IOTF25 = IOTF25 + IOTF30
#'   `WHO+1` = `WHO+1` + `WHO+2`})
#' ob_convertr(age = Age, sex = Sex, from = 'IOTF 25', to = 'WHO +1',
#'   pfrom = IOTF25, pto = `WHO+1`, data = deren, plot = 'compare')
#'
#' ## convert IOTF overweight and obesity prevalence to WHO using
#' ## ob_convertr2 - which is more accurate than ob_convertr
#' ob_convertr2(age = Age, sex = Sex, from = 'IOTF', to = 'WHO',
#'   pfrom = c('IOTF25', 'IOTF30'), pto = c('WHO+1', 'WHO+2'),
#'   data = deren, plot = 'compare')
#'
#' ## extrapolate WHO overweight and obesity prevalence (cutoffs +1 and +2)
#' ## to severe obesity prevalence based on cutoffs +2.5 or +3
#' ob_convertr2(Age, Sex, from = 1:2, to = c(2.5, 3),
#'   pfrom = c('WHO+1', 'WHO+2'), data = deren, report = 'wider')
#'
#' @importFrom forcats fct_inorder fct_collapse fct_relabel
#' @importFrom ggplot2 ggplot xlab ylab geom_path geom_point geom_vline
#'   geom_abline labs scale_x_continuous scale_y_continuous aes
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate rename filter transmute bind_cols select
#'   across pull contains starts_with ends_with if_else n rowwise
#'   group_by ungroup rename_with c_across matches
#' @importFrom purrr map_chr map_dfc map_dfr
#' @importFrom rlang .data enquo quo_name
#' @importFrom stats pnorm qnorm lm cor coef
#' @importFrom tidyr pivot_wider pivot_longer all_of any_of everything
#' @importFrom glue glue
#' @export
ob_convertr <- function(age, sex, from, to, pfrom = NA, pto = NA, data = parent.frame(),
                        report = c('vector', 'wider', 'longer'),
                        plot = c('no', 'density', 'compare')) {

  # longer data format
  data.longer <- function(data) {
    data %>%
      pivot_longer(ends_with(c('_from', '_to')),
                   names_to = c('.value', 'cutoff'), names_sep = '_') %>%
      select(-c(from, to)) %>%
      mutate(cutoff = fct_relabel(fct_inorder(.data$cutoff), ~c(from, to))) %>%
      select(.data$age, .data$sex, .data$cutoff, .data$epto, starts_with('p'),
             .data$dz, .data$zmean, contains('z'), everything())
  }

  # create cutoffs matrix
  cutoffs <- tibble(ref = c('CDC', 'IOTF', 'WHO'),
                    cutoff = list(c(5, 85, '95'),
                                  c(16:17, 18.5, 25, 30, '35'),
                                  c(-2:-1, paste0('+', 1:3))),
                    sex = list(c('boys', 'girls')),
                    f = list(function(cutoff, sex) qnorm(cutoff / 100),
                             function(cutoff, sex) LMS2z(18, cutoff, sex, 'bmi', 'iotf'),
                             function(cutoff, sex) cutoff)) %>%
    unnest(.data$cutoff) %>%
    unnest(.data$sex) %>%
    rowwise %>%
    mutate(z = .data$f(as.numeric(.data$cutoff), .data$sex), # apply function to cutoffs
           cutoff = paste(.data$ref, .data$cutoff),
           cutoff_sex = paste(.data$cutoff, .data$sex))

  # check from and to, age and sex
  from <- unique(toupper(from))
  to <- unique(toupper(to))
  stopifnot('`from` should be length 1' = length(from) == 1L,
            '`to` should be length 1' = length(to) == 1L,
            '`to` is the same as `from`' = to != from,
            '`from` not recognised' = from %in% cutoffs$cutoff,
            '`to` not recognised' = to %in% cutoffs$cutoff,
            '`age` or `sex` missing' = !missing(age) && !missing(sex))

  # create tibble of age sex and pfrom
  data <- if (identical(data, parent.frame())) {
    tibble(age = {{age}}, sex = test_sex({{sex}}), from = {{from}}, to = {{to}},
           pfrom = {{pfrom}}, pto = {{pto}})
  } else {
    data %>%
      transmute(age = {{age}}, sex = test_sex({{sex}}), from = {{from}}, to = {{to}},
                pfrom = {{pfrom}}, pto = {{pto}})
  }
  data <- data %>%
    select(pfrom, where(~!all(is.na(.))))
  stopifnot('`age` or `sex` missing' = all(c('age', 'sex') %in% names(data)))

  # check for out of range prevalences
  if (length(na.omit(data$pfrom)) > 0) {
    if (with(data, min(pfrom * (100 - pfrom), na.rm = TRUE) <= 0)) {
      warning('prevalences <= 0 or >= 100 set missing')
      data <- data %>%
        mutate(pfrom = if_else(pfrom * (100 - pfrom) <= 0, as.numeric(NA), pfrom))
    }
    if (length(na.omit(data$pfrom)) > 0)
      if(max(data$pfrom, na.rm = TRUE) <= 1)
        warning('is prevalence fractional? - should be percentage\n')
  }

  # create meanz and dz and epto
  data <- data %>%
    mutate(age = round(age * 2) / 2,
           age = if_else(age < 2 | age > 18, as.numeric(NA), age),
           n = 1:n()) %>%
           select(.data$age, .data$n, everything()) %>%
    pivot_longer(c(from, to), values_to = 'cutoff', names_to = NULL) %>%
    mutate(z = cutoffs$z[factor(paste(.data$cutoff, .data$sex), levels = cutoffs$cutoff_sex)],
           cor = fct_relabel(fct_inorder(.data$cutoff), ~c(to, from))) %>%
    group_by(.data$cutoff) %>%
    mutate(bmi = LMS2z(.data$age, .data$z, .data$sex, 'bmi', ref_sitar(.data$cutoff), toz = FALSE),
           zrev = LMS2z(.data$age, .data$bmi, .data$sex, 'bmi', ref_sitar(cor))) %>%
    ungroup %>%
    mutate(cutoff = fct_relabel(fct_inorder(.data$cutoff), ~c('from', 'to')),
           zmean = (.data$z + .data$zrev) / 2,
           cor = NULL) %>%
    pivot_wider(names_from = .data$cutoff, values_from = .data$z:.data$zmean) %>%
    mutate(from = from,
           to = to,
           dz = .data$zmean_to - .data$zmean_from,
           epto = .data$dz + qnorm(pfrom / 100) * -sign(.data$z_from),
           epto = pnorm(.data$epto * -sign(.data$z_to)) * 100,
           n = NULL) %>%
    select(.data$age, .data$sex, .data$from, .data$to, .data$epto, starts_with('p'), .data$dz,
           starts_with('zmean'), contains('z'), everything())

  # return data ± plot
  switch(match.arg(plot),
         # return data only
         no = {
           data <- switch(match.arg(report),
                  vector = return(data %>%
                                    pull(.data$epto)),
                  wider = data,
                  longer = data.longer(data)
                  )
           data <- data %>%
             rename_with(~paste(to, 'pred'), .cols = .data$epto) %>%
             rename_with(~quo_name(enquo(pfrom)), .cols = pfrom)
           if ('pto' %in% names(data))
             data <- data %>%
             rename_with(~quo_name(enquo(pto)), .cols = pto)
           return(data)
           },
         # density plot plus data in longer format
         density = {
           data <- data.longer(data) %>%
             group_by(.data$cutoff) %>%
             mutate(LMS = attr(LMS2z(.data$age, .data$z, .data$sex, 'bmi', ref_sitar(.data$cutoff),
                                     toz = FALSE, LMStable = TRUE), 'LMStable')) %>%
             ungroup %>%
             rowwise %>%
             mutate(density = with(.data$LMS,
                                   list(as_tibble(pdLMS(L, M, S, N = 200, plot = FALSE)[c('x', 'density')])))) %>%
             unnest(.data$density) %>%
             mutate(density = drop(.data$density)) %>%
             select(.data$x, .data$density, .data$cutoff, .data$bmi)
           plot <- data %>%
             ggplot(aes(.data$x, .data$density, colour = .data$cutoff)) +
             xlab(expression(body~mass~index~~(kg/m^2))) +
             geom_path() +
             geom_vline(aes(xintercept = .data$bmi, colour = .data$cutoff),
                        data = data %>% nest(data = c(.data$x, .data$density)), linetype = 2)
           print(plot)
           invisible(list(data = data, plot = plot))
         },
         # comparison plot plus data in wider format
         compare = {
           stopifnot('compare plot needs pto set' = !any(is.na(data$pto)))
           data <- data %>%
             select(.data$age, .data$sex, .data$pfrom, .data$pto, .data$epto)
           plot <- ggplot(data, aes(.data$pto, .data$epto, shape = .data$sex)) +
             xlab('observed prevalence (%)') +
             ylab('predicted prevalence (%)') +
             labs(caption = glue('prevalence for {to} predicted from {from}')) +
             geom_point() +
             geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'gray')
           print(plot)
           invisible(list(data = data, plot = plot))
         }
  )
}

utils::globalVariables("where")

#' @rdname ob_convertr
#' @export
ob_convertr2 <- function(age, sex, from, to, pfrom = NA, pto = NA,
                         data = parent.frame(),
                         report = c('vector', 'wider', 'longer'),
                         plot = c('no', 'density', 'compare')) {

  wic <- function(x) {
    if (length(x) > 1)
      return(x)
    x <- match.arg(x, c("WHO", "IOTF", "CDC"))
    switch(x,
           WHO = c('WHO +1', 'WHO +2'),
           IOTF = c('IOTF 25', 'IOTF 30'),
           CDC = c('CDC 85', 'CDC 95'))
  }

  to_p <- function(z) pnorm(-z) * 100

  to_z <- function(p) -qnorm(p / 100)

  co2ref <- function(x) map_chr(strsplit(x, ' '), ~.x[1])

  stopifnot('`age` or `sex` missing' = !missing(age) && !missing(sex))
  if (any(is.na(pfrom))) {
    stopifnot('`pfrom` not found' = match.arg(plot) == 'density')
    pfrom <- setNames(2:1*30, wic(from)) # distinct dummy values
  }
  if (any(is.na(pto)))
    pto <- character()
  if (identical(data, parent.frame())) {
    if (is.null(names(pfrom)))
      names(pfrom) <- wic(from)
    data <- tibble(age = {{age}}, sex = {{sex}}) %>%
      bind_cols(as_tibble(t(pfrom))) %>%
      bind_cols(as_tibble(t(pto)))
    pfrom <- names(pfrom)
    pto <- names(pto)
  } else {
    stopifnot('`pfrom` not found' = all(pfrom %in% names(data)),
              '`pto` not found' = all(pto %in% names(data)))
    data <- data %>%
      select(age = {{age}}, sex = {{sex}}, pfrom, pto) %>%
      mutate(sex = test_sex(.data$sex))
  }

  # check for equal pfrom values
  data <- data %>%
    mutate(across(pfrom[2], ~if_else(get(pfrom[2]) == get(pfrom[1]),
                                 as.numeric(NA), .x)))

  stopifnot('`from` and `to` should both be either character or numeric' =
              identical(mode(from), mode(to)) && (is.character(to) || is.numeric(to)))

  # set up from, to, pto and epto
  pto_na <- length(pto) == 0L
  if (is.character(from)) {
    from <- wic(toupper(from))
    to <- wic(toupper(to))
    stopifnot('`from` should be length 2 or more' = length(from) > 1L,
              '`to` should be length 2 or more' = length(to) > 1L)
    if (pto_na) {
      pto = paste(to, 'to')
      data[, pto] <- t(50 + seq_along(to))
    }
    epto <- paste(to, 'pred')
  } else {
    if (pto_na) {
      pto = paste0('z', sprintf('%.2g', to), '_to')
      data[, pto] <- t(50 + seq_along(to))
    }
    epto <- paste0('z', sprintf('%.2g', to), '_pred')
  }

  stopifnot('`pto` not found' = all(pto %in% names(data)),
            '`pfrom` and `from` must be same length' = length(pfrom) == length(from),
            '`pto` and `to` must be same length' = length(pto) == length(to))

  # set up zfrom and zto
  if (is.character(from)) {
    data <- data %>%
      bind_cols({
        map_dfc(seq_along(from), ~ {
          dotx <- .x # ? bug not recognising .x in rename_with()
          ob_convertr(age, sex, from = from[.x], to = to[.x], get(pfrom[.x]),
                      data = data, report = 'wider') %>%
            select(.data$z_from, .data$z_to, .data$zmean_from, .data$zmean_to) %>%
            rename_with(~ paste0(c('z0from', 'z0to', 'zfrom', 'zto'), dotx))
        })
      })
  } else {
    data[, paste0('zfrom', seq_along(from))] <- t(from)
    data[, paste0('zto', seq_along(to))] <- t(to)
  }

  # fit regression models
  data <- if (identical(from, to)) {
    data[, epto] <- data[, pfrom]
    data$b <- data$b_all <- 1
    data
  } else {
    weights <- c(rep(1, length(from)), rep(0, length(to)))
    data %>%
      rowwise %>%
      mutate(z = list(c_across(c(starts_with(c('z0from', 'z0to'))))),
             x = list(c_across(c(starts_with(c('zfrom', 'zto'))))),
             y = list(c_across(all_of(c(pfrom, pto)))),
             y = list(if_else(.data$y * (100 - .data$y) > 0, .data$y,
                              as.numeric(NA))), # out of range prevalences NA
             models = list(lm(to_z(.data$y) ~ x, weights = weights)),
             b = coef(.data$models)[2],
             models_all = list(lm(to_z(.data$y) ~ x)),
             b_all = coef(.data$models_all)[2],
             r = cor(.data$models_all$model)[1, 2]) %>%
      bind_cols(map_dfr(.$models, ~{to_p(fitted(.x))[!weights]}) %>%
                  rename_with(~all_of(epto))) %>%
      select(-matches(c('from.', 'to.'))) %>%
      mutate(across(all_of(epto), ~if_else(is.na(.data$b), as.numeric(NA), .x))) %>%
      ungroup
  }
  if (pto_na)
    data <- data %>%
    select(-c(all_of(pto), .data$models_all, .data$b_all, .data$r))

  # return data ± plot
  switch(match.arg(plot),
         # return data only
         no = {
           data <- switch(match.arg(report),
                          vector = data %>%
                            select(any_of(epto)),
                          wider = data %>%
                            select(1:2, any_of(epto), all_of(pfrom), any_of(pto), .data$b),
                          longer = data %>%
                            select(1:2, any_of(epto), all_of(pfrom), any_of(pto), .data$b, everything()))
           return(data)
         },
         # density plot plus data in longer format
         density = {
           data <- data %>%
             select(1:2, .data$z) %>%
             unnest(.data$z) %>%
             mutate(cutoff = c(from, to),
                    reference = fct_inorder(factor(co2ref(.data$cutoff))),
                    level = c(rep(c('overweight', 'obesity'), 2)),
                    level = fct_inorder(.data$level)) %>%
             group_by(.data$reference) %>%
             mutate(bmi = LMS2z(.data$age, .data$z, .data$sex, 'bmi', ref_sitar(.data$reference),
                                toz = FALSE, LMStable = TRUE),
                    LMS = attr(.data$bmi, 'LMStable')) %>%
             ungroup %>%
             rowwise %>%
             mutate(density = with(.data$LMS,
                                   list(as_tibble(pdLMS(L, M, S, N = 200, plot = FALSE)[c('x', 'density')])))) %>%
             unnest(.data$density) %>%
             mutate(density = drop(.data$density)) %>%
             select(.data$x, .data$density, .data$reference, .data$level, .data$cutoff, .data$bmi)
           plot <- data %>%
             ggplot(aes(.data$x, .data$density, group = .data$reference, colour = .data$reference)) +
             xlab(expression(body~mass~index~~(kg/m^2))) +
             geom_path(data = data %>% filter(.data$level == 'overweight')) +
             geom_vline(aes(xintercept = .data$bmi, linetype = .data$level, colour = .data$reference),
                        data = data %>% nest(data = c(.data$x, .data$density)))
           print(plot)
           invisible(list(data = data, plot = plot))
         },
         # comparison plot plus data in wider format
         compare = {
           stopifnot('compare plot needs pto set' = !pto_na)
           owob <- c('overweight', 'obesity')
           data <- data %>%
             select(-(.data$z:.data$r)) %>%
             rename_with(~c(t(outer(c('pfrom', 'pto', 'epto'), owob,
                                    paste, sep = '_'))), .cols = c(pfrom, pto, epto)) %>%
             pivot_longer(ends_with(owob), names_to = c('.value', 'level'), names_sep = '_') %>%
             mutate(level = fct_inorder(.data$level))
           plot <- ggplot(data, aes(.data$pto, .data$epto, colour = .data$level, shape = .data$sex)) +
             xlab('observed prevalence (%)') +
             ylab('predicted prevalence (%)') +
             labs(caption = glue('prevalence for {to[1]} / {to[2]} predicted from {from[1]} / {from[2]}')) +
             geom_point() +
             geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'gray')
           print(plot)
           invisible(list(data = data, plot = plot))
         }
  )
}


