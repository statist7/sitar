#' Convert between IOTF, WHO and CDC prevalence rates for child thinness,
#' overweight and obesity
#'
#' Child thinness, overweight and obesity are defined as the child's body mass
#' index (BMI) lying beyond a pre-specified reference cutoff. Three references
#' are compared here: IOTF (International Obesity Task Force), WHO (World Health
#' Organization) and CDC (US Centers for Disease Control and Prevention), each
#' of which have their own cutoffs. \code{ob_convertr} takes age-sex-specific
#' prevalence rates of thinness, overweight or obesity based on one cutoff, and
#' converts them to rates based on a different cutoff, using a novel prediction
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
#' Note 1: the overweight category needs to be analysed as overweight prevalence plus
#' obesity prevalence. To predict overweight excluding obesity, first calculate predicted
#' overweight plus obesity then subtract predicted obesity.
#'
#' Note 2: the category labels are harmonised and not necessarily as originally
#' described.
#'
#' The conversion algorithm exploits the fact that all three references are
#' based on the LMS method, which allows prevalence to be converted to a common
#' BMI centile and z-score scale.
#'
#' The algorithm is actually two distinct algorithms, depending on the number
#' of prevalence rates used for the prediction. In the simpler case
#' (Cole and Lobstein, 2022) just one
#' rate is used, in which case the algorithm is commutative, meaning that
#' converting a prevalence rate from cutoff A to cutoff B and then from B to A
#' returns the original value. For this algorithm \code{from} and \code{to} are
#' the names of the cutoffs, and \code{pfrom} and optionally \code{pto} are vectors.
#'
#' The alternative algorithm  (Cole and Lobstein, 2023) uses two known prevalence rates,
#' typically overweight and obesity based on one reference, and returning the corresponding
#' rates based on another reference. Here \code{from} and \code{to} are
#' the names of the cutoffs as length-2
#' character strings, while \code{pfrom} and optionally \code{pto}
#' are character strings giving the names of the corresponding vector prevalence rates.
#' For convenience the \code{from} or \code{to} names 'CDC', 'IOTF' or 'WHO'
#' expand to the corresponding pairs of cutoffs for overweight and obesity,
#' e.g. 'CDC' expands to c('CDC 85', 'CDC 95').
#'
#' Alternatively the latter algorithm can be used to interpolate or extrapolate
#' to a specified z-score cutoff assuming the same reference for all cutoffs.
#' Here the values of \code{from} and \code{to} are numerical z-score cutoffs,
#' with at least two for \code{from}. See the final example.
#'
#' @param age vector of ages between 2 and 18 years corresponding to each rate.
#' @param sex vector of the sexes corresponding to each rate, coded as either
#'   'boys/girls' or 'male/female' or '1/2' (upper or lower case, and only the
#'   first character considered).
#' @param from name(s) of the BMI cutoff(s) (see Details) on which the prevalence
#' is based.
#' @param to name(s) of the BMI cutoff(s) (see Details) on which to base the
#' predicted prevalence.
#' @param pfrom vector of age-sex-specific percentage prevalence rates
#' based on \code{from}, or the names of two or more such prevalence rates.
#' @param pto optional vector of known percentage prevalence rates
#'   based on \code{to}, or the names of two or more such prevalence rates
#'   (for validation purposes).
#' @param report character controlling the format of the returned data: 'vector'
#'   for the estimated prevalence rates, 'wider' for the working tibble in wide
#'   format, i.e. the \code{from} and \code{to} data side by side, or 'longer'
#'   for the tibble in long format, i.e. two rows per rate, one for \code{from}
#'   and one for \code{to}.
#' @param plot character controlling what if anything is plotted: 'no' for no
#'   plot, 'density' to display the BMI density distributions and cutoffs
#'   corresponding to \code{from} and \code{to}, or 'compare' to display the
#'   predicted prevalence rates plotted against the observed rates in
#'   \code{pto}.
#' @param data data frame containing \code{pfrom}, \code{age}, \code{sex} and
#'   \code{pto}.
#'
#' @return Either the predicted prevalence rates or a plot visualizing the
#'   findings, depending on the \code{report} and \code{plot} settings.
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
#' \url{https://www.who.int/growthref/growthref_who_bull/en/}
#' @examples
#' ## convert 10% IOTF overweight prevalence (cutoff IOTF 25) in 8-year-old boys
#' ## to the overweight prevalence based on WHO, i.e. cutoff WHO +1
#' ob_convertr(age = 8, sex = 'boys', from = 'IOTF 25', to = 'WHO +1', pfrom = 10)
#'
#' ## compare the BMI density functions and cutoffs for IOTF 25 and WHO +1
#' ob_convertr(age = 8, sex = 'boys', from = 'IOTF 25', to = 'WHO +1', pfrom = 10, plot = 'density')
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
#' ## \code{ob_convertr2} which is more accurate than \code{ob_convertr}
#' ob_convertr2(age = Age, sex = Sex, from = 'IOTF', to = 'WHO',
#'   pfrom = c('IOTF25', 'IOTF30'), pto = c('WHO+1', 'WHO+2'), data = deren, report = 'wider')
#'
#' ## extrapolate WHO overweight and obesity prevalence (cutoffs +1 and +2)
#' ## to severe obesity prevalence based on cutoff +3
#' ob_convertr2(Age, Sex, from = 1:2, to = 3,
#'   pfrom = c('WHO+1', 'WHO+2'), data = deren, report = 'wider')
#'
#' @importFrom forcats fct_inorder fct_collapse fct_relabel
#' @importFrom ggplot2 ggplot xlab ylab geom_path geom_point geom_vline
#'   geom_abline scale_x_continuous scale_y_continuous aes
#' @importFrom tibble tibble column_to_rownames as_tibble
#' @importFrom dplyr mutate rename filter transmute bind_rows bind_cols select
#'   across left_join pull contains starts_with ends_with if_else n rowwise
#'   group_by ungroup rename_with c_across
#' @importFrom purrr map_dfc map_dfr
#' @importFrom rlang .data enquo quo quo_is_null sym
#' @importFrom stats pnorm qnorm lm cor coef
#' @importFrom tidyr pivot_wider pivot_longer drop_na all_of any_of everything
#' @export
ob_convertr <- function(age, sex, from, to, pfrom = 50, pto = NA, data = parent.frame(),
                        report = c('vector', 'wider', 'longer'),
                        plot = c('no', 'density', 'compare')) {

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

  # check from and to
  from <- unique(toupper(from))
  to <- unique(toupper(to))
  stopifnot('`from` should be length 1' = length(from) == 1L,
            '`to` should be length 1' = length(to) == 1L,
            '`to` is the same as `from`' = to != from,
            '`from` not recognised' = from %in% cutoffs$cutoff,
            '`to` not recognised' = to %in% cutoffs$cutoff)

  # create tibble of age sex and pfrom
  data <- if (identical(data, parent.frame())) {
    eval.parent(substitute(
      data.frame(age = age, sex = test_sex(sex), pfrom = pfrom, pto = pto)))
  } else {
    age <- enquo(age)
    sex <- enquo(sex)
    pfrom <- enquo(pfrom)
    pto <- enquo(pto)
    data %>%
      transmute(age = !!age, sex = test_sex(!!sex), pfrom = !!pfrom, pto = !!pto)
  }

  # check for out of range prevalences
  data <- data %>%
    mutate(pfrom = if_else(pfrom * (100 - pfrom) <= 0, as.numeric(NA), pfrom))

  if (length(na.omit(data$pfrom) > 0))
    if (max(data$pfrom, na.rm = TRUE) <= 1)
      warning('is prevalence fractional? - should be percentage\n')

  # create meanz and dz and epto
  data <- data %>%
    mutate(from = !!from,
           to = !!to,
           age = round(age * 2) / 2,
           n = 1:n()) %>%
           select(age, n, everything()) %>%
    pivot_longer(c(from, to), values_to = 'cutoff', names_to = NULL) %>%
    mutate(z = cutoffs$z[factor(paste(.data$cutoff, .data$sex), levels = cutoffs$cutoff_sex)],
           cor = fct_relabel(fct_inorder(.data$cutoff), ~c(!!to, !!from))) %>%
    group_by(.data$cutoff) %>%
    mutate(bmi = LMS2z(.data$age, .data$z, .data$sex, 'bmi', ref_sitar(.data$cutoff), toz = FALSE),
           zrev = LMS2z(.data$age, .data$bmi, .data$sex, 'bmi', ref_sitar(cor))) %>%
    ungroup %>%
    mutate(cutoff = fct_relabel(fct_inorder(.data$cutoff), ~c('from', 'to')),
           zmean = (.data$z + .data$zrev) / 2,
           cor = NULL) %>%
    pivot_wider(names_from = .data$cutoff, values_from = .data$z:.data$zmean) %>%
    mutate(dz = .data$zmean_to - .data$zmean_from,
           epto = .data$dz + qnorm(pfrom / 100) * -sign(.data$z_from),
           epto = pnorm(.data$epto * -sign(.data$z_to)) * 100,
           from = !!from,
           to = !!to,
           n = NULL) %>%
    select(age, sex, starts_with('pfrom'), contains('z'), starts_with('bmi'), everything())

  # check for plot
  report <- match.arg(report)
  plot <- match.arg(plot)
  if (plot == 'density') {
    if (nrow(data) > 1L)
      cat('\ndensity plot is best with just one row of data\n')
    report <- 'longer'
  }
  else if (plot == 'compare') {
    stopifnot('compare plot needs pto set' = !any(is.na(data$pto)))
    report <- 'wider'
  }

  # format prevalence as vector or tibble
  data <- switch(report,
                 vector = data %>%
                   pull(.data$epto),
                 wider = data,
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
             group_by(.data$cutoff) %>%
             mutate(LMS = attr(LMS2z(.data$age, .data$z, .data$sex, 'bmi', ref_sitar(.data$cutoff),
                                     toz = FALSE, LMStable = TRUE), 'LMStable')) %>%
             ungroup %>%
             rowwise %>%
             mutate(density = with(.data$LMS, list(as_tibble(pdLMS(L, M, S, plot=F)[c('x', 'density')])))) %>%
             unnest(.data$density) %>%
             mutate(density = drop(.data$density)) %>%
             ggplot(aes(.data$x, .data$density, colour = .data$cutoff)) +
             xlab(expression(body~mass~index~~(kg/m^2))) +
             geom_path() +
             geom_vline(aes(xintercept = .data$bmi, colour = .data$cutoff), data = data, linetype = 2)
         },
         # comparison plot
         compare = {
           # compare true and estimated prevalence (output with report = 'wider')
           ggplot(data, aes(.data$pto, .data$epto, colour = .data$sex)) +
             xlab('observed prevalence (%)') + ylab('predicted prevalence (%)') +
             geom_point() +
             geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'gray') +
             scale_x_continuous(trans='log2') +
             scale_y_continuous(trans='log2')
         }
  )
}

#' @rdname ob_convertr
#' @export
ob_convertr2 <- function(age, sex, from, to, pfrom, pto = NA,
                         data = parent.frame(),
                         report = c('vector', 'wider', 'longer')) {

  wic <- function(type = c("WHO", "IOTF", "CDC")) {
    type <- match.arg(type)
    switch(type,
           WHO = c('WHO +1', 'WHO +2'),
           IOTF = c('IOTF 25', 'IOTF 30'),
           CDC = c('CDC 85', 'CDC 95'))
  }

  to_p <- function(z) pnorm(-z) * 100
  to_z <- function(p) -qnorm(p / 100)

  # read data
  if (identical(data, parent.frame())) {
    data <- tibble(age = !!age, sex = !!sex) %>%
      bind_cols(as_tibble(t(pfrom)))
    age <- quo(data$age)
    sex <- quo(data$sex)
    pfrom <- names(pfrom)
  } else {
    age <- enquo(age)
    sex <- enquo(sex)
    stopifnot('`pfrom` not found' = all(pfrom %in% names(data)))
    if (any(!is.na(pto)))
      stopifnot('`pto` not found' = all(pto %in% names(data)))
  }

  stopifnot('`from` and `to` should both be either character or numeric' =
              identical(mode(from), mode(to)) && (is.character(to) || is.numeric(to)))

  # set up from, to, pto and epto
  pto_na <- any(is.na(pto))
  if (is.character(from)) {
    if (length(from) == 1)
      from <- wic(toupper(from))
    if (length(to) == 1)
      to <- wic(toupper(to))
    stopifnot('`from` should be length 2 or more' = length(from) > 1L,
              '`to` should be length 2 or more' = length(to) > 1L)
    if (pto_na) {
      pto = paste0(to, '_to')
      data[, pto] <- 50
    }
    epto <- paste0(to, '_pred')
  } else {
    if (pto_na) {
      pto = paste0('z', sprintf('%.2g', to), '_to')
      data[, pto] <- 50
    }
    epto <- paste0('z', sprintf('%.2g', to), '_pred')
  }

  stopifnot('`pto` not found' = all(pto %in% names(data)),
            '`pfrom` and `from` must be same length' = length(pfrom) == length(from),
            '`pto` and `to` must be same length' = length(pto) == length(to))

  # set up zfrom and zto
  if (is.character(from)) {
    stopifnot('age and/or sex missing' = !(quo_is_null(age) | quo_is_null(sex)))
    data <- data %>%
      bind_cols({
        map_dfc(seq_along(from), ~ {
          dotx <- .x # ? bug not recognising .x in rename_with()
          ob_convertr(!!age, !!sex, from = from[.x], to = to[.x], !!sym(pfrom[.x]),
                      data = data, report = 'wider') %>%
            select(zmean_from, zmean_to) %>%
            rename_with(~ paste0(c('zfrom', 'zto'), dotx))
        })
      })
  } else {
    data[, paste0('zfrom', seq_along(from))] <- t(from)
    data[, paste0('zto', seq_along(to))] <- t(to)
  }

  # fit regression models
  data <- if (identical(from, to)) {
    data[, epto] <- data[, pfrom]
    data$rc <- data$rc_all <- 1
    data
  } else {
    weights <- c(rep(1, length(from)), rep(0, length(to)))
    models_all <- rc <- rc_all <- r <- NULL # to satisfy R CMD CHECK
    data %>%
      rowwise %>%
      mutate(x = list(c_across(c(starts_with(c('zfrom', 'zto'))))),
             y = list(c_across(all_of(c(pfrom, pto)))),
             y = list(if_else(.data$y * (100 - .data$y) > 0, .data$y,
                              as.numeric(NA))), # out of range prevalences NA
             models = list(lm(to_z(.data$y) ~ x, weights = weights)),
             rc = coef(.data$models)[2],
             models_all = list(lm(to_z(.data$y) ~ x)),
             rc_all = coef(.data$models_all)[2],
             r = cor(.data$models_all$model)[1, 2]) %>%
      bind_cols(map_dfr(.$models, ~{to_p(fitted(.x))[!weights]}) %>%
                  rename_with(~all_of(epto))) %>%
      mutate(across(all_of(epto), ~if_else(is.na(rc), as.numeric(NA), .x))) %>%
      ungroup
  }
  if (pto_na)
    data <- data %>%
    select(-c(all_of(pto), models_all, rc_all, r))

  # report
  report <- match.arg(report)
  switch(report,
         vector = data %>%
           select(any_of(epto)),
         wider = data %>%
           select(1:2, any_of(epto), all_of(pfrom), any_of(pto), rc),
         longer = data %>%
           select(1:2, any_of(epto), all_of(pfrom), any_of(pto), rc, everything()))
}

