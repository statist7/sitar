#' Deren prevalence data on child thinness, overweight and obesity
#'
#' Age-sex-specific prevalence rates of thinness, overweight and obesity in
#' Ukraine children based on body mass index and IOTF, WHO and CDC cut-offs.
#'
#' Note that the overweight prevalences are for overweight excluding obesity,
#' i.e. the prevalence for BMI between the overweight and obesity cutoffs.
#'
#' @name deren
#' @docType data
#' @format A tibble with 22 observations on the following 11 variables:
#'   \describe{
#'   \item{Age}{postnatal age from 7 to 17 completed years}
#'   \item{Sex}{two-level factor - Boys and Girls}
#'   \item{N}{integer - group sample size}
#'   \item{IOTF18.5}{thinness prevalence based on IOTF reference and 18.5 cutoff}
#'   \item{WHO-2}{thinness prevalence based on WHO reference and -2 cutoff}
#'   \item{CDC5}{thinness prevalence based on CDC reference and 5 cutoff}
#'   \item{IOTF25}{overweight prevalence based on IOTF reference and 25 cutoff}
#'   \item{WHO+1}{overweight prevalence based on WHO reference and +1 cutoff}
#'   \item{CDC85}{overweight prevalence based on CDC reference and 85 cutoff}
#'   \item{IOTF30}{obesity prevalence based on IOTF reference and 30 cutoff}
#'   \item{WHO+2}{obesity prevalence based on WHO reference and +2 cutoff}
#'   \item{CDC95}{obesity prevalence based on CDC reference and 95 cutoff} }
#' @references Deren K, Wyszynska J, Nyankovskyy S, Nyankovska O, Yatsula M,
#' Luszczki E, Sobolewski M, Mazur A. 2020. Assessment of body mass index in a
#' pediatric population aged 7-17 from Ukraine according to various
#' international criteria-A cross-sectional study. PLoS ONE 15.
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
#' Health Organization 2007; 85: 660-7.
#'
#' The CDC reference for children aged 2-20 years is: Must A, Dallal GE, Dietz
#' WH. Reference data for obesity: 85th and 95th percentiles of body mass index
#' (wt/ht2) and triceps skinfold thickness. American Journal of Clinical
#' Nutrition 1991; 53: 839-46.

#' @source The values are obtained from Table 2 of Deren et al (2020), recalculated to full accuracy.
#'   \url{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244300}.
#' @keywords datasets
#' @examples
#' ## convert IOTF obesity prevalence to WHO obesity prevalence
#' ## and compare with true WHO obesity prevalence - boys and girls age 7-17
#' data(deren)
#'   ob_convertr(age = Age, sex = Sex, from = 'IOTF 30', to = 'WHO +2',
#'     pfrom = IOTF30, pto = `WHO+2`, data = deren, plot = 'compare')
"deren"
