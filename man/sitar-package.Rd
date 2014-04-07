\name{sitar-package}
\alias{sitar-package}
\docType{package}
\title{SITAR (SuperImposition by Translation And Rotation) growth curve analysis}
\description{SITAR is a method of growth curve analysis, based on \pkg{nlme},
that estimates a single mean growth curve as a regression spline, plus a set of
up to three fixed and random effects (a, b and c) defining how individual growth
curves differ from the mean curve. SITAR stands for SuperImposition by
Translation And Rotation.
 
 The package also contains some utility functions for the LMS method, as used to
 construct growth reference centiles (see \pkg{gamlss}).
}
\details{
\tabular{ll}{
Package: \tab sitar\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-23\cr
License: \tab GPL-2\cr
}
Effect a (or alpha) measures size, and is a random intercept relative to the
 spline curve intercept. Effect b (or beta) measures tempo, the timing of the growth
 process, and reflects a shift on the x scale relative to the mean. Effect c (or gamma)
 is velocity, and indicates how the x scale is stretched or shrunk reflecting the rate
 at which 'time' passes for individuals. The aim is for individual curves, adjusted for
 a, b and c, to lie on top of the mean curve. 
 
 The package creates an object of class \code{sitar}, based on \pkg{nlme},
 representing the nonlinear mixed-effects model fit. Generic functions such as
 \code{print}, \code{plot} and \code{summary} have methods to show the results
 of the fit. The functions \code{resid}, \code{coef}, \code{fitted},
 \code{fixed.effects}, and \code{random.effects} can be used to extract some of
 its components. The functions \code{devadj} and \code{varexp} compare
 respectively the deviance and variance explained of a series of models, taking
 into account any transformations of the y variable. Functions \code{plotclean},
 \code{velout}, \code{codeplot} and \code{zapvelout} are useful to clean the
 data file.
}
\author{ Tim Cole \email{tim.cole@ucl.ac.uk} }
\references{
The idea of SITAR growth curve analysis arose from the paper by Beath (2007) and
was first described in Cole et al (2010). The other references describe
applications of SITAR to a variety of data forms.

Beath KJ. Infant growth modelling using a shape invariant model with random
efffects. Stat Med 2007;26:2547-64.

Cole TJ, Cortina Borja M, Sandhu J, et al. Nonlinear growth generates age
changes in the moments of the frequency distribution: the example of height in
puberty. Biostatistics 2008;9:159-71.

Cole TJ, Donaldson MD, Ben-Shlomo Y. SITAR--a useful instrument for growth curve
analysis. Int J Epidemiol 2010;39:1558-66.

Gault EJ, Perry RJ, Cole TJ, et al. Effect of oxandrolone and timing of pubertal
induction on final height in Turner's syndrome: randomised, double blind,
placebo controlled trial. BMJ 2011;342:d1980.

Johnson L, Llewellyn CH, van Jaarsveld CHM, et al. Genetic and environmental
influences on infant growth: prospective analysis of the Gemini twin birth
cohort. PLoS ONE 2011;6:e19918.

Prentice A, Dibba B, Sawo Y, et al. The effect of prepubertal calcium carbonate
supplementation on the age of peak height velocity in Gambian adolescents. Am J
Clin Nutr 2012;96:1042-50.

Dean MC, Cole TJ. Human life history evolution explains dissociation between the
timing of tooth eruption and peak rates of root growth. PLoS ONE 2013;8:e54534.
}
\keyword{ package }
