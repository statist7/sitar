% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDV.R
\name{getDV}
\alias{getDV}
\title{Extract individual estimates of peak velocity and age at peak velocity}
\usage{
getDV(object, x = "apv")
}
\arguments{
\item{object}{a SITAR model.}

\item{x}{the age of interest, specified either as 'apv', the mean age at peak
velocity, or 'ato', the mean age at take-off, or a numerical age.}
}
\value{
A tibble with one row per individual and five columns, where id, age
  and distance have names corresponding to those used in the SITAR call:
  \item{id}{subject id.} \item{age}{subject's age corresponding to \code{x}.}
  \item{distance}{subject's distance at specified age.}
  \item{velocity}{subject's velocity at specified age.}
  \item{missing}{logical where TRUE means subject's specified age lies outside
  their measurement range.}
}
\description{
\code{getDV} takes a point on the mean fitted velocity curve, defined
typically as the age at peak velocity or an equivalent landmark age, and maps
it onto the corresponding point on individual velocity curves taking into
account their timing and intensity.
}
\details{
SITAR is a shape-invariant model, so if there is a turning point on the mean
velocity curve, e.g. a peak in puberty, then there will be a corresponding
turning point on the velocity curves of all individuals, at an age depending
on their timing and intensity. This applies even to individuals who lack
measurements at that age, have stopped earlier or started later, so their
growth curve is incomplete and lacks the turning point. The returned variable
\code{missing} flags such individuals.

Note that 'D' and 'V' in \code{getDV} correspond to the plot options for
individual Distance and Velocity curves.
}
\examples{

data(heights)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

# fit sitar model
model <- sitar(x = log(age), y = height, id = id, data = heights, df = 5)
(dv <- getDV(model))
(id_missing <- dv \%>\%
  filter(missing) \%>\%
  pull(id))

# plot individual velocity curves
ggplot(plot_V(model), aes(age, height, group = id, colour = id)) +
    theme(legend.position = 'inside',
          legend.position.inside = c(0.9, 0.6)) +
    geom_line() +

# add individual peak velocities
geom_point(aes(y = velocity), data = dv) +

# highlight subjects 2, 7 and 12 with dashed lines
# despite incomplete curves their peak velocities are estimated
geom_line(data = . \%>\% filter(id \%in\% id_missing), linetype = 2, colour = 'white') +
            geom_point(aes(y = velocity), data = dv \%>\%
                         filter(id \%in\% id_missing), shape = 1, size = 3)

}
\author{
Tim Cole \email{tim.cole@ucl.ac.uk}
}
