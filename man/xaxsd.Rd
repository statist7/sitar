\name{xaxsd}
\alias{xaxsd}
\alias{yaxsd}
\title{xaxsd and yaxsd}
\description{
Implements par('xaxs') and par('yaxs') option 'd'.
}
\usage{
xaxsd(usr = par()$usr[1:2])
yaxsd(usr = par()$usr[3:4])
}
\arguments{
  \item{usr}{a length-2 vector defining the length of the x-axis or y-axis.}
}
\details{
Implements par('xaxs') and par('yaxs') option 'd'
e.g. plot(..., xlim=xaxsd(), ylim=yaxsd()).
}
\value{
By default returns xlim/ylim args to match current setting of par()$usr, i.e. previous plot scales. Specifying \code{usr} gives scales with the usr args at the extremes.
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (usr = par()$usr[1:2]) 
(usr + mean(usr) * 0.08)/1.08
}
