\name{funcall}
\alias{funcall}
\title{Function call with optional inverse}
\description{
Applies an expression to vector v, optionally inverting the expression first. For example if the expression is log, funcall returns log(v) if inverse is FALSE, and exp(v) if inverse is TRUE.}
\usage{funcall(v, vcall, inverse = FALSE)}
\arguments{
  \item{v}{vector}
  \item{vcall}{expression}
  \item{inverse}{logical}
}
\details{Inverse covers functions log, exp, sqrt, ^, *, /, +, -.}
\value{Returns a vector of length v.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
