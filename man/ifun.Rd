\name{ifun}
\alias{ifun}
\title{Invert an expression defining a data transformation}
\description{
Enables a transformed variable to be back-transformed, e.g. for plotting purposes,
by inverting the transforming expression.
}
\usage{
ifun(fun)
}
\arguments{
  \item{fun}{a language object defining the expression to be inverted.}
}
\details{
The expression \code{fun} involves a single variable named \code{varname} such that \code{fun = f(g(h((get(varname)))))}, where each function \code{f(.)} involves up to one function of \code{varname} and up to one numerical constant (including \code{pi}). The number of such functions is in principle unlimited.

\code{ifun} then returns \code{h^-1(g^-1(f^-1((eval(fun)))))}, which ensures that the expression is invertible if and only if the individual functions are invertible.

When modelling with transformed variables, \code{ifun} enables the \code{x} and/or
\code{y} variables to be automatically back-transformed to their original scales
for plotting, simply by inverting the \code{x} and \code{y} transformations in the
model call. \code{\link{plot.sitar}} uses \code{ifun} for this purpose, via the
arguments \code{xfun} and \code{yfun}.
}
\value{
A list of length two:
  \item{fn}{the inverse function with argument \code{x} which applied to \code{eval(fun)} returns \code{get(varname)}.}
  \item{varname}{the name of the variable in \code{fun}.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{
\code{\link{plot.sitar}}
}
\examples{
## simplest case
(ifun(quote(age)))

(age <- 1:9)

## more complex case
## transform age to log age since conception
(fun <- quote(log(age + 0.75)))

## create inverse function
(ifun(fun))

##  transform age
(transformed.age <- eval(fun))

## ifun back-transforms transformed.age to age
(ifun(fun)$fn(transformed.age))

## ludicrously complex case
## including trigonometry etc
(fun <- quote((exp(sin(pi * log10(age + 0.75)/2) - 1)^4)))
(ifun(fun))
(transformed.age <- eval(fun))

## ta da!
(ifun(fun)$fn(transformed.age))
}
