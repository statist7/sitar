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
The expression \code{fun} involves a single mention of a single variable named \code{varname} such that
\code{fun} is of the form {f(g(h(...(get(varname)))))}.

\code{ifun} then returns \code{h^-1(g^-1(f^-1(...(eval(fun)))))}.
Invertible arithmetic, logarithmic, exponential and power functions are recognised, but not trigonometric or hyperbolic.

When modelling with transformed variables, \code{ifun} enables the x and/or y variables to be automatically back-transformed to their original scales for plotting, simply by inverting the x and y transformations in the model call.
}
\value{
A length 2 list:
  \item{fn}{the inverse function with argument \code{x} which applied to \code{eval(fun)} returns \code{get(varname)}.}
  \item{varname}{the name of the variable extracted from \code{fun}.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\examples{
## simple case
ifun(quote(age))

## more complex case - transform age to log age since conception
fun <- quote(log(age + 0.75))

## create inverse function
(results <- ifun(fun))

##  apply inverse function to transformed age
age <- 1:10
(aget <- eval(fun))
(results$fn(aget))
}
