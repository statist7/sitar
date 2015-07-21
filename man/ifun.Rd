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
The expression \code{fun} involves a single variable named \code{varname} such that
\code{fun = f(g(h(...(get(varname)))))}. \code{ifun} returns \code{h^-1(g^-1(f^-1(...(fun))))}.
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

## more complex case - define age transformation
fun <- quote(log(age + 0.75))

## create inverse function
(results <- ifun(fun))

##  apply inverse function to transformed age
age <- 1:10
(results$fn(eval(fun)))
}
