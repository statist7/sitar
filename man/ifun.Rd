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

When modelling with transformed variables, \code{ifun} enables the \code{x} and/or
\code{y} variables to be automatically back-transformed to their original scales
for plotting, simply by inverting the \code{x} and \code{y} transformations in the
model call. \code{\link{plot.sitar}} uses \code{ifun} for this purpose, via the
functions \code{xfun} and \code{yfun}.
}
\value{
A length 2 list:
  \item{fn}{the inverse function with argument \code{x} which applied to \code{eval(fun)} returns \code{get(varname)}.}
  \item{varname}{the name of the variable extracted from \code{fun}.}
}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{
\code{\link{plot.sitar}}
}
\examples{
## simplest case
(ifun(quote(age)))

(age <- 1:10)

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
(fun <- quote((sqrt(sin(log1p(age+0.75)/2)+1)^4)))
(ifun(fun))
(transformed.age <- eval(fun))

## ta da!
(ifun(fun)$fn(transformed.age))
}
