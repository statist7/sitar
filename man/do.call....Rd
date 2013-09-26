\name{do.call...}
\alias{do.call...}
\title{Extended version of do.call to handle repeated arguments}
\description{
The function \code{do.call} fails if arguments are duplicated in the call. \code{do.call...} includes \code{...} as an argument and checks for duplicate arguments. Those in \code{...} overwrite earlier arguments.}
\usage{
do.call...(what, args, ..., quote = FALSE, envir = parent.frame())
}
\arguments{
  \item{what}{
either a function or a non-empty character string naming the function to be called.}
  \item{args}{
a named list of arguments to the function call.}
  \item{\dots}{
a further list of arguments, some of which may also appear in args.}
  \item{quote}{
a logical value indicating whether to quote the arguments.}
  \item{envir}{
an environment within which to evaluate the call. This will be most useful if \code{what} is a character string and the arguments are symbols or quoted expressions.}
}
\value{The result of the (evaluated) function call.}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
