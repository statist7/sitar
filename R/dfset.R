#' Find degrees of freedom for a natural spline curve to minimise BIC or AIC
#'
#' \code{dfset} fits a natural cubic spline for a range
#' of degrees of freedom, and returns the df minimising the BIC or AIC.
#'
#' @param x vector of x coordinates.
#' @param y vector of y coordinates.
#' @param data data frame containing \code{x} and \code{y}.
#' @param FUN function to be minimised (e.g. BIC or AIC).
#' @param df vector of degrees of freedom to be searched.
#' @param plot logical controlling plotting of FUN versus df.
#' @param \dots parameters to pass to \code{plot}.
#' @return Optimal degrees of freedom.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' dfset(age, height, heights, FUN=BIC, plot=TRUE)
#' dfset(age, height, heights, FUN=function(a) AIC(a, k=1))
#' @export dfset
  dfset <- function(x, y, data=parent.frame(), FUN=BIC, df=1:15, plot=FALSE, ...) {
    mc <- match.call()
    x <- eval(mc$x, data)
    y <- eval(mc$y, data)
    . <- vapply(df, function(i) {
      obj <- try(lm(y ~ ns(x, df=i)), silent=TRUE)
      if (inherits(obj, 'lm'))
        FUN(obj)
      else
        NA
    }, 0)
    if (plot) {
      plot(. ~ df, ..., ylab=deparse(substitute(FUN)))
      abline(h=df[which.min(.)], lty=3)
    }
    df[which.min(.)]
  }
