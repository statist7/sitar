\name{zapvelout}
\alias{zapvelout}
\title{zapvelout}
\description{
\code{zapvelout} sets to missing values of \code{y} flagged as outliers by \code{velout}.}
\usage{
zapvelout(outliers, icode)
}
\arguments{
  \item{outliers}{dataframe returned by \code{\link{velout}}.}
  \item{icode}{code(s) to be zapped. See Details.}
}
\details{
Codes range from 0 (normal) to 8, where 4 and 6 are conventional outliers (see \code{\link{velout}}.}
\value{
Returns dataframe \code{data} from \code{velout} with values of \code{y} set to \code{NA} corresponding to the specified \code{icode}(s).}
\author{Tim Cole \email{tim.cole@ucl.ac.uk}}
\seealso{\code{\link{codeplot}}, \code{\link{velout}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (outliers, icode) 
{
    data <- get(attr(outliers, "data"))
    if (dim(data)[1] != dim(outliers)[1]) 
        stop("numbers of rows in data and outliers differ")
    y <- names(outliers)[3]
    zap <- outliers$code \%in\% icode
    cat(sum(zap), y, "values set missing\n")
    data[outliers$code \%in\% icode, y] <- NA
    invisible(data)
  }
}
