#' Select equally spaced ages from a vector of ages
#'
#' \code{timegap} indexes elements in a vector of ages such that the indexed
#' ages are spaced integer multiples of a time interval apart, to within a given
#' tolerance. \code{timegap.id} is a wrapper to apply \code{timegap} within levels
#' of factor \code{id}. The selected ages can then be split into age groups the
#' specified time interval wide, ensuring that (virtually) every subject
#' has at most one measurement per interval.
#'
#' \code{timegap} calculates all possible differences between pairs of ages,
#' expresses them as integer multiples of \code{gap}, restricts them to
#' those within tolerance and identifies those providing the longest sequences.
#' For sequences of the same length, those with the smallest standard deviation
#' of successive differences (modulo the time interval) are selected.
#'
#' @aliases timegap timegap.id diffid
#' @param age vector of ages.
#' @param id factor of subject ids.
#' @param data data frame optionally containing \code{age} and \code{id}.
#' @param gap numeric, the required positive time gap between selected ages.
#' @param tol numeric, the positive tolerance around the gap (default \code{0.1 * gap}).
#' @param multiple logical, whether or not to return multiple solutions
#' when found (default FALSE).
#' @param lag an integer indicating which lag to use.
#' @param differences an integer indicating the order of the difference.
#' @param sort a logical indicating whether to first sort by id and age.
#' @param keepNA a logical indicating whether to keep generated NAs.
#' @return With \code{timegap}, for unique solutions, or multiple solutions with
#' \code{multiple FALSE}, a vector of indices named with \code{age}. With
#' \code{timegap.id} the subject vectors are returned invisibly, concatenated.
#'
#' With \code{multiple TRUE}, where there are multiple solutions
#' they are returned as a named matrix.
#'
#' \code{diffid} returns \code{diff(age)} applied within \code{id}.
#' With \code{keepNA} TRUE a suitable number of \code{NA}s are added at the end,
#' while if FALSE all \code{NA}s are omitted.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#'
#' ## bin age into 1-year groups by id
#' ## gives multiple measurements per id per year
#' with(heights, table(floor(age), id))
#'
#' ## now select heights measured multiples of 1 year apart
#' (tg1 <- timegap.id(age, id, heights, 1))
#'
#' ## no more than one measurement per id per year
#' with(heights[tg1, ], table(floor(age), id))
#'
#' ## most time intervals close to 1 year
#' summary(diffid(age, id, heights[tg1, ], lag=1))
#' @export timegap
timegap <- function(age, gap, tol=0.1*gap, multiple=FALSE) {
# age, a vector of ages
# gap, target time gap between measurements
# tol, tolerance criterion defaults to 0.1 * gap
# multiple, if TRUE and more than one optimal solution, return all (default FALSE)
# returns a vector of selected positions in age
  nt <- length(age)
  to <- order(age)
  age <- na.omit(age[to])
  to <- order((1:nt)[to])
  nt <- length(age)
  if (nt < 2 || gap <= 0 || tol <= 0) return(integer(0))
# symmetric matrix of differences
  tt <- matrix(abs(c(outer(age, age, '-'))), nt, nt, dimnames=list(age, age))
# res and int of matrix
  intmat <- tt %/% gap
  resmat <- tt - intmat * gap
# trap rounding error
  resmat[isTRUE(all.equal(resmat, 0))] <- 0
  intmat <- intmat * (resmat < tol)
# select best choices from many
  for (i in 1:nt) {
    intcol <- intmat[, i]
    rescol <- resmat[, i]
    for (j in unique(intcol)) {
      choose <- intcol == j
      if (sum(choose) > 0) {
        resmat[choose, i] <- 0
        resmat[choose, i][which.min(rescol[choose])] <- 1
      }
    }
  }
# compare length of solutions
  margin <- apply(resmat, 2, sum)
# if length only 1 return null
  if (max(margin) == 1) return(integer(0))
# compare longest solutions (at least 3 long to calculate sd)
  if (max(margin) > 2 & sum(margin == max(margin)) > 1) {
# penalise with sd of differences
    margin <- margin - apply(resmat, 2, function(z) {
      tt <- (age[as.logical(z)])
      tt <- tt[order(tt)]
      sd(abs((diff(tt) + gap/2) %% gap - gap/2))
    })
    margin[is.na(margin)] <- 0
# if still multiple longest solutions, simplify matrix
    if (sum(margin == max(margin)) > 1 & multiple) {
      resmat <- resmat[, margin==max(margin), drop=FALSE]
      resmat <- resmat[apply(resmat, 1, sum) > 0, , drop=FALSE]
      resmat <- unique(resmat, MARGIN=2)
      nt <- ncol(resmat)
      if (nt > 1) {
        dimnames(resmat)[[2]] <- 1:nt
        return(t(resmat))
      }
      margin <- 1
    }
  }
# return first optimal solution
  return(which((resmat[, which.max(margin)] == 1)[to]))
}

#' @rdname timegap
#' @export
timegap.id <- function(age, id, data=parent.frame(), gap, tol=0.1*gap, multiple=FALSE) {
  xd <- data.frame(lapply(match.call()[3:2], eval, data))
  xd$id <- factor(xd$id)
  bylist <- with(xd, by(xd, id, function(z)
    list(timegap(z$age, gap=gap, tol=tol, multiple=multiple), nrow(z))
  ))
  count <- 0
  subset <- integer()
  for (i in 1:length(bylist)) {
    subset <- c(subset, bylist[[i]][[1]] + count)
    count <- count + bylist[[i]][[2]]
  }
  invisible(subset)
}

#' @rdname timegap
#' @export
diffid <- function(age, id, data=parent.frame(), lag=1, differences=1, sort=FALSE, keepNA=FALSE) {
  mc <- match.call()
  xd <- cbind(as.integer(factor(eval(mc$id, data))), eval(mc$age, data))
  if (sort) xd <- xd[order(xd[, 1], xd[, 2]), ]
  xd <- diff(xd, lag=lag, differences=differences)
  xd[xd[, 1] != 0, 2] <- NA
  xd <- xd[, -1]
  if (keepNA) c(xd, rep(NA, lag * differences))
    else na.omit(xd)
}
