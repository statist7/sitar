#' Select equally spaced ages from a vector of ages
#'
#' \code{timegap} indexes elements in a vector of ages such that the indexed
#' ages are spaced integer multiples of a time interval apart, to within a given
#' tolerance. \code{timegap.id} is a wrapper to apply \code{timegap} within levels
#' of factor \code{id}. The selected ages can then be split into age groups the
#' specified interval wide, ensuring that each subject has just one measurement
#' per interval.
#'
#' \code{timegap} calculates all possible differences between pairs of ages,
#' expresses them as integer multiples of \code{gap}, restricts them to
#' those within tolerance and identifies those providing the longest sequence.
#' For solutions of the same length, those with the smallest standard deviation
#' of successive differences are selected.
#'
#' With \code{timegap}, for unique solutions, or multiple solutions with
#' \code{multiple FALSE},
#' a vector of indices the same length as, and named with, \code{age}.
#' Where there are multiple solutions with \code{multiple TRUE},
#' they are returned as a named matrix. With \code{timegap.id} the subject
#' vectors are returned invisibly, concatenated into a single vector.
#'
#' @aliases timegap timegap.id diff.id
#' @param age vector of ages.
#' @param id factor of subject ids.
#' @param data frame optionally containing \code{age} and \code{id}.
#' @param gap numeric, the required time gap between selected ages.
#' @param tol numeric, the tolerance around the gap (default \code{0.1 * gap}).
#' @param multiple logical, whether or not to return multiple solutions
#' when found (default FALSE).
#' @return index of selected ages, named by value of \code{age}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @examples
#' data(heights)
#' ## select heights measured multiples of 1 year apart
#' (timegap.id(age, id, heights, 1))
#'
#' ## select heights measured multiples of 4 months apart
#' (timegap.id(age, id, heights, 1/3))
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
  resmat[abs(resmat) < .Machine$double.eps * 10] <- 0
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

#' @rdname timegap.id
#' @export
timegap.id <- function(age, id, data=parent.frame(), gap, tol=0.1*gap, multiple=FALSE) {
  bylist <- with(data, by(data, id, function(z)
    list(timegap(z$age, gap=gap, tol=tol, multiple=multiple), nrow(z))
  ))
  count <- 0
  subset <- integer(0)
  for (i in 1:length(bylist)) {
    subset <- c(subset, bylist[[i]][[1]] + count)
    count <- count + bylist[[i]][[2]]
  }
  invisible(subset)
}

#' @rdname diff.id
#' @export
diff.id <- function(age, id, data=parent.frame()) {
  on.exit(detach(data))
  eval(parse(text='attach(data)'))
  age <- setNames(age, id)
  xd <- diff(age)
  xd[diff(as.numeric(id)) != 0] <- NA
  c(xd, NA)
}
