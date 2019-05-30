#' Plot SITAR model
#'
#' \code{plot} and \code{lines} methods for objects of class \code{sitar},
#' providing various flavours of plot of the fitted growth curves. Also helper
#' functions to return the data for plotting, e.g. with \code{ggplot2}.
#'
#' For options involving both distance curves (options 'dcDua') and velocity curves
#' (options 'vV') the velocity curve plot (with right axis) can be annotated with
#' \code{par} parameters given as a named list called \code{y2par}.
#' To suppress the legend that comes with it set \code{legend = NULL}.
#'
#' The transformations \code{xfun} and \code{yfun} are applied to the x and y
#' variables after inverting any transformations applied in the original SITAR
#' call. So for example if \code{y = log(height)} in the SITAR call, then \code{yfun}
#' is applied to \code{height}. Thus the default \code{yfun = I} has the effect of
#' inverting the transformation. This is achieved by setting
#' \code{yfun = yfun(ifun(x$call.sitar$y))}.
#' For no transformation set \code{yfun = NULL}.
#'
#' The helper functions \code{plot_d}, \code{plot_v}, \code{plot_D},
#' \code{plot_V}, \code{plot_u}, \code{plot_a} and \code{plot_c}
#' correspond to the seven plot \code{option}s defined by their last letter,
#' and return the data for plotting as a \code{tibble}, e.g. for use with \code{ggplot2}.
#'
#' The \code{trim} option allows unsightly long line segments to be omitted
#' from plots with options 'a' or 'u'. It ranks the line segments on the basis
#' of the age gap (dx) and the distance of the midpoint of the line from the
#' mean curve (dy) using the formula \code{abs(dx)/mad(dx) + abs(dy)/mad(dy)}
#' and omits those with the largest values.

#' @aliases plot.sitar lines.sitar plot_d plot_v plot_D plot_V plot_u
#'  plot_a plot_c
#' @param x object of class \code{sitar}.
#' @param opt character string containing a subset of letters corresponding to
#' the options: 'd' for fitted Distance curve, 'v' for fitted Velocity curve,
#' 'c' for fitted Crosssectional distance curve, 'D' for individual fitted
#' Distance curves, 'V' for individual fitted Velocity curves, 'u' for
#' Unadjusted individual growth curves, and 'a' for Adjusted individual growth
#' curves. Options 'dvcDV' give spline curves, while 'ua' give data curves made
#' up as line segments. If both distance and velocity curves are specified, the
#' axis for the velocity curve appears on the right side of the plot (y2), and
#' a legend identifying the distance and velocity curves is provided.
#' @param labels optional character vector containing plot labels for \code{x},
#' \code{y} and \code{y} velocity from the original SITAR model. The three
#' elements can alternatively be provided via parameters
#' \code{xlab}, \code{ylab} and \code{vlab}. The latter take precedence.
#' Default labels are the names of \code{x} and \code{y}, and
#' "\code{y} velocity", suitably adjusted to reflect any back-transformation
#' via \code{xfun} and \code{yfun}.
#' @param apv optional logical specifying whether or not to calculate the age
#' at peak velocity from the velocity curve. If TRUE, age at peak velocity is
#' calculated as the age when the second derivative of the fitted curve changes
#' sign (after applying \code{xfun} and/or \code{yfun}). Age at peak velocity
#' is marked in the plot with a vertical dotted line, and its value, along with
#' peak velocity, is printed and returned. NB their standard errors can be
#' obtained using the bootstrap with the function \code{apv_se}.
#' @param xfun optional function to be applied to the x variable prior to
#' plotting (default I, see Details).
#' @param yfun optional function to be applied to the y variable prior to
#' plotting (default I, see Details).
#' @param subset optional logical vector of length \code{x} defining a subset
#' of \code{data} rows to be plotted, for \code{x} and \code{data} in the
#' original \code{sitar} call.
#' @param ns scalar defining the number of points for spline curves
#' (default 101).
#' @param abc vector of named values of random effects a, b and c used to
#' define an individual growth curve, e.g. abc=c(a=1, c=-0.1). Alternatively a
#' single character string defining an \code{id} level whose random effect
#' values are used. If \code{abc} is set, \code{level} is ignored. If
#' \code{abc} is NULL (default), or if a, b or c values are missing, values of
#' zero are assumed.
#' @param trim number (default 0) of long line segments to be excluded from plot
#' with option 'u' or 'a'. See Details.
#' @param add optional logical defining if the plot is pre-existing (TRUE) or
#' new (FALSE). TRUE is equivalent to using \code{lines}.
#' @param nlme optional logical which set TRUE plots the model as an
#' \code{nlme} object, using \code{plot.nlme} arguments.
#' @param returndata logical defining whether to plot the data (default FALSE)
#' or just return the data for plotting (TRUE). See Value.
#' @param \dots Further graphical parameters (see \code{par}) may also be
#' supplied as arguments, e.g. line
#' type \code{lty}, line width \code{lwd}, and colour \code{col}. For the
#' velocity (y2) plot \code{y2par} can be used (see Details).
#' @param xlab optional label for x axis
#' @param ylab optional label for y axis
#' @param vlab optional label for v axis (velocity)
#' @param xlim optional x axis limits
#' @param ylim optional y axis limits
#' @param vlim optional v axis limits
#' @param legend optional list of arguments for legend with distance-velocity plots
#' @return If \code{returndata} is FALSE returns invisibly a list of (up to) three objects:
#' \item{usr}{value of \code{par('usr')} for the main plot.}
#' \item{usr2}{the value of \code{par('usr')} for the velocity (y2) plot.}
#' \item{apv}{if argument \code{apv} is TRUE a named list giving the age at
#' peak velocity (apv) and peak velocity (pv) from the fitted velocity curve,
#' either overall or (with options D or V, invisibly) for all subjects.}
#' If \code{returndata} is TRUE (which it is with the helper functions) returns
#' invisibly either a tibble or named list of tibbles,
#' containing the data to be plotted. The helper functions each return a tibble.
#' The variable names are '.x', '.y' and
#' (for curves grouped by subject) '.id'. Note that '.x' and '.y' are returned
#' after applying \code{xfun} and \code{yfun}. Hence if for example \code{x = log(age)}
#' in the SITAR call then '.x' corresponds by default to \code{age}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{mplot}},
#' \code{\link{plotclean}}, \code{\link{ifun}}, \code{\link{apv_se}}
#' @keywords aplot
#' @examples
#'
#' ## fit sitar model
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=4)
#'
#' ## draw fitted distance and velocity curves
#' ## with velocity curve in blue
#' ## adding age at peak velocity (apv)
#' plot(m1, y2par=list(col='blue'), apv=TRUE)
#'
#' ## bootstrap standard errors for apv and pv
#' \dontrun{
#' res <- apv_se(m1, nboot=20, plot=TRUE)
#' }

#' ## draw individually coloured growth curves adjusted for random effects
#' ## using same x-axis limits as for previous plot
#' plot(m1, opt='a', col=id, xlim=xaxsd())
#'
#' ## add mean curve in red
#' lines(m1, opt='d', col='red', lwd=2)
#'
#' ## add mean curve for a, b, c = -1 SD
#' lines(m1, opt='d', lwd=2, abc=-sqrt(diag(getVarCov(m1))))
#'
#' ## draw fitted height distance curves coloured by subject, using ggplot
#' \dontrun{
#' require(ggplot2)
#' ggplot(plot_D(m1), aes(.x, .y, colour=.id)) +
#' labs(x='age', y='height') +
#' geom_line(show.legend=FALSE)
#' }

#' @importFrom grDevices xy.coords
#' @importFrom graphics plot axis identify legend lines locator par text title mtext abline
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate rename filter
#' @importFrom rlang .data as_label
#' @importFrom glue glue
#' @export
plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun=I, yfun=I, subset=NULL,
                       ns=101, abc=NULL, trim=0, add=FALSE, nlme=FALSE,
                       returndata=FALSE, ...,
                       xlab=NULL, ylab=NULL, vlab=NULL,
                       xlim=c(NA, NA), ylim=c(NA, NA), vlim=c(NA, NA),
                       legend=list(x='topleft', inset=0.04, bty='o')) {

  plotaxes <- function(dv, xlab, ylab, vlab, xlim, ylim, vlim, ...)
#	dv = 1 for distance, 2 for velocity, 3 for distance and velocity
#	returns par()$usr for d and v
  {
#	args for y1 and y2
    if (dv == 2) {
      ylab <- vlab
      ylim <- vlim
    } else if (dv == 3) {
      mar <- par()$mar
      mar[4] <- mar[2]
      par(mar=mar)
    }
    dots <- match.call(expand.dots=FALSE)$...
    ARG <- if (!is.null(dots))
      lapply(as.list(dots), eval, data, parent.frame())
    ARG1 <- ARG[names(ARG) != 'y2par']
    ARG2 <- ARG[['y2par']]
    if ('las' %in% names(ARG1))
      ARG2$las <- ARG1$las

#	plot x & y1 axes
    do.call('plot', c(list(x=xlim, y=ylim, type='n', xlab=xlab, ylab=ylab), ARG1))
    #	save x & y1 axis limits
    xy <- list()
    xy$usr <- par('usr')
#	optionally add right axis
    if (dv == 3) {
      par(new=TRUE)
      plot(xlim, vlim, type='n', bty='n', ann=FALSE, axes=FALSE)
      localaxis <- function(..., col, bg, pch, cex, lty, lwd) axis(...)
      do.call('localaxis', c(list(side=4), ARG2))
      mtext(vlab, 4, par('mgp')[1])
#	save y2 axis limits
      xy$usr2 <- par('usr')
      eval(parse(text=".par.usr2 <<- par('usr')"))
#	reset axis limits
      par(usr=xy$usr)
    } else {
# save null y2 axis limits
      eval(parse(text=".par.usr2 <<- NULL"))
    }
    invisible(xy)
  }

  v2d <- function(ylim, vlim) {
# converts velocity to distance
    rc <- lm(ylim ~ vlim)$coef
    function(v) rc[[1]] + rc[[2]] * v
  }

  xseq <- function(x, n=ns) {
# n is the number of points across the x range
    rx <- range(x, na.rm=TRUE)
    seq(rx[1], rx[2], length.out=n)
  }

  distance <- velocity <- function(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, ns=ns) {
# generate x values across the range to plot mean spline curve
    dvt <- as_label(match.call()[[1]])
    dvt <- as.numeric(dvt == 'velocity')
    .x <- getCovariate(model)[subset]
    .x <- xseq(.x, ns)
    newdata <- tibble(.x)
    if (sum(subset) < length(subset)) attr(newdata, 'subset') <- subset
    level <- ifelse(is.null(abc), 0, 1)
    . <- tibble(
      .x=xfun(.x),
      .y=predict(model, newdata, level=level, deriv=dvt, abc=abc, xfun=xfun, yfun=yfun)
    )
  }

  Distance <- Velocity <- function(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, ns=ns) {
# generate x and id values across the range to plot spline curves
    dvt <- as_label(match.call()[[1]])
    dvt <- as.numeric(dvt == 'Velocity')
    .x <- getCovariate(model)[subset]
    .id <- getGroups(model)[subset]
    npt <- ns / diff(range(.x))
    if (is.null(abc)) {
      . <- by(tibble(.x, .id), .id, function(z) {
        xrange <- range(z$.x)
        nt <- ceiling(npt * diff(xrange))
        tibble(
          .x=xseq(xrange, nt),
          .id=z$.id[[1]]
        )
      })
      . <- do.call('rbind', .)
    } else {
      . <- tibble(
        .x=xseq(.x, ns),
        .id=.id[[1]],
      )
    }
    . <- mutate(.,
                .y=predict(model, ., deriv=dvt, abc=abc, xfun=xfun, yfun=yfun),
                .x=xfun(.x)
    )
    .[, c('.x', '.y', '.id')]
  }

  unadjusted <- function(model, subset=subset, xfun=xfun, yfun=yfun, trim=trim) {
# unadjusted individual curves
    data <- tibble(
      .x=getCovariate(model),
      .y=getResponse(model),
      .id=getGroups(model)) %>%
      filter(subset)
    data <- trimlines(model, data, level=1, trim) %>%
      mutate(.x=xfun(.data$.x),
             .y=yfun(.data$.y))
    data
  }

  adjusted <- function(model, subset=subset, xfun=xfun, yfun=yfun, trim=trim) {
# adjusted individual curves
    data <- as_tibble(xyadj(model)) %>%
      mutate(.id=getGroups(model)) %>%
      rename(.x=.data$x,
             .y=.data$y) %>%
      filter(subset)
    data <- trimlines(model, data, level=0, trim) %>%
      mutate(.x=xfun(.data$.x),
             .y=yfun(.data$.y))
    data
  }

  trimlines <- function(model, data, level, trim) {
# if midpoint of line segment is far from mean curve, or age gap is large,
# insert NA row to data to omit line segment
    if (trim == 0)
      return(data)
    data <- with(data, data[order(.id, .x), ]) # sort data
    extra <- as_tibble(diff(as.matrix(data[, 1:2]))) # diff for .x and .y
    extra$.id <- data$.id[-1] # save .id
    did <- diff(as.integer(data$.id)) # 1+ = change of .id, 0 = .id
    extra$dx <- extra$.x # save dx age gap
    extra[, 1:2] <- data[-1, 1:2] - extra[, 1:2] / 2 # midpoints for .x and .y
    extra <- extra[!did, ] # restrict to same .id
    extra$ey <- predict(model, extra, level=level) # predicted y value at midpoint
    extra <- extra %>%
      mutate(dy=abs(.data$.y - .data$ey), # gap between line segment and mean curve dy
             xy=.data$dx / mad(.data$dx) + .data$dy / mad(.data$dy)) # add scaled dx and dy
    outliers <- order(extra$xy, decreasing=TRUE)[1:trim] # identify outliers
    extra <- extra[outliers, 1:3] # trim
    extra$.y <- NA
    data <- rbind(data, extra)
    with(data, data[order(.id, .x), ]) # sort data
  }

  crosssectional <- function(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, ns=ns) {
# fixed effect mean curve
    x <- getCovariate(model)[subset]
    x <- xseq(x, ns) - model$xoffset
    . <- tibble(
      .y=yfun(predict(model$ns, tibble(x))),
      .x=xfun(x + model$xoffset)
    )
    .[, c('.x', '.y')]
  }

  dolegend <- function(ARG1, ARG2, legend) {
# add legend
    parlu <- function(...) par(do.call('par', list(...)))
    if (is.null(ARG2$lty))
      ARG2$lty <- 2
    llc <- sapply(c('lty', 'lwd', 'col'), function(i) {
      p12 <- rep(par()[[i]], 2)
      for (j in 1:2) {
        k <- list(ARG1[i], ARG2[i])[[j]]
        if (!is.null(k[[1]]) && length(k[[1]]) == 1)
          p12[j] <- parlu(k)[[1]]
      }
      p12
    })
    do.call('legend', c(legend, list(lty=llc[, 'lty'], lwd=llc[, 'lwd'], col=llc[, 'col'])))
  }

  getlab <- function(lab, label, call, fun) {
    if (!is.null(lab))
      return(lab)
    if (label != '' && !is.na(label) && !is.null(label))
      return(label)
    if (fun == 'velocity')
      return(paste(call, 'velocity'))
    if (fun == 'NULL')
      return(deparse(call))
    icall <- ifun(call)
    lab <- attr(icall, 'varname')
    if (fun == 'I')
      return(lab)
    labfun <- ifelse(grepl('\\(', fun), paste0('(', fun, ')'), fun)
    paste0(labfun, '(', lab, ')')
  }

  getfun <- function(fun, call) {
    if (is.null(fun))
      return(function(x) x)
    function(x) fun(ifun(call)(x))
  }

# main code ---------------------------------------------------------------

  if (nlme)
    do.call('plot.lme', as.list(match.call()[-1]))
  model <- x
  data <- getData(model)
# check data versus model
  if (nrow(data) != model$dims$N)
    stop(glue('lengths of data ({nrow(data)}) and model ({model$dims$N}) do not match'))
  mcall <- model$call.sitar
#	extract list(...)
  ccall <- match.call(expand.dots=FALSE)
#	subset to plot model
  subset <- eval(ccall$subset, data, parent.frame())
  if (is.null(subset))
    subset <- rep_len(TRUE, model$dims$N)
# ... args
  dots <- ccall$...
  ARG <- if (!is.null(dots))
    lapply(as.list(dots), eval, data[subset, ], parent.frame())

  options   <- c('d', 'c', 'u', 'a', 'D', 'v', 'V')
  optnames  <- c('distance', 'crosssectional', 'unadjusted', 'adjusted', 'Distance', 'velocity', 'Velocity')
  optaxis   <- c( 1,   1,   1,   1,   1,   2,   2 ) # default y1=1, y2=2
  optmult   <- c( FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, TRUE ) # multiple curves
  optsmooth <- c( TRUE,  TRUE,  FALSE, FALSE, TRUE,  TRUE,  TRUE ) # spline curves
  opts <- unique(na.omit(match(unlist(strsplit(opt, '')), options)))
  if (length(opts) == 0)
    stop('option(s) not recognised')
  dv <- range(optaxis[opts])
  dv <- min(dv) + diff(dv) * 2 # 1=d, 2=v, 3=dv

# create missing labels
  if (missing(labels))
    labels <- vector('character', 3)

#	get axis labels
  xlab <- getlab(xlab, labels[1], mcall$x, paste(deparse(substitute(xfun)), collapse=''))
  ylab <- getlab(ylab, labels[2], mcall$y, paste(deparse(substitute(yfun)), collapse=''))
  vlab <- getlab(vlab, labels[3], ylab, 'velocity')

# get xfun and yfun
  xfun <- getfun(xfun, mcall$x)
  yfun <- getfun(yfun, mcall$y)

# generate list of data frames for selected options
  data <- lapply(opts, function(i) {
    if (optsmooth[[i]])
      do.call(optnames[[i]], list(model=model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, ns=ns))
    else
      do.call(optnames[[i]], list(model=model, subset=subset, xfun=xfun, yfun=yfun, trim=trim))
  })

# return data?
  if (returndata){
    if (length(data) == 1)
      data <- data[[1]]
    else
      names(data) <- options[opts]
    return(invisible(data))
  }

# extract axis ranges and plot axes
  if (!add) {
    if (any(is.na(xlim)))
      xlim <- range(vapply(data, function(z) range(z$.x, na.rm=TRUE), numeric(2)))
    if (any(is.na(ylim)) && any(optaxis[opts] == 1))
      ylim <- range(vapply(data[optaxis[opts] == 1], function(z) range(z$.y, na.rm=TRUE), numeric(2)))
    if (any(is.na(vlim)) && any(optaxis[opts] == 2))
      vlim <- range(vapply(data[optaxis[opts] == 2], function(z) range(z$.y, na.rm=TRUE), numeric(2)))
    xy <- do.call('plotaxes',
                  c(list(dv=dv, xlab=xlab, ylab=ylab, vlab=vlab,
                          xlim=xlim, ylim=ylim, vlim=vlim), ARG))
# add legend
    if (dv == 3 && !is.null(legend)) {
      legend[['legend']] <- c(ylab, vlab)
      dolegend(ARG[names(ARG) != 'y2par'], ARG$y2par, legend)
    }
  }
# else retrieve axis ranges
  else {
    xy <- list()
    xy$usr <- par('usr')
    xlim <- xaxsd()
    ylim <- yaxsd()
    if (exists('.par.usr2') && !is.null(.par.usr2)) {
      dv <- 3
      vlim <- yaxsd(.par.usr2[3:4])
      xy$usr2 <- .par.usr2
    } else if (dv == 3)
        stop('right y axis not set up')
  }

# plot curves
  lapply(1:length(opts), function(i) {
    opt <- opts[[i]]
    . <- data[[i]]
    ARG0 <- ARG
# data frame extended, extend ARG
    if (optmult[opt] && nrow(.) != model$dims$N && !is.null(dots)) {
      names(.) <- unlist(as.list(mcall[2:(length(.) + 1)]))
      ARG0 <- lapply(as.list(dots), eval, ., parent.frame())
    }
# select distance or velocity axis
    if (optaxis[opt] == 1 || dv < 3) {
      fun <- I
      ARG0 <- ARG0[names(ARG0) != 'y2par']
    } else {
      fun <- v2d(ylim, vlim)
      ARG0 <- ARG0[['y2par']]
      if (is.null(ARG0[['lty']]))
        ARG0[['lty']] <- 2
    }
    if (optmult[opt])
      do.call("mplot", c(list(x=.[[1]], y=fun(.[[2]]), id=.[[3]], add=TRUE), ARG0))
    else
      xy <- do.call("lines", c(list(x=.[[1]], y=fun(.[[2]])), ARG0))
  })

# save and print vertical line(s) at age of peak velocity
  if (apv) {
# single curve
      xy$apv <- with(velocity(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, ns=ns),
                     setNames(getPeak(.x, .y), c('apv', 'pv')))
    print(signif(xy$apv, 4))
    if (any(optmult[opts])) {
# multiple curves
      ids <- levels(factor(getGroups(model)[subset]))
      if (body(xfun) == as.name('x') &&
          body(yfun) == as.name('x')) {
# x and y untransformed
        apv1 <- xyadj(model, xy$apv[1], y=0, id=ids, tomean=FALSE)$x # subject-specific APVs
        pv1 <- xy$apv[2] * exp(ranef(model)$c) # subject-specific PVs
        xy$apv <- data.frame(apv=apv1, pv=pv1)
      }
# xfun or yfun set
      else {
        xt <- xseq(getCovariate(model)[subset])
        newdata <- data.frame(.x=xt)
        if (!is.null(abc))
          ids <- factor(1, labels=ids[[1]])
        . <- vapply(ids, function(z) {
          newdata$.id <- z
          vt <- predict(object=model, newdata=newdata, deriv=1, abc=abc, xfun=xfun, yfun=yfun)
          getPeak(xfun(xt), vt)
        }, numeric(2))
        xy$apv <- setNames(data.frame(t(.)), c('apv', 'pv'))
      }
    }
# plot apv
    do.call('abline', list(v=unlist(xy$apv['apv']), lty=3))
  }
# return xy
  invisible(xy)

}
#############################
#
#	lines.sitar
#
#############################

#' @rdname plot.sitar
#' @export
lines.sitar <- function (x, ...)
{
  mcall <- match.call()
  mcall[[1]] <- as.name("plot")
  mcall[['add']] <- TRUE
  eval(mcall, parent.frame())
}

#' @rdname plot.sitar
#' @export
plot_d <- function (x, ...)
{
  mcall <- match.call()
  opt <- deparse(mcall[[1]])
  opt <- substring(opt, nchar(opt), nchar(opt))
  mcall[[1]] <- as.name("plot")
  mcall[["opt"]] <- opt
  mcall[['returndata']] <- TRUE
  eval(mcall, parent.frame())
}

#' @rdname plot.sitar
#' @export
plot_v <- function (x, ...) {}

#' @rdname plot.sitar
#' @export
plot_D <- function (x, ...) {}

#' @rdname plot.sitar
#' @export
plot_V <- function (x, ...) {}

#' @rdname plot.sitar
#' @export
plot_u <- function (x, ...) {}

#' @rdname plot.sitar
#' @export
plot_a <- function (x, ...) {}

#' @rdname plot.sitar
#' @export
plot_c <- function (x, ...) {}

body(plot_v) <- body(plot_D) <- body(plot_V) <- body(plot_u) <- body(plot_a) <- body(plot_c) <- body(plot_d)
