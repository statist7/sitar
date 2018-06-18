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
#' The helper functions \code{plot_d}, \code{plot_v}, \code{plot_D},
#' \code{plot_V}, \code{plot_u}, \code{plot_a} and \code{plot_c}
#' correspond to the seven plot \code{option}s defined by their last letter,
#' and return the data for plotting, e.g. for use with \code{ggplot2}.
#'
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
#' peak velocity, is printed and returned.
#' @param xfun optional function to be applied to the x variable prior to
#' plotting. Defaults to NULL, which translates to \code{ifun(x$call.sitar$x)}
#' and inverts any transformation applied to x in the original SITAR model
#' call. To plot on the transformed scale set \code{xfun} to \code{I}.
#' @param yfun optional function to be applied to the y variable prior to
#' plotting. Defaults to NULL, which translates to \code{ifun(x$call.sitar$y)}
#' and inverts any transformation applied to y in the original SITAR model
#' call. To plot on the transformed scale set \code{yfun} to \code{I}.
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
#' @param add optional logical defining if the plot is pre-existing (TRUE) or
#' new (FALSE). TRUE is equivalent to using \code{lines}.
#' @param nlme optional logical which set TRUE plots the model as an
#' \code{nlme} object, using \code{plot.nlme} arguments.
#' @param returndata logical defining whether to plot the data (default FALSE)
#' or just return the data for plotting (TRUE).
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
#' either overall or (with options D or V) for all subjects.}
#' If \code{returndata} is TRUE (which it is with the helper functions) returns
#' invisibly either a tibble or named list of tibbles,
#' containing the data to be plotted. The helper functions each return a tibble.
#' The variable names are '.x', '.y' and
#' (for curves grouped by subject) '.id'. Note that '.x' and '.y' are returned
#' after applying \code{xfun} and \code{yfun}. Hence if for examplex \code{x = log(age)}
#' in the original \code{sitar} call then '.x' corresponds by default to \code{age}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{mplot}},
#' \code{\link{plotclean}}, \code{\link{ifun}}
#' @keywords aplot
#' @examples
#'
#' ## fit sitar model
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=4)
#'
#' ## draw fitted distance and velocity curves
#' ## with velocity curve in blue
#' ## adding age at peak velocity
#' plot(m1, y2par=list(col='blue'), apv=TRUE)
#'
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
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @export
plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun=NULL, yfun=NULL, subset=NULL,
                       ns=101, abc=NULL, add=FALSE, nlme=FALSE,
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
    if ('las' %in% names(ARG1)) ARG2$las <- ARG1$las

#	plot x & y1 axes
    do.call('plot', c(list(x=xlim, y=ylim, type='n', xlab=xlab, ylab=ylab), ARG1))
    #	save x & y1 axis limits
    xy <- list()
    xy$usr <- par('usr')
#	optionally add right axis
    if (dv == 3) {
      par(new=TRUE)
      plot(xlim, vlim, type='n', bty='n', ann=FALSE, axes=FALSE)
      do.call('axis', list(side=4))
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

  distance <- velocity <- function(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, n=ns) {
# generate x values across the range to plot mean spline curve
    dvt <- rlang::quo_name(match.call()[[1]])
    dvt <- as.numeric(dvt == 'velocity')
    .x <- getCovariate(model)[subset]
    .x <- xseq(.x, n)
    . <- tibble::tibble(
      .y=predict(model, tibble::tibble(.x), level=0, deriv=dvt, abc=abc, xfun=xfun, yfun=yfun),
      .x=xfun(.x)
    )
    if (length(.x[subset]) != length(.x)) attr(., 'subset') <- subset
    .[, c('.x', '.y')]
  }

  Distance <- Velocity <- function(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, n=ns) {
# generate x and id values across the range to plot spline curves
    dvt <- rlang::quo_name(match.call()[[1]])
    dvt <- as.numeric(dvt == 'Velocity')
    .x <- getCovariate(model)[subset]
    .id <- getGroups(model)[subset]
    npt <- n / diff(range(.x))
    . <- by(tibble::tibble(.x, .id), .id, function(z) {
      xrange <- range(z$.x)
      nt <- floor(npt * diff(xrange)) + 1
      tibble::tibble(
        .x=xseq(xrange, nt),
        .id=z$.id[[1]]
      )
    })
    . <- do.call('rbind', .)
    . <- dplyr::mutate(.,
                .y=predict(model, ., deriv=dvt, abc=abc, xfun=xfun, yfun=yfun),
                .x=xfun(.x)
    )
    .[, c('.x', '.y', '.id')]
  }

  unadjusted <- function(model, subset=subset, xfun=xfun, yfun=yfun) {
# unadjusted individual curves
    tibble::tibble(
      .x=xfun(getCovariate(model)[subset]),
      .y=yfun(getResponse(model)[subset]),
      .id=getGroups(model)[subset]
    )
  }

  adjusted <- function(model, subset=subset, xfun=xfun, yfun=yfun) {
# adjusted individual curves
    . <- xyadj(model)
    tibble::tibble(
      .x=xfun(.$x)[subset],
      .y=yfun(.$y)[subset],
      .id=getGroups(model)[subset]
    )
  }

  crosssectional <- function(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, n=ns) {
# fixed effect mean curve
    x <- getCovariate(model)[subset]
    x <- xseq(x, n) - model$xoffset
    . <- tibble::tibble(
      .y=yfun(predict(model$ns, tibble::tibble(x))),
      .x=xfun(x + model$xoffset)
    )
    .[, c('.x', '.y')]
  }

  dolegend <- function(ypar, y2par, legend) {
# add legend
    parlu <- function(...) par(do.call('par', list(...)))
    llc <- sapply(c('lty', 'lwd', 'col'), function(i) {
      p12 <- rep(par()[[i]], 2)
      if (!is.null(ypar[[i]]))
        p12[1] <- parlu(ypar[i])[[1]]
      if (!is.null(y2par[[i]]))
        p12[2] <- parlu(y2par[i])[[1]]
      p12
    })
    do.call('legend', c(legend, list(lty=llc[, 'lty'], lwd=llc[, 'lwd'], col=llc[, 'col'])))
  }

  if (nlme) {
    do.call('plot.lme', as.list(match.call()[-1]))
  }
  else {
    model <- x
    data <- getData(model)
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
    ARG1 <- ARG[names(ARG) != 'y2par']
    ARG2 <- ARG[['y2par']]

    options   <- c('d', 'c', 'u', 'a', 'D', 'v', 'V')
    optnames  <- c('distance', 'crosssectional', 'unadjusted', 'adjusted', 'Distance', 'velocity', 'Velocity')
    optaxis   <- c( 1,   1,   1,   1,   1,   2,   2 ) # default y1=1, y2=2
    optmult   <- c( FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, TRUE ) # multiple curves
    optsmooth <- c( TRUE,  TRUE,  FALSE, FALSE, TRUE,  TRUE,  TRUE ) # spline curves
    optDV     <- c( FALSE, FALSE, FALSE, FALSE, TRUE,  FALSE, TRUE ) # extended curves
    opts <- unique(na.omit(match(unlist(strsplit(opt, '')), options)))
    if (length(opts) == 0)
      stop('option(s) not recognised')
    dv <- range(optaxis[opts])
    dv <- min(dv) + diff(dv) * 2 # 1=d, 2=v, 3=dv

# create missing labels
    if (missing(labels))
      labels <- vector('character', 3)
    else if (length(labels) < 3)
      labels <- c(labels, '', '')
    if (labels[3] == '' && labels[2] != '')
      labels[3] <- paste(labels[2], 'velocity')

#	if xlab not specified replace with label or x name (depending on xfun)
    if (is.null(xlab)) {
      if (labels[1] == '') {
        xlab <- if (!is.null(xfun))
          paste0('(', deparse(substitute(xfun)), ')(', deparse(mcall$x), ")")
        else
          attr(ifun(mcall$x), 'varname')
        labels[1] <- xlab
      }
      else
        xlab <- labels[1]
    }

#	if ylab not specified replace with label or y name (depending on yfun)
    if (is.null(ylab)) {
      if (labels[2] == '') {
        ylab <- if (!is.null(yfun))
          paste0('(', deparse(substitute(yfun)), ')(', deparse(mcall$y), ")")
        else
          attr(ifun(mcall$y), 'varname')
        labels[2] <- ylab
      }
      else
        ylab <- labels[2]
    }

#	if vlab not specified replace with label or v name
    if (is.null(vlab)) {
      if (labels[3] == '')
        labels[3] <- paste(labels[2], 'velocity')
      vlab <- labels[3]
    }

# derive xfun and yfun
    if (is.null(xfun))
      xfun <- ifun(mcall$x)
    if (is.null(yfun))
      yfun <- ifun(mcall$y)

# generate list of data frames for selected options
    data <- lapply(opts, function(i) {
      if (optsmooth[[i]])
        do.call(optnames[[i]], list(model=model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, n=ns))
      else
        do.call(optnames[[i]], list(model=model, subset=subset, xfun=xfun, yfun=yfun))
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
        xlim <- range(vapply(data, function(z) range(z$.x), 0:1/2))
      if (any(is.na(ylim)) && any(optaxis[opts] == 1))
        ylim <- range(vapply(data[optaxis[opts] == 1], function(z) range(z$.y), 0:1/2))
      if (any(is.na(vlim)) && any(optaxis[opts] == 2))
        vlim <- range(vapply(data[optaxis[opts] == 2], function(z) range(z$.y), 0:1/2))
      xy <- do.call('plotaxes',
                    c(list(dv=dv, xlab=xlab, ylab=ylab, vlab=vlab,
                            xlim=xlim, ylim=ylim, vlim=vlim), ARG))
# add legend
      if (dv == 3) {
        # default dotted line for velocity curve
        if (is.null(ARG$y2par$lty)) {
          ARG$y2par$lty <- 2
          ARG2 <- ARG$y2par
        }
        if (!is.null(legend)) {
          legend[['legend']] <- labels[2:3]
          dolegend(ARG1, ARG2, legend)
        }
      }
    }
# else retrieve axis ranges
    else {
      xy <- list()
      xy$usr <- par('usr')
      xlim <- xaxsd()
      ylim <- yaxsd()
      vlim <- .par.usr2
      if (!is.null(vlim)) {
        xy$usr2 <- vlim
        vlim <- yaxsd(.par.usr2[3:4])
        dv <- 3
      } else if (dv == 3)
          stop('right y axis not set up')
    }

# plot curves
    lapply(1:length(opts), function(i) {
      opt <- opts[[i]]
      . <- data[[i]]
      ARG0 <- ARG
      # if D or V and dots, extend ARG
      if (optDV[opt] && !is.null(dots)) {
        names(.) <- unlist(as.list(mcall[2:(length(.)+1)]))
        ARG0 <- lapply(as.list(dots), eval, ., parent.frame())
      }
      # select distance or velocity axis
      if (optaxis[opt] == 1 || dv < 3) {
        fun <- I
        ARG0 <- ARG0[names(ARG0) != 'y2par']
      } else {
        fun <- v2d(ylim, vlim)
        ARG0 <- ARG0[['y2par']]
      }
      if (optmult[opt])
        do.call("mplot", c(list(x=.[[1]], y=fun(.[[2]]), id=.[[3]], add=TRUE), ARG0))
      else
        xy <- do.call("lines", c(list(x=.[[1]], y=fun(.[[2]])), ARG0))
    })

# save and print vertical line(s) at age of peak velocity
    if (apv) {
# single curve
        xy$apv <- with(velocity(model, subset=subset, abc=abc, xfun=xfun, yfun=yfun, n=ns),
                       setNames(getPeakTrough(.x, .y), c('apv', 'pv')))
      print(signif(xy$apv, 4))
      if (any(optmult[opts])) {
# multiple curves
        if (body(ifun(mcall$x)) == as.name('x') &&
            body(ifun(mcall$y)) == as.name('x')) {
# x and y untransformed
          apv1 <- xyadj(model, xy$apv[1], tomean=FALSE)$x # subject-specific APVs
          pv1 <- xy$apv[2] * exp(ranef(model)$c) # subject-specific PVs
          xy$apv <- data.frame(apv=apv1, pv=pv1)
        }
# xfun or yfun set
        else {
          xt <- xseq(getCovariate(model)[subset])
          newdata <- data.frame(.x=xt)
          . <- vapply(levels(factor(getGroups(model)[subset])), function(z) {
            newdata$.id <- z
            vt <- predict(object=model, newdata=newdata, deriv=1, abc=abc, xfun=xfun, yfun=yfun)
            getPeakTrough(xfun(xt), vt)
          }, c(0, 0))
          xy$apv <- setNames(data.frame(t(.)), c('apv', 'pv'))
        }
      }
# plot apv
      ARG2$lty <- 3
      do.call('abline', c(list(v=unlist(xy$apv['apv'])), ARG2))
    }
# return xy
    invisible(xy)
  }
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
