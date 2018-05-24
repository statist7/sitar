#' Plot SITAR model
#'
#' \code{plot} and \code{lines} methods for objects of class \code{sitar},
#' providing various flavours of plot of the fitted growth curves.
#'
#' For option 'dv' (the default) the velocity curve plot (with right axis) can
#' be annotated with \code{par} parameters given as a named list called
#' \code{y2par}. To suppress the legend that comes with it set \code{xlegend =
#' NULL}.
#'
#' @aliases plot.sitar lines.sitar
#' @param x object of class \code{sitar}.
#' @param opt character string containing a subset of letters corresponding to
#' the options: 'd' for fitted Distance curve, 'v' for fitted Velocity curve,
#' 'e' for fitted fixed Effects distance curve, 'D' for individual fitted
#' Distance curves, 'V' for individual fitted Velocity curves, 'u' for
#' Unadjusted individual growth curves, and 'a' for Adjusted individual growth
#' curves. Options 'dveDV' give spline curves, while 'ua' give data curves made
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
#' of \code{data} rows to be plotted.
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
#' @param \dots Further graphical parameters (see \code{par}) may also be
#' supplied as arguments, e.g. axis labels \code{xlab} and \code{ylab}, line
#' type \code{lty}, line width \code{lwd}, and color \code{col}. For the
#' velocity (y2) plot \code{y2par} can be used (see Details).
#' @return Returns invisibly a list of three objects:
#' \item{ss}{\code{smooth.spline} object corresponding to the fitted distance
#' curve.} \item{usr}{value of \code{par('usr')} for the main plot.}
#' \item{usr2}{the value of \code{par('usr')} for the velocity (y2) plot.}
#' \item{apv}{if argument \code{apv} is TRUE a named list giving the age at
#' peak velocity (apv) and peak velocity (pv) from the fitted velocity curve.}
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{mplot}},
#' \code{\link{plotclean}}, \code{\link{y2plot}}, \code{\link{ifun}}
#' @keywords aplot
#' @examples
#'
#' ## fit sitar model
#' m1 <- sitar(x=age, y=height, id=id, data=heights, df=7)
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
#' @importFrom grDevices xy.coords
#' @importFrom graphics axis identify legend lines locator par text title mtext abline
#' @export
plot.sitar <- function(x, opt="dv", labels, apv=FALSE, xfun=NULL, yfun=NULL, subset=NULL,
                       ns=101, abc=NULL, add=FALSE, nlme=FALSE, ...,
                       xlab=NULL, ylab=NULL, vlab=NULL,
                       xlim=c(NA, NA), ylim=c(NA, NA), vlim=c(NA, NA)) {

  plotaxes <- function(dv, xlab, ylab, vlab, xlim, ylim, vlim)
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
    #	plot x & y1 axes
    plot(xlim, ylim, type='n', xlab=xlab, ylab=ylab)
    #	save x & y1 axis limits
    xy <- list()
    xy$usr <- par('usr')
    #	optionally add right axis
    if (dv == 3) {
      par(new=TRUE)
      plot(xlim, vlim, type='n', bty='n', ann=FALSE, axes=FALSE)
      axis(4)
      mtext(vlab, 4, par('mgp')[1])
      #	save y2 axis limits
      xy$usr2 <- par('usr')
      eval(parse(text=".par.usr2 <<- par('usr')"))
      #	reset axis limits
      par(usr=xy$usr)
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

  # stackage <- function(x, id, n=ns) {
  #   # generate x and id values across the x range to plot spline curves
  #   npt <- n / diff(range(x))
  #   xid <- by(data.frame(x=x, id=id), id, function(z) {
  #     nt <- floor(npt * diff(range(z$x))) + 1
  #     data.frame(x=seq(min(z$x), to=max(z$x), length.out=nt), id=rep.int(z$id[[1]], nt))
  #   })
  #   df <- xid[[1]][FALSE, ]
  #   for (dft in xid) df <- rbind(df, dft)
  #   df
  # }

  distance <- velocity <- function(model, n=ns) {
    dvt <- rlang::quo_name(match.call()[[1]])
    dvt <- as.numeric(dvt == 'velocity')
    .x <- getCovariate(model)[subset]
    .x <- xseq(.x, n)
    . <- tibble(
      .y=predict(model, tibble(.x), level=0, deriv=dvt, xfun=xfun, yfun=yfun),
      .x=xfun(.x)
    )
    if (length(.x[subset]) != length(.x)) attr(., 'subset') <- subset
    .
  }

  Distance <- Velocity <- function(model, n=ns) {
    # generate x and id values across the x range to plot spline curves
    dvt <- rlang::quo_name(match.call()[[1]])
    dvt <- as.numeric(dvt == 'Velocity')
    .x <- getCovariate(model)[subset]
    .id <- getGroups(model)[subset]
    npt <- n / diff(range(.x))
    . <- by(tibble(.x, .id), .id, function(z) {
      xrange <- range(z$.x)
      nt <- floor(npt * diff(xrange)) + 1
      tibble(
        .x=xseq(xrange, nt),
        .id=z$.id[[1]]
      )
    })
    . <- do.call('rbind', .)
    . <- mutate(.,
                .y=predict(model, ., deriv=dvt, xfun=xfun, yfun=yfun),
                .x=xfun(.x)
    )
    invisible(.)
  }

  unadjusted <- function(model) {
    tibble(
      .x=xfun(getCovariate(model)),
      .y=yfun(getResponse(model)),
      .id=getGroups(model)
    )
  }

  adjusted <- function(model) {
    . <- xyadj(model)
    tibble(
      .x=xfun(.$x),
      .y=yfun(.$y),
      .id=getGroups(model)
    )
  }

  effect <- function(model) {
    x <- getCovariate(model)[subset]
    x <- xseq(x, ns) - model$xoffset
    . <- tibble(
      .y=yfun(predict(model$ns, tibble(x))),
      .x=xfun(x + model$xoffset)
    )
    .
  }

  if (nlme) {
    do.call('plot.lme', as.list(match.call()[-1]))
  }
  else {
    model <- x
    data <- getData(model)
    mcall <- model$call.sitar
    x <- getCovariate(model)
    y <- getResponse(model)
    id <- getGroups(model)
    nf <- length(fitted(model))
    if (nf != length(y))
      stop(paste0('model (length=', nf, ') incompatible with data (nrows=', length(y), ')'))
    #	extract list(...)
    ccall <- match.call(expand.dots=FALSE)
    #	subset to plot model
    subset <- eval(ccall$subset, data, parent.frame())
    if (is.null(subset)) subset <- rep_len(TRUE, nf)
    # ... args
    dots <- ccall$...
    ARG <- if(!is.null(dots))
      lapply(as.list(dots), eval, data, parent.frame())
    else
      NULL

    options <- c('d', 'e', 'u', 'a', 'D', 'v', 'V')
    optnames <- c('distance', 'effect', 'unadjusted', 'adjusted', 'Distance', 'velocity', 'Velocity')
    optaxis <- c( 1,   1,   1,   1,   1,   2,   2 ) # default y1=1, y2=2
    optmult <- c( FALSE,   FALSE,   TRUE,   TRUE,   TRUE,   FALSE,   TRUE ) # multiple curves
    optsmooth <- c( TRUE,   TRUE,   FALSE,   FALSE,   TRUE,   TRUE,   TRUE ) # spline curves
    opts <- unique(na.omit(match(unlist(strsplit(opt, '')), options)))
    if (length(opts) == 0) stop('option(s) not recognised')
    dv <- range(optaxis[opts])
    dv <- min(dv) + diff(dv) * 2 # 1=d, 2=v, 3=dv
    mult <- any(optmult[opts]) # multiple curves

# create missing labels
    if (missing(labels))
      labels <- vector('character', 3)
    else if (length(labels) < 3)
      labels <- c(labels, '', '')
    if (labels[3] == '' && labels[2] != '')
      labels[3] <- paste(labels[2], 'velocity')

#	if xlab not specified replace with label or x name (depending on xfun)
    if (is.null(xlab)) {
      if(labels[1] == '') {
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
      if(labels[2] == '') {
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
      if(labels[3] == '')
        labels[3] <- paste(labels[2], 'velocity')
      vlab <- labels[3]
    }

    # set up ARG for left and right axes
    ARG1 <- ARG[names(ARG) != 'y2par']
    ARG2 <- ARG[['y2par']]
# default dotted line for velocity curve
  	if (is.null(ARG2$lty)) ARG2$lty <- 2

# derive xfun and yfun
    if (is.null(xfun))
      xfun <- ifun(mcall$x)
    if (is.null(yfun))
      yfun <- ifun(mcall$y)

# generate list of data frames for selected options
    data <- lapply(opts, function(i) {
      if (optsmooth[[i]])
        do.call(optnames[[i]], list(model=model, n=ns))
      else
        do.call(optnames[[i]], list(model=model))
    })

# extract axis ranges
    for (i in 1:length(opts)) {
      xlim <- range(xlim, data[[i]]$.x, na.rm=TRUE)
      if (optaxis[opts[[i]]] == 1)
        ylim <- range(ylim, data[[i]]$.y, na.rm=TRUE)
      else
        vlim <- range(vlim, data[[i]]$.y, na.rm=TRUE)
    }

# plot axes
    if (!add)
      xy <- plotaxes(dv, xlab, ylab, vlab, xlim, ylim, vlim)
    else
      xy <- .par.usr2

# plot curves
    lapply(1:length(opts), function(i) {
      opt <- opts[[i]]
      if (optaxis[opt] == 2 && dv == 3) {
        fun <- v2d(ylim, vlim)
        ARG <- ARG2
      } else {
        fun <- I
        ARG <- ARG1
      }
      if (optmult[opt])
        do.call("mplot", c(list(x=data[[i]]$.x, y=fun(data[[i]]$.y), id=data[[i]]$.id, add=TRUE), ARG))
      else
        xy <- do.call("lines", c(list(x=data[[i]]$.x, y=fun(data[[i]]$.y)), ARG))
    })

    if (TRUE) return(invisible(xy))

    # #	plot y vs t by subject
    # if (grepl("u", opt)) {
    #   xt <- x
    #   yt <- y
    #   do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG1))
    #   add <- TRUE
    # }
    #
    # #	plot fitted distance or velocity curves by subject
    # for (o in c('D', 'V')) {
    #   oV <- as.numeric(o == 'V')
    #   if (grepl(o, opt)) {
    #     newdata <- setNames(stackage(x[subset], id[subset]), c('.x', '.id'))
    #     newdata <- cbind(newdata, y=predict(model, newdata=newdata, xfun=xfun, yfun=yfun, deriv=oV))
    #     ylab <- labels[2 + oV]
    #     # adjust ARG=id to ARG=newid
    #     for (i in 1:length(ARG)) {
    #       arg <- ARG[[i]]
    #       if (length(arg) == length(id) &&
    #           length(lvl <- unlist(tapply(id, arg, unique))) == nlevels(id))
    #         ARG[[i]] <- lvl[newdata[, 2]]
    #     }
    #     do.call("mplot", c(list(x=xfun(newdata[, 1]), y=newdata[, 3], id=newdata[, 2],
    #                             data=newdata, add=add), ARG))
    #     add <- TRUE
    #   }
    # }
    #
    # #	plot fitted distance and velocity curves
    # if (grepl("d", opt) || grepl("v", opt) || apv) {
    #   xt <- xseq(x[subset])
    #   newdata <- data.frame(.x=xt)
    #   # if subset, flag for predict
    #   if (!identical(subset, rep_len(TRUE, nf))) attr(newdata, 'subset') <- subset
    #
    #   # distance and velocity
    #   xt <- xfun(xt)
    #   yt <- yfun(predict(object=model, newdata=newdata, level=0, abc=abc))
    #   vt <- predict(object=model, newdata=newdata, level=0, deriv=1, abc=abc, xfun=xfun, yfun=yfun)

    # save and print apv
    if (apv) {
      xy$apv <- setNames(getPeakTrough(xt, vt), c('apv', 'pv'))
      print(signif(xy$apv, 4))
    }

    #   # plot d &| v curve(s)
    #   if (grepl("d", opt) && grepl("v", opt))
    #     xy <- do.call("y2plot", c(list(x=xt, y1=yt, y2=vt, labels=labels, add=add, xy=xy), ARG))
    #   else {
    #     if (grepl("v", opt)) yt <- vt
    #     xy <- do.call("y2plot", c(list(x=xt, y1=yt, add=add, xy=xy), ARG1))
    #   }
    #   add <- TRUE
    # }
    #
    # #	plot fixed effects distance curve
    # if (grepl("e", opt)) {
    #   xt <- xseq(x[subset])
    #   yt <- predict(model$ns, newdata=data.frame(x=xt - model$xoffset))
    #   ox <- order(xt)
    #   xy <- do.call("y2plot", c(list(x=xfun(xt[ox]), y1=yfun(yt[ox]), add=add, xy=xy), ARG1))
    #   add <- TRUE
    # }
    #
    # #	plot adjusted y vs adjusted t by subject
    # if (grepl("a", opt)) {
    #   yt <- xyadj(model)
    #   xt <- yt$x
    #   yt <- yt$y
    #   do.call("mplot", c(list(x=xfun(xt), y=yfun(yt), id=id, subset=subset, add=add), ARG1))
    #   add <- TRUE
    # }

# plot vertical line(s) at age of peak velocity
    if (apv) {
      # multiple curves
      if (mult) {
        # x and y untransformed
        if (sum(x - xfun(x)) == 0 && sum(y - yfun(y)) == 0) {
          apv1 <- xyadj(model, xy$apv[1], tomean=FALSE)$x # subject-specific APVs
          pv1 <- xy$apv[2] * exp(ranef(model)$c) # subject-specific PVs
          xy$apv <- data.frame(apv=apv1, pv=pv1)
        }
        # xfun or yfun set
        else {
          xt <- xseq(x[subset])
          newdata <- data.frame(.x=xt)
          . <- vapply(levels(factor(id[subset])), function(z) {
            newdata$.id <- z
            vt <- predict(object=model, newdata=newdata, deriv=1, xfun=xfun, yfun=yfun)
            getPeakTrough(xfun(xt), vt)
          }, c(0, 0))
          xy$apv <- setNames(data.frame(t(.)), c('apv', 'pv'))
        }
      }
      if (is.null(ARG2$lty)) ARG2$lty <- 3
      if (add) do.call('abline', c(list(v=unlist(xy$apv['apv'])), ARG$y2par))
    }
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

