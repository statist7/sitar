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
#' variables after back-transforming any transformations in the original SITAR
#' call. So for example if \code{y = log(height)} in the SITAR call, then \code{yfun}
#' is applied to \code{height}. Thus the default \code{yfun = identity} has the effect of
#' back-transforming the SITAR call transformation - this is achieved by setting
#' \code{yfun = yfun(ifun(x$call.sitar$y))}.
#' For no transformation set \code{yfun = NULL}. The same applies to \code{xfun}.
#'
#' For models that include categorical fixed effects (e.g. \code{a.formula = ~sex + region})
#' the options 'dv' plot mean curves for each distinct group. Any continuous (as opposed
#' to grouped) fixed effect variables are set to their mean values in the plots, to ensure that the mean curves are
#' smooth. Setting \code{design} allows the grouping variables to be selected, e.g. \code{design = ~sex},
#' and \code{design = ~1} gives a single mean curve. The resulting plots can
#' be formatted with \code{par} in the usual way, indexed either by the individual grouping
#' variables (e.g. \code{sex} or \code{region} in the example) or the subject
#' factor \code{id} which indexes all the distinct plots.
#'
#' The helper functions \code{plot_d}, \code{plot_v}, \code{plot_D},
#' \code{plot_V}, \code{plot_u}, \code{plot_a} and \code{plot_c}
#' correspond to the seven plot \code{option}s defined by their last letter,
#' and return the data for plotting as a \code{tibble}, e.g. for use with
#' \code{ggplot2}. Setting \code{returndata = TRUE} works similarly
#' but handles multiple \code{option}s, returning a list of tibbles corresponding
#' to each specified \code{option}.
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
#' from positive to negative (after applying \code{xfun} and/or \code{yfun}). Age at peak velocity
#' is marked in the plot with a vertical dotted line, and its value, along with
#' peak velocity, is printed and returned. NB their standard errors can be
#' obtained using the bootstrap with the function \code{apv_se}. Values of \code{apv}
#' for individual subjects or groups are also returned invisibly.
#' @param xfun optional function to be applied to the x variable prior to
#' plotting (default identity, see Details).
#' @param yfun optional function to be applied to the y variable prior to
#' plotting (default identity, see Details).
#' @param subset optional logical vector of length \code{x} defining a subset
#' of \code{data} rows to be plotted, for \code{x} and \code{data} in the
#' original \code{sitar} call.
#' @param ns scalar defining the number of points for spline curves
#' (default 101).
#' @param design formula defining the variables to use to group data for multiple
#' mean distance and/or velocity curves (\code{opt = 'dv'}). By default includes
#' all the categorical variables named in \code{a.formula}, \code{b.formula},
#' \code{c.formula} and \code{d.formula}.
#' @param abc vector of named values of random effects a, b, c and d used to
#' define an individual growth curve, e.g. abc = c(a = 1, c = -0.1). Alternatively a
#' single character string defining an \code{id} level whose random effect
#' values are used. If \code{abc} is set, \code{level} is ignored. If
#' \code{abc} is NULL (default), or if a, b, c or d values are missing, values of
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
#' containing the data to be plotted. The helper functions each return a tibble
#' where the first three variables are '.x', '.y' and '.id', plus
#' variable '.groups' for curves grouped by \code{design}) and other covariates in the model.
#' Note that '.x' and '.y' are returned
#' after applying \code{xfun} and \code{yfun}. Hence if for example \code{x = log(age)}
#' in the SITAR call then '.x' corresponds by default to \code{age}.
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @seealso \code{\link{mplot}},
#' \code{\link{plotclean}}, \code{\link{ifun}}, \code{\link{apv_se}}
#' @keywords aplot
#' @examples
#'
#' ## fit sitar model
#' m1 <- sitar(x = age, y = height, id = id, data = heights, df = 4)
#'
#' ## draw fitted distance and velocity curves
#' ## with velocity curve in blue
#' ## adding age at peak velocity (apv)
#' plot(m1, y2par = list(col = 'blue'), apv = TRUE)
#'
#' ## bootstrap standard errors for apv and pv
#' \dontrun{
#' res <- apv_se(m1, nboot = 20, plot = TRUE)
#' }

#' ## draw individually coloured growth curves adjusted for random effects
#' ## using same x-axis limits as for previous plot
#' plot(m1, opt = 'a', col = id, xlim = xaxsd())
#'
#' ## add mean curve in red
#' lines(m1, opt = 'd', col = 'red', lwd = 2)
#'
#' ## add mean curve for a, b, c = -1 SD
#' lines(m1, opt = 'd', lwd = 2, abc = -sqrt(diag(getVarCov(m1))))
#'
#' ## use subset to plot mean curves by group
#' ## compare curves for early versus late menarche
#' heights <- within(sitar::heights, {
#'   men <- abs(men)
#'     late <- factor(men > median(men))
#'   })
#' # fit model where size and timing differ by early vs late menarche
#' m2 <- sitar(log(age), height, id, heights, 5,
#'   a.formula = ~late, b.formula = ~late)
#' ## early group
#' plot(m2, subset = late == FALSE, col = 4, lwd = 3,
#'   y2par = list(col = 4, lwd = 2), ylim = range(heights$height))
#' ## late group
#' lines(m2, subset = late == TRUE, col = 2, lwd = 3,
#'   y2par = list(col = 2, lwd = 2))
#' ## add legend
#' legend('right', paste(c('early', 'late'), 'menarche'),
#'   lty = 1, col = c(4, 2), inset = 0.04)
#'
#' ## alternatively plot both groups together
#' plot(m2, lwd = 3, col = late, y2par = list(lwd = 3, col = late))
#' legend('right', paste(c('early', 'late'), 'menarche'),
#'   lwd = 3, col = 1:2, inset = 0.04)

#' ## draw fitted height distance curves coloured by subject, using ggplot
#' \dontrun{
#' require(ggplot2)
#' ggplot(plot_D(m1), aes(.x, .y, colour = .id)) +
#' labs(x = 'age', y = 'height') +
#' geom_line(show.legend = FALSE)
#' }

#' @importFrom grDevices xy.coords
#' @importFrom graphics plot axis identify legend lines locator par text title mtext abline
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr mutate rename filter nest_by slice
#' @importFrom tidyr expand_grid
#' @importFrom rlang .data as_label %||%
#' @importFrom glue glue
#' @export
plot.sitar <- function(x, opt="dv", labels=NULL, apv=FALSE, xfun=identity, yfun=identity, subset=NULL,
                       ns=101, design=NULL, abc=NULL, trim=0, add=FALSE, nlme=FALSE,
                       returndata=FALSE, ...,
                       xlab=NULL, ylab=NULL, vlab=NULL,
                       xlim=c(NA, NA), ylim=c(NA, NA), vlim=c(NA, NA),
                       legend=list(x='topleft', inset=0.04, bty='o')) {

  plotaxes <- function(data, dv, xlab, ylab, vlab, xlim, ylim, vlim, ...)
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

  distance <- velocity <- function(model, subset=subset, xfun=xfun, yfun=yfun, abc=abc, ns=ns,
                                   design=design) {
# generate x values across the range to plot mean spline curve
    dvt <- as_label(match.call()[[1]])
    dvt <- as.numeric(dvt == 'velocity')
    level <- as.numeric(!is.null(abc))
    design_names <- all.vars(design)
    .x <- xseq(getCovariate(model), ns)
    newdata <- tibble()
    if (length(design_names) > 0L) {
      newdata <- getData(model)[subset, ] %>%
        select(all_of(design_names)) %>%
        select(!where(is.numeric)) %>%
        unique()
      if (length(newdata) > 0)
        newdata <- newdata %>%
          mutate(.groups = factor(paste(!!!lapply(names(.), as.name), sep = '_')), .before = 1) %>%
          expand_grid(.x, .)
    }
    if (length(newdata) == 0)
      newdata <- tibble(.x)
    if (sum(subset) < length(subset))
      attr(newdata, 'subset') <- subset
    newdata <- newdata %>%
      mutate(.y = predict(model, ., level=level, deriv=dvt, abc=abc, xfun=xfun, yfun=yfun), .after = .x,
             .x = xfun(.x),
             .id = getGroups(model)[1])
    newdata
  }

  Distance <- Velocity <- function(model, subset=subset, xfun=xfun, yfun=yfun, abc=abc, ns=ns,
                                   design=design) {
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
      select(-.data$v) %>%
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

  crosssectional <- function(model, subset=subset, xfun=xfun, yfun=yfun, abc=abc, ns=ns,
                             design=design) {
# fixed effect mean curve
    x <- getCovariate(model)[subset]
    x <- xseq(x, ns) - model$xoffset
    . <- tibble(
      .y = yfun(predict(model$ns, tibble(x))),
      .x = xfun(x + model$xoffset),
      .id = getGroups(model)[1]
    )
    .[, c('.x', '.y', '.id')]
  }

  dolegend <- function(ARG1, ARG2, legend) {
# add legend
    parlu <- function(...) par(do.call('par', list(...)))
    ARG2$lty <- ARG2$lty %||% 2
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
    if (fun == 'identity')
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
  subset <- eval(ccall$subset, data, parent.frame()) %||% rep_len(TRUE, model$dims$N)
# ... args
  dots <- ccall$...
  ARG <- if (!is.null(dots))
    lapply(as.list(dots), eval, data[subset, ], parent.frame())

  # derive design
  design <- design %||% as.list(mcall)[grepl('.formula', names(mcall))] %>%
    asOneFormula()
  stopifnot('design should be a formula' = is.language(design))

  # create labels
  labels <- labels %||% vector('character', 3)

  #	get axis labels
  xlab <- getlab(xlab, labels[1], mcall$x, paste(deparse(substitute(xfun)), collapse=''))
  ylab <- getlab(ylab, labels[2], mcall$y, paste(deparse(substitute(yfun)), collapse=''))
  vlab <- getlab(vlab, labels[3], ylab, 'velocity')

  # get xfun and yfun
  xfun <- getfun(xfun, mcall$x)
  yfun <- getfun(yfun, mcall$y)

  pt <- tibble(
    options   = c('d', 'c', 'u', 'a', 'D', 'v', 'V'),
    optnames  = c('distance', 'crosssectional', 'unadjusted', 'adjusted', 'Distance', 'velocity', 'Velocity'),
    optdv   = c( 1,   1,   1,   1,   1,   2,   2 ), # distance or velocity
    optmult   = c( FALSE, FALSE, TRUE,  TRUE,  TRUE,  FALSE, TRUE ), # multiple curves
    optsmooth = c( TRUE,  TRUE,  FALSE, FALSE, TRUE,  TRUE,  TRUE ) # spline curves
  )

  pt <- pt %>%
    slice(opt %>% strsplit('') %>% unlist %>% match(options) %>% unique)
  stopifnot('no options recognised' = nrow(pt) > 0L)
  dv <- pt %>% pull(.data$optdv) %>% range %>% unique %>% sum # 1 = d, 2 = v, 3 = dv

  # names of x y id vars
  xyid <- setNames(1:3, map_chr(mcall[2:4], all.vars))
  nid <- as.name(names(xyid)[3])

  # generate tibble of data frames for selected options
  pt <- pt %>%
    rowwise %>%
    mutate(optdv = if_else(.data$optdv == 2 & dv == 3, 3, .data$optdv), # velocity on y2 axis
           data = ifelse(.data$optsmooth,
                         list(do.call(.data$optnames, list(model=model, subset=subset, xfun=xfun, yfun=yfun,
                                                     abc=abc, ns=ns, design=design))),
                         list(do.call(.data$optnames, list(model=model, subset=subset, xfun=xfun, yfun=yfun,
                                                     trim = trim)))),
           xlim = list(range(data$.x, na.rm = TRUE)),
           ylim = list(range(data$.y, na.rm = TRUE)),
           groups = '.groups' %in% names(data),
           data = if_else(.data$groups, list(data %>% select(-3)), list(data)),
           optmult = if_else(.data$groups, TRUE, .data$optmult),
           data = list(data %>% rename(all_of(xyid))))

  # return data?
  if (returndata){
    data <- structure(pt %>% pull(data), names = pt %>% pull(options))
    if (length(data) == 1)
      data <- data[[1]]
    return(invisible(data))
  }

  # extract axis ranges and plot axes
  if (!add) {
    if (any(is.na(xlim)))
      xlim <- pt %>% pull(xlim) %>% unlist %>% range(., na.rm = TRUE)
    if (any(is.na(ylim)) && any(pt %>% pull(.data$optdv) == 1))
      ylim <- pt %>% filter(.data$optdv == 1) %>% pull(ylim) %>% unlist %>% range(., na.rm = TRUE)
    if (any(is.na(vlim)) && any(pt %>% pull(.data$optdv) > 1))
      vlim <- pt %>% filter(.data$optdv > 1) %>% pull(ylim) %>% unlist %>% range(., na.rm = TRUE)

    xy <- do.call('plotaxes',
                  c(list(data = pt %>% pull(data), dv = dv, xlab = xlab, ylab = ylab, vlab = vlab,
                         xlim = xlim, ylim = ylim, vlim = vlim), ARG))
    # add legend
    if (dv == 3 && !is.null(legend)) {
      legend[['legend']] <- c(ylab, vlab)
      dolegend(ARG[names(ARG) != 'y2par'], ARG$y2par, legend)
    }
    # else retrieve axis ranges
  } else {
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
  dots1 <- dots[names(dots) != 'y2par']
  dots2 <- dots[['y2par']]
  dots2 <- as.list(dots2)[-1]
  dots2[['lty']] <- dots2[['lty']] %||% 2
  pt <- pt %>%
    rowwise %>%
    mutate(ARG1 = list(lapply(as.list(dots1), eval, data)),
           ARG2 = list(lapply(as.list(dots2), eval, data)),
           ARG = ifelse(.data$optdv == 3, list(.data$ARG2), list(.data$ARG1)),
           fun = list(ifelse(.data$optdv == 3, v2d(!!ylim, !!vlim), identity)),
           plot = ifelse(.data$optmult,
                         list(do.call("mplot", c(list(x = data[[1]], y = .data$fun(data[[2]]), id = data[[3]], add = TRUE), ARG))),
                         list(do.call("lines", c(list(x = data[[1]], y = .data$fun(data[[2]])), ARG)))))

  # save and print vertical line(s) at age of peak velocity
  if (apv) {
    # single curve
    xy$apv <- with(velocity(model, subset=subset, xfun=xfun, yfun=yfun, abc=abc, ns=ns, design=~1),
                   setNames(getPeak(.x, .y), c('apv', 'pv')))
    print(signif(xy$apv, 4))
    # multiple smooth curves or abc
    pt <- pt %>%
      ungroup %>%
      filter(.data$optmult & .data$optsmooth) %>%
      select(-starts_with('opt')) %>%
      slice(1)
    # multiple dv groups or id groups
    level <- pt %>% pull(.data$groups) %>% `!` %>% as.numeric
    if (nrow(pt) > 0)
      xy$apv <- pt$data[[1]] %>%
        nest_by({{nid}}, .keep = TRUE) %>%
        mutate(vel = list(predict(model, data, level = level, deriv = 1, abc=abc, xfun = xfun, yfun = yfun)),
               xy = list(getPeak(xfun(data[[1]]), .data$vel) %>% t %>% as_tibble)) %>%
        select(-c(.data$data:.data$vel)) %>%
        unnest(xy) %>%
        ungroup %>%
        select({{nid}}, apv = x, pv = .data$y)

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
