#' Fit SITAR growth curve model
#'
#' SITAR is a method of growth curve analysis, based on \pkg{nlme}, that
#' summarises a set of growth curves with a mean growth curve as a regression
#' spline, plus a set of up to four fixed and random effects (a, b, c and d)
#' defining how individual growth curves differ from the mean curve.
#'
#' The SITAR model usually has up to three random effects (a, b and c), termed
#' size, timing and intensity respectively. \code{df} sets the degrees of
#' freedom for the mean spline curve, taking values from 1 (i.e. linear)
#' upwards. In addition there is a random effect for the slope, d, which is
#' fitted when \code{df = 0}, and combined with a, it provides the classic
#' random intercept random slope model, which is similar to the 1 df spline
#' model. In addition d can be fitted, along with a, b and c, to extend SITAR to
#' model variability in the adult slope of the growth curve.
#'
#' \code{xoffset} allows the origin of \code{x} to be varied, while
#' \code{bstart} specifies the starting value for \code{b}, both of which can
#' affect the model fit and particularly \code{b}. The values of \code{bstart},
#' \code{knots} and \code{bounds} are offset by \code{xoffset} for fitting
#' purposes, and similarly for fixed effect \code{b}.
#'
#' The formulae \code{a.formula}, \code{b.formula}, \code{c.formula} and
#' \code{d.formula} allow for cov.names and can include functions and
#' interactions. \code{\link{make.names}} is used to ensure that the names of
#' the corresponding model terms are valid. The modified names, not the original
#' names, need to be specified in \code{predict.sitar}.
#'
#' \code{update} updates the model by taking the \code{object} call, adding any
#' new parameters and replacing changed ones. Where feasible the fixed and
#' random effects of the model being updated are suitably modified and passed
#' via the \code{start} argument.
#'
#' The mean curve is a natural cubic spline, fitted originally with
#' \code{splines::ns} but now by default with \code{splines2::nsk}, where the
#' choice is set by \code{stype}. Both give the same mean curve, but \code{nsk}
#' is faster than \code{ns} and its fixed effect coefficients have a direct
#' interpretation: added to fixed effect \code{a} they give the values of the
#' mean curve at the \code{df} right-most knots, while \code{a} is the value at
#' the left-most knot. The knots themselves depend on fixed effects \code{b} and
#' \code{c} (where fitted) and are accessed via \code{attr(object, "knots")}.
#'
#' @aliases sitar update.sitar
#' @param x vector of ages.
#' @param y vector of measurements.
#' @param id factor of subject identifiers.
#' @param data data frame containing variables \code{x}, \code{y} and \code{id}.
#' @param df degrees of freedom for cubic regression spline (0 or more, see
#'   Details).
#' @param knots vector of values for knots (default \code{df} quantiles of
#'   \code{x} distribution).
#' @param fixed character string specifying a, b, c, d fixed effects (default
#'   \code{random} or the subset of "a + b + c + d" within \code{random}).
#' @param random character string specifying a, b, c, d random effects (default
#'   \code{"a+b+c"}). Alternatively \code{nlme} formula e.g. \code{"list(id =
#'   pdDiag(a + b + c ~ 1))"}.
#' @param pdDiag logical which if TRUE fits a diagonal random effects covariance
#'   matrix, or if FALSE (default) a general covariance matrix.
#' @param a.formula formula for fixed effect a (default \code{~ 1}).
#' @param b.formula formula for fixed effect b (default \code{~ 1}).
#' @param c.formula formula for fixed effect c (default \code{~ 1}).
#' @param d.formula formula for fixed effect d (default \code{~ 1}).
#' @param d.adjusted logical defining x scale for random effect d (default FALSE
#'   means x unadjusted, TRUE means x adjusted with random effects b and c).
#' @param stype version of natural spline curve to be fitted, default
#'   \code{"nsk"} or \code{"ns"}. See Details.
#' @param bounds span of \code{x} for regression spline, or fractional extension
#'   of range (default 0.04).
#' @param start optional numeric vector of initial estimates for the fixed
#'   effects, or list of initial estimates for the fixed and random effects (see
#'   \code{\link[nlme]{nlme}}).
#' @param xoffset optional value of offset for \code{x} (either "mean"
#'   (default), "apv" or value).
#' @param bstart optional starting value for fixed effect \code{b} (either
#'   "mean", "apv" or value (default \code{xoffset})).
#' @param returndata logical which if TRUE causes the model matrix to be
#'   returned, or if FALSE (default) the fitted model. Setting returndata TRUE
#'   is useful in conjunction with \code{\link{subset}} and \code{\link{subsample}} for
#'   simulation purposes.
#' @param verbose optional logical value to print information on the evolution
#'   of the iterative algorithm (see \code{\link[nlme]{nlme}}).
#' @param correlation optional \code{corStruct} object describing the
#'   within-group correlation structure (see \code{\link[nlme]{nlme}}).
#' @param weights optional \code{varFunc} object or one-sided formula describing
#'   the within-group heteroscedasticity structure (see \code{\link[nlme]{nlme}}).
#' @param subset optional expression indicating the subset of the rows of data
#'   that should be used in the fit (see \code{\link[nlme]{nlme}}).
#' @param method character string, either "REML" or "ML" (default) (see
#'   \code{\link[nlme]{nlme}}).
#' @param na.action function for when the data contain NAs (see
#'   \code{\link[nlme]{nlme}}).
#' @param control list of control values for the estimation algorithm (see
#'   \code{\link[nlme]{nlme}}) (default \code{nlmeControl(returnObject = TRUE)}).
#' @param keep.data logical to control saving \code{data} as part of the model
#'   object (default TRUE).
#' @param object object of class \code{sitar}.
#' @param \dots further parameters for \code{update} consisting of any of the
#'   above \code{sitar} parameters.
#' @param evaluate logical to control evaluation. If TRUE (default) the
#'   expanded \code{update} call is passed to \code{sitar} for evaluation, while
#'   if FALSE the expanded call itself is returned.
#'
#' @return An object inheriting from class \code{sitar} representing the
#'   nonlinear mixed-effects model fit, with all the components returned by
#'   \code{nlme} (see \code{\link[nlme]{nlmeObject}} for a full description) plus the
#'   following components: \item{fitnlme}{the function returning the predicted
#'   value of \code{y}.} \item{data}{copy of \code{data} (if \code{keep.data}
#'   true).}
#' \item{constants}{data frame of mean a-b-c-d values for unique combinations
#'   of covariates (excluding \code{x}).}
#'   \item{call.sitar}{the internal \code{sitar} call that produced the object.}
#'   \item{xoffset}{the value of \code{xoffset}.}
#' \item{ns}{the \code{lm} object providing starting values for the B-spline
#'  curve.}
#'
#'   In addition attributes \code{stype} and \code{knots} are set.
#'
#'   Generic functions such as \code{print}, \code{plot}, \code{anova} and
#'   \code{summary} have methods to show the results of the fit. The functions
#'   \code{residuals}, \code{coef}, \code{fitted}, \code{fixed.effects},
#'   \code{random.effects}, \code{predict}, \code{getData}, \code{getGroups},
#'   \code{getCovariate} and \code{getVarCov} can be used to extract some of its
#'   components.
#'
#' @author Tim Cole \email{tim.cole@@ucl.ac.uk}
#' @keywords package nonlinear regression models
#' @examples
#'
#' data(heights)
#' ##  fit simple model
#' (m1 <- sitar(x=age, y=height, id=id, data=heights, df=5))
#'
#' ##  relate random effects to age at menarche (with censored values +ve)
#' ##  both a (size) and b (timing) are positively associated with age at menarche
#' (m2 <- update(m1, a.formula = ~abs(men), b.formula = ~abs(men), c.formula = ~abs(men)))
#' @import nlme
#' @importFrom glue glue
#' @importFrom splines ns
#' @importFrom splines2 nsk
#' @importFrom stats AIC BIC as.formula coef cor fitted lm logLik mad
#'   model.frame model.matrix na.fail na.omit pnorm predict qnorm quantile
#'   residuals sd setNames smooth.spline spline update update.formula
#' @export
sitar <-
  function(x,
           y,
           id,
           data,
           df,
           knots,
           fixed = NULL,
           random = 'a + b + c',
           pdDiag = FALSE,
           a.formula =  ~ 1,
           b.formula =  ~ 1,
           c.formula =  ~ 1,
           d.formula =  ~ 1,
           d.adjusted = FALSE,
           stype = c('nsk', 'ns'),
           bounds = 0.04,
           start,
           xoffset = 'mean',
           bstart = xoffset,
           returndata = FALSE,
           verbose = FALSE,
           correlation = NULL,
           weights = NULL,
           subset = NULL,
           method = 'ML',
           na.action = na.fail,
           control = nlmeControl(msMaxIter = 100, returnObject = TRUE),
           keep.data = TRUE) {

    stype <- match.arg(stype)
    splinefun <- get(stype)

    b.origin <- function(b) {
      if (b == 'mean')
        return(mean(x))
      if (b == 'apv') {
        spline.lm <- lm(y ~ splinefun(x, knots = knots, Bound = bounds))
        return(getPeak(x, predict(
          smooth.spline(x, fitted(spline.lm)), x, deriv = 1
        )$y)[1])
      }
      if (!is.numeric(b))
        stop('must be "mean" or "apv" or numeric')
      b
    }

    # get data
    mcall <- match.call()
    data <- eval.parent(mcall$data)
    subset <- eval(mcall$subset, data)
    if (!is.null(subset))
      data <- data[subset, ]
    x <- eval(mcall$x, data)
    y <- eval(mcall$y, data)

    # get df, knots and bounds
    if (missing(df) &&
        missing(knots))
      stop("either df or knots must be specified")
    if (!missing(df) &&
        !missing(knots)) {
      if (df == 0)
        knots <- NULL
      else
        warning("both df and knots specified - df redefined from knots\n")
    }
    if (!missing(knots)) {
      if (!identical(range(c(knots, x)), range(x)))
        stop("knots outside x range")
      knots <- sort(knots)
      df <- length(knots) + 1
      mcall$df <- NULL
    } else {
      df <- round(df)
      if (df < 0)
        stop("df must be 0 or more")
      knots <- if (df > 1)
        quantile(x, (1:(df - 1)) / df)
      else
        numeric(0)
    }
    if (nrow(data) <= df)
      stop("too few data to fit spline curve")
    if (length(bounds) == 1)
      bounds <- range(x) + abs(bounds) * c(-1, 1) * diff(range(x))
    else if (length(bounds) != 2)
      stop("bounds should be length 1 or 2")

    xoffset <- b.origin(xoffset)
    bstart <- b.origin(bstart) - xoffset
    x <- x - xoffset
    knots <- knots - xoffset
    bounds <- bounds - xoffset

    # extract and format random and fixed
    extract <- function(x, set = letters[1:4]) {
      x <- paste(as.character(x), collapse = ' ')
      x <- unique(strsplit(x, '[^a-z]')[[1]]) # extract words
      paste(x[x %in% set], collapse = ' + ') # combine selected
    }

    # if random contains ~ extract effects from formula
    if (any(grepl('~', random))) {
      fullrandom <- ifelse(is.character(random),
                           random,
                           deparse(random))
    } else {
      fullrandom <- ifelse(pdDiag,
                           glue('list(id = pdDiag({random} ~ 1))'),
                           NA)
    }
    random <- extract(random)

    # default fixed effects
    if (is.null(fixed))
      # drop d from fixed if fixed missing and d in random and df > 1
      if (grepl('d', random) && df > 1)
        fixed <- extract(random, letters[1:3])
      else
        fixed <- random
    # else format fixed
    else
      fixed <- extract(fixed)

    #	force fixed effect and intercept for a
    fixed <- extract(paste('a', fixed))
    a.formula <- update.formula(a.formula, ~. + 1)

    # adjust fixed and random for df = [01] and set start
    if (df > 0) { # fit spline
      ss <- paste0('s', 1:df)
    #	if start missing get start values for ss and a
      spline.lm <- lm(y ~ splinefun(x, knots = knots, Bound = bounds))
      if (nostart <- missing(start))
        start <- coef(spline.lm)[c(2:(df + 1), 1)]
    # if df = 1 exclude b and d
      if (df == 1) {
        fixed <- 'a' # s1 = c
        random <- extract(random, c('a', 'c'))
      }
    } else { # 0 df for spline so fit y ~ x
      ss <- character(0)
      #	if start missing get start values for [ad]
      spline.lm <- lm(y ~ x)
      if (nostart <- missing(start))
        start <- coef(spline.lm)
      fixed <- 'a + d'
      if (!'random' %in% names(mcall))
        random <- fixed
      else
        random <- extract(random, c('a', 'd'))
    }
    if (nchar(random) == 0)
      stop('no random effects')

    # set up args for fitnlme
    fix <- fixed
    fixed <- ss
    random.names <- all.names(as.formula(paste0('~', random)), functions = FALSE)
    pars <- c('x', random.names, ss)

    #	set up model elements for a, b, c and d
    model <- setNames(nm = letters[1:4])
    model[!names(model) %in% random.names] <- ''

    # set up fulldata
    id <- eval(mcall$id, data)
    subset <- 1:length(x)
    fulldata <- data.frame(x, y, id, subset)

    mm.formula <- ~1
    cmm <- matrix(nrow = nrow(data), ncol = 0)
    for (l in names(model)) {
      formula <- update.formula(get(paste(l, 'formula', sep = '.')), ~.)
      if ((!grepl(l, fix) && formula == ~1) || formula == ~1 - 1)
        next
      if (formula == ~1)
        mm.intercept <- TRUE
      else {
        if (formula != mm.formula) {
          mm <- model.matrix(formula, data)
          if (nrow(mm) < length(y))
            stop('Missing values in data')
          mm.formula <- formula
          mm.intercept <- colnames(mm)[1] == '(Intercept)'
          # ensure names are valid
          colnames(mm) <- make.names(colnames(mm), unique = TRUE)
          # omit constant columns
          mm <- mm[, apply(mm, 2, sd) > 0, drop = FALSE]
        }
        if (exists('mm'))
          for (i in 1:ncol(mm)) {
            var <- colnames(mm)[i]
            rc <- paste(l, var, sep = '.')
            pars <- c(pars, rc)
            fixed <- c(fixed, rc)
            if (nostart)
              start <- c(start, 0)
            model[l] <- paste0(model[l], '+', rc, '*', var)
            if (!var %in% pars) {
              pars <- c(pars, var)
              cmm <- cbind(cmm, mm[, i, drop = FALSE])
            }
          }
      }
      if (mm.intercept) {
        fixed <- c(fixed, l)
        if (!l %in% pars) {
          pars <- c(pars, l)
          model[l] <- paste0(l, model[l])
        }
        if (nostart) {
          if (l == 'b')
            start <- c(start, bstart)
          else if ((l == 'c' || l == 'd') && grepl(l, fix) && df > 0)
            start <- c(start, 0)
        }
      }
    }
    # centre covariate columns
    cmm <- scale(cmm, scale = FALSE)
    # combine with data
    fulldata <- cbind(fulldata, cmm)
    # save covariate means for predict
    attr(fulldata, 'scaled:center') <- attr(cmm, 'scaled:center')

    # 	if (!is.null(weights)) {
    #     if (is.list(weights)) form <- asOneFormula(lapply(weights, function(z) attr(z, 'formula')))
    #       else form <- attr(weights, 'formula')
    #     if (!is.null(form)) {
    #       wt <- as.data.frame(model.frame(form, data))
    #       wt.names <- names(wt)[!names(wt) %in% names(fulldata)]
    #       wt <- as.data.frame(wt[, wt.names])
    #       names(wt) <- wt.names
    #       fulldata <- cbind(fulldata, wt)
    # 	  }
    # 	}

    if (returndata)
      return(invisible(fulldata))
    if (nostart)
      names(start) <- fixed
    fixed <- paste(fixed, collapse = '+')
    pars <- paste(pars, collapse = ',')
    ss <- paste(ss, collapse = ',')

    #	combine model elements
    cglue <- function(x, start, end) {
      ifelse (x == '', x, glue(start, x, end))
    }

    nsa <- cglue(model['a'], '', '')
    nsb <- cglue(model['b'], '-(', ')')
    nsc <- cglue(model['c'], ')*exp(', '')
    dx <- ifelse(d.adjusted, 'ex', 'x') # set scale for d
    nsd <- cglue(model['d'], '+(', ')*{dx}')
    ex <- glue('(x{nsb}{nsc})')
    spline <- glue(cglue(ss, '+rowSums((cbind(', glue(')*', stype, '(ex,k=knots,B=bounds)))')))
    # expand fixed and if necessary random
    fixed <- glue('{fixed} ~ 1')
    random <- ifelse (is.na(fullrandom),
                      glue('{random} ~ 1 | id'),
                      fullrandom)

    #	code to parse
    fitcode <- glue(
      "fitenv <- new.env()\n",
      "fitenv$fitnlme <- function(<<pars>>) {\n",
      "ex <- <<ex>>\n",
      "<<nsa>><<spline>><<nsd>>\n",
      "}\n",
      "on.exit(detach(fitenv))\n",
      "attach(fitenv)\n",
      "nlme(y ~ fitnlme(<<pars>>),",
      "fixed = <<fixed>>,",
      "random = <<random>>,\n",
      "data = fulldata,",
      "start = start, correlation = correlation,",
      "weights = weights, subset = subset, method = method,\n",
      "na.action = na.action, control = control, verbose = verbose)",
      .open = "<<",
      .close = ">>"
    )

    #	print values
    if (verbose) {
      cat('\nconstructed code', fitcode, sep = '\n')
      cat(
        '\ndf',
        df,
        'bstart',
        bstart,
        'xoffset',
        xoffset,
        '\nknots\n',
        knots,
        '\nbounds\n',
        bounds
      )
      if (is.list(start)) {
        cat('\nstarting values\n  fixed effects\n', start$fixed)
        if (!is.null(start$random)) {
          cat('\n  random effects\n')
          print(start$random)
        }
      }
      else
        cat('\nstarting values\n', start, '\n')
    }

    #	save fitted model
    nlme.out <- eval(parse(text = fitcode))
    class(nlme.out) <- c('sitar', class(nlme.out))
    # save fitnlme
    nlme.out$fitnlme <- fitenv$fitnlme
    # save data
    if (keep.data)
      nlme.out$data <- data
    # save constants
    model <- lapply(model, function(x)
      if (x == '')
        NULL
      else
        as.formula(paste0('~', x)))
    model <- model[!sapply(model, is.null)]
    cov.names <- unique(unlist(lapply(names(model), function(x)
      all.names(get(paste0(x, '.formula')), functions = FALSE)
    )))
    random.mean <- if (is.data.frame(ranef(nlme.out)))
      apply(ranef(nlme.out), 2, mean)
    else # fudge for random as list
      setNames(rep(0, length(random.names)), random.names)
    fixed <- fixef(nlme.out)
    fixed[random.names[!random.names %in% names(fixed)]] <- 0
    fixed[random.names] <- fixed[random.names] + random.mean
    newdata <- data.frame(fulldata, t(fixed))
    model <- as.data.frame(lapply(model, function(x) eval(x[[2]], newdata)))
    if (length(cov.names) > 0L) {
      newdata <- setNames(as.data.frame(lapply(cov.names, function(x) with(data, get(x)))),
                          cov.names)
      newdata <- cbind(newdata, model)
      for (i in rev(seq.int(cov.names)))
        newdata <- newdata[order(newdata[, i]), ]
    } else
      newdata <- model
    newdata <- unique(newdata)
    rownames(newdata) <- 1:nrow(newdata)
    nlme.out$constants <- newdata
    # save call
    nlme.out$call.sitar <- mcall
    # save xoffset
    nlme.out$xoffset <- xoffset
    # save spline type
    attr(nlme.out, 'stype') <- stype
    # save knots on original scale to match nsk fixed effects
    fe <- fixef(nlme.out)
    fe[setdiff(letters[2:3], names(fe))] <- 0
    aknots <- sort(c(knots, bounds))
    aknots <- ifun(mcall$x)(aknots / exp(fe[['c']]) + fe[['b']] + xoffset)
    attr(nlme.out, 'knots') <- aknots
    # save ns curve
    nlme.out$ns <- spline.lm
    # save d.adjusted
    if ('d' %in% random.names)
      attr(nlme.out, 'd.adjusted') <- d.adjusted
    #   if (exists('start.')) rm(start., inherits=TRUE)
    nlme.out
  }

#' @rdname sitar
#' @method update sitar
#' @export
update.sitar <- function (object, ..., evaluate = TRUE)
{
  mcall <- object$call.sitar
  if (is.null(mcall))
    stop("need an object with call.sitar component")
  extras <- as.list(match.call(expand.dots = FALSE)$...)
  mcall$start <- NULL
  if (exists('data', object)) {
    mcall$data <- NULL
    mcall$subset <- NULL
  }
  #	expand formulae
  if (any(grep('formula', names(extras)))) {
    for (n in paste0(letters[1:4], '.formula')) {
      if (!is.null(extras[[n]]) && !is.null(mcall[[n]]))
        extras[[n]] <- update.formula(mcall[[n]], extras[[n]])
    }
  }
  # update args
  mcall[names(extras)] <- extras
  # check if mean curve is/was a spline with 2+ df
  fo <- fixef(object)
  ro <- ranef(object)
  spline2 <- !is.na(fo['s2']) && (is.null(extras$df) || eval.parent(extras$df) >= 2)
  # if df = 0 drop knots & bounds
  if (!is.null(mcall$df) && eval.parent(mcall$df) == 0) {
    mcall$knots <- mcall$bounds <- NULL
  }
  #	add start arg if none of these args specified and is/was spline curve with 2+ df
  if (!sum(pmatch(
    names(extras),
    c(
      "x",
      "y",
      "id",
      "fixed",
      "random",
      "a.formula",
      "b.formula",
      "c.formula",
      "d.formula",
      "start",
      "returndata"
    ),
    0
  )) && spline2) {
    start. <- list(fixed = fo, random = ro)
    # update start if any of these args specified
    if (sum(pmatch(
      names(extras),
      c(
        'data',
        'subset',
        'df',
        'knots',
        'bounds',
        'xoffset',
        'bstart'
      ),
      0
    ))) {
      # get data
      if (any(c('data', 'subset') %in% names(mcall))) {
        data <- if (is.null(mcall$data))
          object$data
        else
          eval.parent(mcall$data)
        subset <- eval(mcall$subset, data)
        if (!is.null(subset))
          data <- data[subset, ]
      } else
        data <- getData(object)
      x <- eval(mcall$x, data)
      xoffset <- object$xoffset
      if (is.null(xoffset))
        xoffset <- mean(x)
      x <- x - xoffset
      df <- object$ns$rank - 1
      knots <- attr(object$ns$model$splinefun, 'knots')
      bounds <- attr(object$ns$model$splinefun, 'Boundary.knots')
      # update random effects
      if (any(c('data', 'subset') %in% names(mcall))) {
        id <- factor(eval(mcall$id, data))
        levels.obj <- levels(getGroups(object))
        if (!identical(levels(id), levels.obj)) {
          #	omit random effects for missing levels in id
          start.$random <-
            start.$random[idcheck <- levels.obj %in% levels(id),]
          cat(length(levels.obj) - sum(idcheck),
              'subjects omitted\n')
          #	add zero random effects for new levels in id
          newid <- !levels(id) %in% levels.obj
          if (sum(newid) > 0) {
            newre <- matrix(
              0,
              nrow = sum(newid),
              ncol = dim(ro)[2],
              dimnames = list(levels(id)[newid], dimnames(ro)[[2]])
            )
            start.$random <- rbind(start.$random, newre)
            cat(sum(newid), 'subjects added\n')
          }
        }
      }
      #	update fixed effects
      if (length(fo) > df + 1)
        fixed.extra <- names(fo)[(df + 2):length(fo)]
      else
        fixed.extra <- NULL
      # new arg xoffset
      if (!is.null(extras$xoffset)) {
        xoffset.t <- xoffset
        xoffset <- eval.parent(extras$xoffset)
        xoffset.t <- xoffset - xoffset.t
        x <- x - xoffset.t
        knots <- knots - xoffset.t
        bounds <- bounds - xoffset.t
      }
      # new arg knots
      if (!is.null(extras$knots)) {
        knots <- eval.parent(extras$knots) - xoffset
        df <- length(knots) + 1
         mcall$df <- NULL
     }
      # new arg df
      else if (!is.null(extras$df)) {
        df <- eval.parent(extras$df)
        mcall$knots <- NULL
        knots <- quantile(x, 1:(df - 1) / df)
      }
      # new arg bounds
      if (!is.null(extras$bounds)) {
        bounds <- eval.parent(extras$bounds)
        if (length(bounds) == 1)
          bounds <- range(x) + abs(bounds) * c(-1, 1) * diff(range(x))
        else
          bounds <- bounds - xoffset
      }
      #	get spline start values
      stype <- extras$stype %||% attr(object, 'stype') %||% 'ns'
      splinefun <- get(stype)
      # define missing b and c fixed effects
      fo[setdiff(letters[2:3], names(fo))] <- 0
      # refit spline to update coefficients
      spline.lm <-lm(predict(object, data, level = 0) ~
                       splinefun(((x - fo['b']) * exp(fo['c'])),
                                 knots = knots, Bound = bounds))
      start.$fixed <-
        c(coef(spline.lm)[c(2:(df + 1), 1)],
          fo[names(fo) %in% fixed.extra])
      # new arg bstart
      if (!is.null(extras$bstart) && !is.null(start.$fixed['b'])) {
        bstart <- eval.parent(extras$bstart)
        if (bstart == 'mean')
          bstart <- mean(x)
        else
          bstart <- bstart - xoffset
        start.$fixed['b'] <- bstart
      }
    }
    #	save start. object
    assign('start.', start., parent.frame())
    mcall$start <- quote(start.)
  }
  if (evaluate) {
    # if data stored in object$data then data not in mcall
    # add data to mcall and put copy of data in globalenv
    if (!'data' %in% names(mcall)) {
      mcall$data <- quote(.data.)
      assign('.data.', object$data, parent.frame())
    }
    eval.parent(mcall)
  }
  else
    mcall
}
