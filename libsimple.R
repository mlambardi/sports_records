require(dplyr, quietly = T)
require(R6, quietly = T)

# An R6 class in R that implements a general model
model <- R6Class(
  classname = "model",
  private = list(
    # par stores estimated parameters
    par = NULL,
    # opts stores additional options
    opts = NULL
  ),
  public = list(
    # getters
    get_par = function() private$par,
    get_opts = function() private$opts,
    # method used to instantiate when model$new is called
    # requires trainingdata and model options, or par as an alternative
    initialize = function(trainingdata, ..., par) {
      if (missing(par)) {
        self$train(trainingdata = trainingdata, ...)
      } else {
        private$opts <- list(...)
        private$par <- par
      }
    },
    # setter
    train = function(trainingdata, ...) {
      private$opts <- list(...)
      private$par <- self$est(trainingdata = trainingdata, ...)
    },
    # estimates the model parameter
    est = function(trainingdata, ...) {
      stop("must implement MLE")
    },
    # generates a synthetic dataset based on the estimated model
    sim = function() {
      stop("must implement DGP")
    },
    # generates makes nboot bootstrap estimates of the same model based on as many synthetic datasets
    estboot = function(nboot) {
      lapply(seq_len(nboot), function(foo) {
        myclone <- self$clone()
        args <- c(list(trainingdata = self$sim()), private$opts)
        do.call(myclone$train, args = args)
        myclone
      })
    },
    # estimative cumulative distribution function (CDF)
    p = function(q) {
      do.call(private$cdf, args = purrr::list_modify(private$par, q = q))
    },
    # estimative quantile function (QF)
    q = function(p) {
      do.call(private$qf, args = purrr::list_modify(private$par, p = p))
    },
    # Treats itself as a bootstrap estimate, to be combined with MLE
    # When averaging over multiple bootstraps, one gets the probability-calibrated QF
    # to be inverted to get the probability-calibrated CDF (see the "calibration" R6 class)
    qpq = function(p, ml) {
      ml$q(self$p(ml$q(p)))
    },
    # Treats itself as a bootstrap estimate, to be combined with MLE
    # When averaging over multiple bootstraps, one gets the quantile-calibrated CDF
    # to be inverted to get the quantile-calibrated QF (see the "calibration" R6 class)
    pqp = function(q, ml) {
      ml$p(self$q(ml$p(q)))
    }
  )
)

# An R6 class in R that implements the estimative and the calibrated predictions
# requires a model (see the dedicated R6 class above)
calibration <- R6Class(
  classname = "calibration",
  private = list(
    # base estimate
    ml = NULL,
    # derived bootstrap estimates
    bts = NULL
  ),
  public = list(
    # getters
    get_ml = function() private$ml,
    get_bts = function() private$bts,
    # returns the base and bootstrap estimates as a dataframe, for inspection only
    get_par = function() {
      bind_rows(
        as.data.frame(private$ml$get_par()) |>
          mutate(isMLE = T, bootstrap = as.integer(NA)),
        lapply(private$bts, \(b) as.data.frame(b$get_par())) |>
          bind_rows() |>
          mutate(isMLE = F, bootstrap = as.integer(1:n()))
      )
    },
    # must provide a estimated model "est" and tell how many "nboot" bootstraps are wanted
    # As an alternative, can provide an R6 class extending the stub class "model" provided above
    initialize = function(est, nboot, par, mod) {
      if (missing(par)) {
        private$ml <- est
        private$bts <- est$estboot(nboot = nboot)
      } else {
        private$ml <- mod$new(par = as.list(filter(par, isMLE)))
        private$bts <- filter(par, !isMLE) %>%
          split(., seq_len(nrow(.))) |>
          lapply(\(p, mod) mod$new(par = as.list(p)), mod = mod)
      }
    },
    # estimative prediction, only uses the base estimate
    qest = function(p) private$ml$q(p),
    pest = function(q) private$ml$p(q),
    # probability-calibrated predictions
    # qCP implements the probability-calibrated QF, to be inverted to obtain calibrated probabilities, see pCP later on
    qCP = function(p) {
      # this averages over the bootstraps' qpq quantities
      purrr::reduce(
        .x = private$bts,
        .f = function(prev, bt, p, ml) {
          prev + bt$qpq(p, ml = ml)
        },
        p = p,
        ml = private$ml,
        .init = numeric(length(p))
      ) /
        length(private$bts)
    },
    # pCQ implements the quantile-calibrated CDF, to be inverted to obtain calibrated quantiles, see qCQ later on
    pCQ = function(q) {
      # this averages over the bootstraps' pqp quantities
      purrr::reduce(
        .x = private$bts,
        .f = function(prev, bt, q, ml) {
          prev + bt$pqp(q, ml = ml)
        },
        q = q,
        ml = private$ml,
        .init = numeric(length(q))
      ) /
        length(private$bts)
    },
    # computes bounds for qCQ, to ease inversion
    qCQ_limits = function(p) {
      # estimative QF is included in the bounds, so one can fall back to estimative in line of principle
      qest <- private$ml$q(p)
      purrr::reduce(
        .x = private$bts,
        .f = function(prev, bt, p, ml) {
          qs <- bt$qpq(p, ml = ml)
          mutate(prev, qinf = pmin(qinf, qs), qsup = pmax(qsup, qs))
        },
        p = p,
        ml = private$ml,
        # this init ensures estimative can be fallen back to
        .init = data.frame(qinf = qest, qsup = qest)
      )
    },
    # computes bounds for pCP, to ease inversion
    pCP_limits = function(q) {
      # estimative CDF is included in the bounds, so one can fall back to estimative in line of principle
      pest <- private$ml$p(q)
      purrr::reduce(
        .x = private$bts,
        .f = function(prev, bt, q, ml) {
          ps <- bt$pqp(q, ml = ml)
          mutate(prev, pinf = pmin(pinf, ps), psup = pmax(psup, ps))
        },
        q = q,
        ml = private$ml,
        # this init ensures estimative can be fallen back to
        .init = data.frame(pinf = pest, psup = pest)
      )
    },
    # pCP simpler than qCQ, so addressing pCP first
    # calibrated probabilities, inverse of qCP, obtained via bisection starting from pCP_limits
    pCP = function(q, tol = 1e-7) {
      plims <- self$pCP_limits(q)
      pcand <- qcand <- rep(as.numeric(NA), length(q))
      for (i in seq_len(max(1, ceiling(log(tol, 0.5))))) {
        pcand[] <- 0.5 * (plims$pinf + plims$psup)
        qcand[] <- self$qCP(pcand)
        plims$pinf[] <- ifelse(qcand <= q, pcand, plims$pinf)
        plims$psup[] <- ifelse(qcand >= q, pcand, plims$psup)
      }
      rowMeans(plims)
    },
    # calibrated quantiles, inverse of pCQ, obtained via bisection starting from pCQ_limits
    qCQ = function(p, tol = 1e-7, maxit = 100) {
      lims <- self$qCQ_limits(p) |>
        mutate(
          # handles numeric infinites
          infinity = !is.finite(qinf) | !is.finite(qsup),
          warn = is.finite(qinf) != is.finite(qsup),
          qinf = ifelse(warn & !is.finite(qsup), qsup, qinf),
          qsup = ifelse(warn & !is.finite(qinf), qinf, qsup),
          pinf = self$pCQ(qinf),
          psup = self$pCQ(qsup),
          needsupdate = T,
          qcand = as.numeric(NA),
          pcand = as.numeric(NA),
          iters = 0
        )
      if (any(lims$warn)) {
        warning("entering a grey zone")
      }
      # ensuring max error on quantile is proportional to some scale parameter of the distribution
      # for invariance purposes
      tol <- abs(tol)
      tolq <- tol * abs(diff(private$ml$q(c(0.25, 0.75))))
      maxit <- ceiling(max(maxit, 1)) # at least one iteration, to make non-NA q's and p's
      for (it in seq_len(maxit)) {
        if (!any(lims$needsupdate)) break
        lims <- lims |>
          mutate(
            iters = iters + needsupdate, # +0 or +1
            qcand = ifelse(needsupdate, 0.5 * (qinf + qsup), qcand),
            # repeated evaluations not that costly, to be fixed in the future
            pcand = self$pCQ(qcand),
            qinf = ifelse(needsupdate & (pcand <= p), qcand, qinf),
            qsup = ifelse(needsupdate & (pcand >= p), qcand, qsup),
            needsupdate = !infinity &
              (abs(qsup - qinf) > tolq) &
              (abs(psup - pinf) > tol)
          )
      }
      lims <- lims |>
        mutate(
          warn = warn | (needsupdate & (it >= maxit)),
          qcand = 0.5 * (qinf + qsup)
        )
      q <- lims$qcand
      attr(q, "warn") <- lims$warn
      attr(q, "iters") <- lims$iters
      q
    }
  )
)

# for the paper's analyses:
source("libsimpleARGumbel.R")
