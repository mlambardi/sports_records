source("lib.R")

# Generates a sequence from an AR(1) with EVD marginals
# The last value of the sequence can be supplied
rargumbel <- function(n, m = 0, s = 1, r = 0, last) {
  se <- sqrt(1 - r^2)
  y <- rnorm(n)
  if (!missing(last)) {
    y[n] <- qnorm(EnvStats::pevd(last, location = m, scale = s))
  }
  for (i in (n-1):1) {
    y[i] <- r * y[i + 1] + se * y[i]
  }
  EnvStats::qevd(pnorm(y), location = m, scale = s)
}

# Makes a "model" class
argumbel <- R6Class(
  classname = "AR-Gumbel model",
  inherit = model,
  private = list(
    cdf = function(q, m, s, r, last, ...) {
      last <- qnorm(EnvStats::pevd(last, location = m, scale = s))
      qz <- qnorm(EnvStats::pevd(q, location = m, scale = s))
      pnorm(qz, r * last, sqrt(1 - r^2))
    },
    qf = function(p, m, s, r, last, ...) {
      last <- qnorm(EnvStats::pevd(last, location = m, scale = s))
      qz <- qnorm(p, r * last, sqrt(1 - r^2))
      EnvStats::qevd(pnorm(qz), location = m, scale = s)
    }
  ),
  public = list(
    est = function(trainingdata, indep = F) {
      y <- trainingdata
      n <- length(y)
      # marginal composite likelihood for mu, sigma
      ml <- as.list(EnvStats::eevd(y)$parameters)
      ml <- setNames(ml[c("location", "scale")], c("m", "s"))
      ml$r <- if (indep) {
        0
      } else {
        # semi-parametric estimation of Spearman copula
        2 * sin(pi / 6 * as.vector(cor(y[-1], y[-n], method = "spearman")))
      }
      ml$last <- tail(y, 1)
      ml$n <- n
      ml$indep <- indep
      ml
    },
    sim = function() {
      v <- private$par
      rargumbel(n = v$n, m = v$m, s = v$s, r = v$r, last = v$last)
    }
  )
)

# Expected quantile loss for a given quantile, under a Gaussian DGP
expected_quantile_loss_normal <- function(qpred, p, mean = 0, sd = 1) {
  plower <- pnorm(qpred, mean, sd)
  pupper <- 1 - plower
  elower <- qpred - truncnorm::etruncnorm(b = qpred, mean = mean, sd = sd)
  eupper <- truncnorm::etruncnorm(a = qpred, mean = mean, sd = sd) - qpred
  (1 - p) * plower * elower + p * pupper * eupper
}
