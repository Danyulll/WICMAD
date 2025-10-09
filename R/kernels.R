
# Kernels (with +Bias variants) ------------------------------------------

make_kernels <- function(add_bias_variants = TRUE) {
  k_sqexp <- function(t, l_scale) {
    D2 <- dist_rows(t)^2
    exp(-0.5 * D2 / (l_scale^2))
  }
  k_mat32 <- function(t, l_scale) {
    D  <- dist_rows(t); r <- D / l_scale; a <- sqrt(3) * r
    (1 + a) * exp(-a)
  }
  k_mat52 <- function(t, l_scale) {
    D  <- dist_rows(t); r <- D / l_scale; a <- sqrt(5) * r
    (1 + a + 5 * r^2 / 3) * exp(-a)
  }
  k_periodic <- function(t, l_scale, period) {
    D <- dist_rows(t)
    exp(- 2 * sin(base::pi * D / period)^2 / (l_scale^2))
  }

  base <- list(
    list(
      name="SE", fun=function(t, par) k_sqexp(t, par$l_scale),
      pnames=c("l_scale"),
      prior = function(par) stats::dgamma(par$l_scale, 2, 2, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 2, 2)),
      prop_sd = list(l_scale=0.20)
    ),
    list(
      name="Mat32", fun=function(t, par) k_mat32(t, par$l_scale),
      pnames=c("l_scale"),
      prior = function(par) stats::dgamma(par$l_scale, 2, 2, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 2, 2)),
      prop_sd = list(l_scale=0.20)
    ),
    list(
      name="Mat52", fun=function(t, par) k_mat52(t, par$l_scale),
      pnames=c("l_scale"),
      prior = function(par) stats::dgamma(par$l_scale, 2, 2, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 2, 2)),
      prop_sd = list(l_scale=0.20)
    ),
    list(
      name="Periodic", fun=function(t, par) k_periodic(t, par$l_scale, par$period),
      pnames=c("l_scale","period"),
      prior = function(par) stats::dgamma(par$l_scale, 3, 2, log=TRUE) + stats::dbeta(par$period, 5, 5, log=TRUE),
      pstar = function() list(l_scale = stats::rgamma(1, 3, 2), period = stats::rbeta(1, 5, 5)),
      prop_sd = list(l_scale=0.20, period=0.20)
    )
  )

  if (!add_bias_variants) return(base)

  k_bias <- function(t, s0) { P <- nloc(t); (s0^2) * matrix(1, P, P) }

  make_bias_variant <- function(kcfg) {
    nm <- paste0(kcfg$name, "+Bias")
    fun_bias <- function(t, par) {
      Kbase <- kcfg$fun(t, par)
      Kbase + k_bias(t, par$s0)
    }
    prior_bias <- function(par) {
      kcfg$prior(par) + dnorm(log(pmax(par$s0, 1e-9)), mean = -3, sd = 0.75, log = TRUE)
    }
    pstar_bias <- function() { th <- kcfg$pstar(); th$s0 <- exp(rnorm(1, -3, 0.75)); th }
    prop_sd_bias <- kcfg$prop_sd; prop_sd_bias$s0 <- 0.25
    list(
      name   = nm,
      fun    = fun_bias,
      pnames = c(kcfg$pnames, "s0"),
      prior  = prior_bias,
      pstar  = pstar_bias,
      prop_sd= prop_sd_bias
    )
  }

  bias_variants <- lapply(base, make_bias_variant)
  c(base, bias_variants)
}
