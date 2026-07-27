###############################################################################
## Probability of Success (PoS) for a Phase 3 study, based on Phase 2 data
## -----------------------------------------------------------------------------
## Method: FREQUENTIST CONDITIONAL POWER.
##   The Phase 2 estimate of the treatment effect is treated as the assumed
##   "true" effect for planning. PoS = the statistical power of the planned
##   Phase 3 test evaluated at that Phase 2 effect.
##
## For every endpoint the script returns TWO estimates of PoS:
##   (1) analytic  -- a closed-form power formula (fast), and
##   (2) simulated -- a Monte Carlo estimate that generates Phase 3 data sets
##                    and runs the actual test, as a cross-check.
##
## Endpoints supported:
##   - continuous       (two-sample mean difference; t-test)
##   - binary           (two proportions; stratified Miettinen-Nurminen)
##   - time-to-event    (log hazard ratio; log-rank / Schoenfeld)
##
## TIME-TO-EVENT DESIGN MODEL (updated, gsDesign2-style structure):
##   * PROPORTIONAL HAZARDS: a single scalar hr (treat vs ctrl).
##   * PIECEWISE-CONSTANT ENROLLMENT rates:
##       enroll_rate = data.frame(duration, rate)
##   * PIECEWISE-EXPONENTIAL failure and dropout rates (control arm):
##       fail_rate = data.frame(duration, fail_rate, dropout_rate)
##     The treatment arm hazard is fail_rate * hr per period (PH);
##     dropout is shared by both arms. The last period is extended
##     indefinitely (as in gsDesign2).
##   Expected events are computed in closed form (piecewise-exponential
##   integrals; no numeric integration), following the architecture of
##   gsDesign2: a single-arm core (.expected_event_1arm, cf.
##   gsDesign2::expected_event) plus a two-arm wrapper (expected_events,
##   cf. gsDesign2::pw_info) that splits enrollment by allocation and sets
##   the experimental hazard to fail_rate * hr.
##   The LEGACY interface (n_treat, n_ctrl, lambda_ctrl, hr, follow_time,
##   enroll_time) still works and now also accepts a scalar dropout_rate.
##
## Each endpoint has:
##   - a summary extractor  ph2_summary_*()   (raw Phase 2 data -> summary stats)
##   - a PoS function        pos_*()          (summary stats + Phase 3 design -> PoS)
## and a single dispatcher   prob_success()   ties them together.
##
## Only base R is required for continuous and binary. Time-to-event simulation
## uses the 'survival' package (analytic TTE PoS needs no extra package).
###############################################################################

suppressWarnings(suppressMessages({
  has_survival <- requireNamespace("survival", quietly = TRUE)
}))

## ---------------------------------------------------------------------------
## Small internal helpers
## ---------------------------------------------------------------------------

## Critical value for a one- or two-sided test at level alpha.
.crit <- function(alpha, alternative = c("two.sided", "one.sided")) {
  alternative <- match.arg(alternative)
  if (alternative == "two.sided") qnorm(1 - alpha / 2) else qnorm(1 - alpha)
}

## Normal-approximation power given a standardized effect (effect / SE under the
## alternative), the null SE, and the design. Kept generic so all three
## endpoints share one power engine.
##   z_effect = effect / se_alt         (signal-to-noise under alternative)
##   se_ratio = se_null / se_alt        (accounts for null vs alt variance, e.g. binary)
.normal_power <- function(z_effect, alpha, alternative, se_ratio = 1) {
  zc <- .crit(alpha, alternative)
  upper <- pnorm(z_effect - zc * se_ratio)          # reject in the effect direction
  if (alternative == "two.sided") {
    lower <- pnorm(-z_effect - zc * se_ratio)        # reject in the opposite tail
    upper + lower
  } else {
    upper
  }
}

## ===========================================================================
## 1. CONTINUOUS ENDPOINT  (two-sample mean difference)
## ===========================================================================

## Turn raw Phase 2 arm vectors into the summary stats PoS needs.
ph2_summary_continuous <- function(y_treat, y_ctrl) {
  n_t <- length(y_treat); n_c <- length(y_ctrl)
  sd_pooled <- sqrt(((n_t - 1) * var(y_treat) + (n_c - 1) * var(y_ctrl)) /
                      (n_t + n_c - 2))
  list(delta = mean(y_treat) - mean(y_ctrl),  # observed mean difference
       sd    = sd_pooled,                      # pooled within-arm SD
       n_treat = n_t, n_ctrl = n_c)
}

## Phase 3 PoS for a continuous endpoint.
##   delta   : assumed effect (Phase 2 mean difference, treat - ctrl)
##   sd      : assumed common within-arm SD (from Phase 2)
##   n_treat, n_ctrl : PLANNED Phase 3 per-arm sample sizes
pos_continuous <- function(delta, sd, n_treat, n_ctrl = n_treat,
                           alpha = 0.05,
                           alternative = c("two.sided", "one.sided"),
                           nsim = 10000, seed = NULL) {
  alternative <- match.arg(alternative)
  se_alt <- sd * sqrt(1 / n_treat + 1 / n_ctrl)   # SE of the mean difference in P3
  df     <- n_treat + n_ctrl - 2
  
  ## Analytic: exact under the t-test via the non-central t distribution.
  ncp <- delta / se_alt
  if (alternative == "two.sided") {
    tc <- qt(1 - alpha / 2, df)
    pos_analytic <- pt(-tc, df, ncp) + (1 - pt(tc, df, ncp))
  } else {
    tc <- qt(1 - alpha, df)
    pos_analytic <- 1 - pt(tc, df, ncp)
  }
  
  ## Simulation: generate P3 arms and run the real t-test.
  pos_sim <- NA_real_
  if (!is.null(nsim) && nsim > 0) {
    if (!is.null(seed)) set.seed(seed)
    hits <- 0L
    for (i in seq_len(nsim)) {
      yt <- rnorm(n_treat, mean = delta, sd = sd)   # ctrl mean fixed at 0 (only diff matters)
      yc <- rnorm(n_ctrl,  mean = 0,     sd = sd)
      p  <- t.test(yt, yc, var.equal = TRUE,
                   alternative = if (alternative == "one.sided") "greater" else "two.sided")$p.value
      if (p < alpha) hits <- hits + 1L
    }
    pos_sim <- hits / nsim
  }
  
  structure(list(
    endpoint = "continuous", method = "conditional power",
    inputs = list(delta = delta, sd = sd, n_treat = n_treat, n_ctrl = n_ctrl,
                  alpha = alpha, alternative = alternative),
    effect_size_d = delta / sd, se_phase3 = se_alt,
    pos_analytic = pos_analytic, pos_simulated = pos_sim, nsim = nsim
  ), class = "pos_result")
}

## ===========================================================================
## 2. BINARY ENDPOINT  (risk difference, STRATIFIED MIETTINEN-NURMINEN test)
## ---------------------------------------------------------------------------
## The Phase 3 analysis is a stratified Miettinen-Nurminen (M&N) score test of
## the common risk difference  RD = p_treat - p_ctrl  against a null value
## 'margin' (margin = 0 for superiority; a negative margin for non-inferiority
## when the event is a "good" outcome).
##
## Per stratum k the M&N score statistic uses the RESTRICTED MLE of the two
## proportions under the constraint p1 - p2 = margin, and a variance with the
## small-sample factor N_k/(N_k - 1). Strata are combined with weights w_k
## (Mantel-Haenszel by default):
##
##     Z = sum_k w_k (phat1_k - phat2_k - margin)
##         ------------------------------------------
##            sqrt( sum_k w_k^2 * Vtilde_k )
##
## Stratification is OPTIONAL: pass scalar rates for an unstratified test, or
## length-K vectors of per-stratum rates for a K-stratum test.
## ===========================================================================

## Restricted (constrained) MLE of (p1, p2) under H0: p1 - p2 = delta.
## Vectorised over strata. For delta = 0 the solution is the pooled proportion
## in closed form; otherwise each stratum is solved by 1-D optimisation of the
## constrained binomial log-likelihood (exact, no fragile cubic).
.mn_rmle <- function(x1, n1, x2, n2, delta) {
  if (abs(delta) < 1e-12) {                    # superiority: pooled proportion
    p <- (x1 + x2) / (n1 + n2)
    return(list(p1 = p, p2 = p))
  }
  p2 <- mapply(function(a, na, b, nb) {
    lo <- max(0, -delta) + 1e-9
    hi <- min(1, 1 - delta) - 1e-9
    if (lo >= hi) return((lo + hi) / 2)
    negll <- function(q) {                     # -log-lik, free parameter q = p2
      p1 <- q + delta
      -(a * log(p1) + (na - a) * log(1 - p1) +
          b * log(q)  + (nb - b) * log(1 - q))
    }
    optimize(negll, c(lo, hi), tol = 1e-10)$minimum
  }, x1, n1, x2, n2)
  list(p1 = p2 + delta, p2 = p2)
}

## M&N null variance per stratum at the restricted MLE (with N/(N-1) factor).
.mn_var0 <- function(rm, n1, n2) {
  Nk <- n1 + n2
  (rm$p1 * (1 - rm$p1) / n1 + rm$p2 * (1 - rm$p2) / n2) * (Nk / (Nk - 1))
}

## Split a per-arm sample size across strata according to prevalences 'prop',
## returning integer counts that sum exactly to n_arm.
strat_sizes <- function(n_arm, prop) {
  prop <- prop / sum(prop)
  raw  <- prop * n_arm
  k    <- floor(raw)
  rem  <- n_arm - sum(k)
  if (rem > 0) {
    ord <- order(raw - k, decreasing = TRUE)
    k[ord[seq_len(rem)]] <- k[ord[seq_len(rem)]] + 1L
  }
  as.integer(k)
}

## Raw Phase 2 data (0/1 vectors, OR event counts) -> summary proportions.
ph2_summary_binary <- function(y_treat = NULL, y_ctrl = NULL,
                               x_treat = NULL, n_treat = NULL,
                               x_ctrl = NULL, n_ctrl = NULL) {
  if (!is.null(y_treat)) {          # 0/1 vectors supplied
    list(p_treat = mean(y_treat), p_ctrl = mean(y_ctrl),
         n_treat = length(y_treat), n_ctrl = length(y_ctrl))
  } else {                          # event counts supplied
    list(p_treat = x_treat / n_treat, p_ctrl = x_ctrl / n_ctrl,
         n_treat = n_treat, n_ctrl = n_ctrl)
  }
}

## Per-STRATUM Phase 2 summary from raw data.
##   response : 0/1 outcome
##   group    : arm label; 'treat_label' identifies the treatment arm
##   stratum  : the stratification factor. If it is continuous, pass it through
##              bin_continuous() first; ordinal/nominal factors are used as-is.
## Returns per-stratum treatment/control rates and the observed stratum
## prevalence (a natural default for how to split the Phase 3 sample).
ph2_summary_binary_strat <- function(response, group, stratum, treat_label) {
  stratum <- as.factor(stratum)
  lv <- levels(stratum)
  out <- lapply(lv, function(s) {
    idx <- stratum == s
    is_t <- idx & group == treat_label
    is_c <- idx & group != treat_label
    c(p_treat = mean(response[is_t]), p_ctrl = mean(response[is_c]),
      prevalence = mean(idx))
  })
  m <- do.call(rbind, out)
  list(stratum = lv, p_treat = m[, "p_treat"], p_ctrl = m[, "p_ctrl"],
       strata_prop = m[, "prevalence"])
}

## Bin a continuous stratifier into K categories (equal-frequency quantiles by
## default) so it can be used as a stratum factor.
bin_continuous <- function(x, k = 3, method = c("quantile", "equalwidth")) {
  method <- match.arg(method)
  br <- if (method == "quantile")
    unique(quantile(x, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE))
  else seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = k + 1)
  cut(x, breaks = br, include.lowest = TRUE)
}

## Phase 3 PoS for a binary endpoint via stratified Miettinen-Nurminen.
##   p_treat, p_ctrl : assumed response rates. Scalars => unstratified;
##                     length-K vectors => per-stratum rates (K strata).
##   n_treat, n_ctrl : TOTAL planned per-arm sample size (scalars). Split across
##                     strata using strata_prop.
##   strata_prop     : stratum prevalences (length K, need not sum to 1).
##                     NULL => equal split when K > 1. Use the Phase 2 observed
##                     prevalences (ph2_summary_binary_strat) or population rates.
##   margin          : null risk difference (0 = superiority; e.g. -0.10 for a
##                     10-point non-inferiority margin on a good outcome).
##   weights         : stratum weighting, "MH" (Mantel-Haenszel) or "inv_var".
pos_binary <- function(p_treat, p_ctrl, n_treat, n_ctrl = n_treat,
                       strata_prop = NULL, margin = 0,
                       weights = c("MH", "inv_var"),
                       alpha = 0.05,
                       alternative = c("two.sided", "one.sided"),
                       nsim = 10000, seed = NULL) {
  alternative <- match.arg(alternative)
  weights     <- match.arg(weights)
  K <- length(p_treat)
  stopifnot(length(p_ctrl) == K)
  
  ## Per-stratum sample sizes from prevalences.
  if (K == 1) {
    n1 <- n_treat; n2 <- n_ctrl
  } else {
    if (is.null(strata_prop)) strata_prop <- rep(1 / K, K)
    n1 <- strat_sizes(n_treat, strata_prop)
    n2 <- strat_sizes(n_ctrl,  strata_prop)
  }
  
  RD    <- p_treat - p_ctrl
  V_alt <- p_treat * (1 - p_treat) / n1 + p_ctrl * (1 - p_ctrl) / n2   # unconstrained var
  w <- if (weights == "MH") n1 * n2 / (n1 + n2) else 1 / V_alt         # stratum weights
  zc <- .crit(alpha, alternative)
  
  ## --- Analytic (conditional power) -----------------------------------------
  ## Treat the assumed true rates as "expected" data to get the null variance.
  rm0  <- .mn_rmle(n1 * p_treat, n1, n2 * p_ctrl, n2, margin)
  V0   <- .mn_var0(rm0, n1, n2)
  effect  <- sum(w * (RD - margin))                 # E[numerator] under the truth
  se_null <- sqrt(sum(w^2 * V0))
  se_alt  <- sqrt(sum(w^2 * V_alt))
  pos_analytic <- .normal_power(z_effect = effect / se_alt, alpha = alpha,
                                alternative = alternative,
                                se_ratio = se_null / se_alt)
  
  ## --- Simulation: run the actual stratified M&N test -----------------------
  pos_sim <- NA_real_
  if (!is.null(nsim) && nsim > 0) {
    if (!is.null(seed)) set.seed(seed)
    hits <- 0L
    for (i in seq_len(nsim)) {
      x1 <- rbinom(K, n1, p_treat)
      x2 <- rbinom(K, n2, p_ctrl)
      rm <- .mn_rmle(x1, n1, x2, n2, margin)
      Vt <- .mn_var0(rm, n1, n2)
      num <- sum(w * (x1 / n1 - x2 / n2 - margin))
      den <- sqrt(sum(w^2 * Vt))
      z   <- if (den > 0) num / den else 0
      reject <- if (alternative == "two.sided") abs(z) > zc else z > zc
      if (isTRUE(reject)) hits <- hits + 1L
    }
    pos_sim <- hits / nsim
  }
  
  structure(list(
    endpoint = "binary",
    method = sprintf("conditional power / stratified Miettinen-Nurminen (%s weights)",
                     weights),
    inputs = list(n_treat = n_treat, n_ctrl = n_ctrl, n_strata = K,
                  margin = margin, alpha = alpha, alternative = alternative),
    per_stratum = data.frame(p_treat = p_treat, p_ctrl = p_ctrl,
                             n_treat = n1, n_ctrl = n2, RD = RD, weight = w),
    risk_difference = sum(w * RD) / sum(w),   # weighted overall RD
    se_phase3 = se_alt / sum(w),              # SE of the weighted RD estimate (RD scale)
    pos_analytic = pos_analytic, pos_simulated = pos_sim, nsim = nsim
  ), class = "pos_result")
}

## ===========================================================================
## 3. TIME-TO-EVENT ENDPOINT  (hazard ratio, log-rank test)
## ---------------------------------------------------------------------------
## Design model (gsDesign2-style, PROPORTIONAL HAZARDS):
##   enroll_rate : data.frame(duration, rate)      -- piecewise-constant
##                 enrollment rates; total N = sum(duration * rate).
##   fail_rate   : data.frame(duration, fail_rate, dropout_rate)
##                 -- piecewise-exponential CONTROL-arm failure hazard and
##                 dropout hazard. Treatment hazard = fail_rate * hr (PH);
##                 dropout is common to both arms. The LAST period is
##                 extended indefinitely.
##   total_duration : calendar time from start of enrollment to data cutoff.
## Helper constructors (compatible with gsDesign2::define_*):
##   define_enroll_rate(duration, rate)
##   define_fail_rate(duration, fail_rate, dropout_rate)
## ===========================================================================

define_enroll_rate <- function(duration, rate) {
  stopifnot(length(duration) == length(rate), all(duration > 0), all(rate >= 0))
  data.frame(duration = duration, rate = rate)
}

define_fail_rate <- function(duration, fail_rate, dropout_rate = 0) {
  n <- length(duration)
  stopifnot(length(fail_rate) %in% c(1L, n), length(dropout_rate) %in% c(1L, n),
            all(duration > 0), all(fail_rate >= 0), all(dropout_rate >= 0))
  data.frame(duration = duration,
             fail_rate = rep_len(fail_rate, n),
             dropout_rate = rep_len(dropout_rate, n))
}

## Raw Phase 2 survival data -> estimated hazard ratio (Cox).
## group: coded so that the reference is control; hr is treat vs ctrl.
ph2_summary_tte <- function(time, status, group, treat_label = NULL) {
  if (!has_survival) stop("The 'survival' package is required for ph2_summary_tte().")
  group <- as.factor(group)
  if (!is.null(treat_label)) group <- relevel(group, ref = setdiff(levels(group), treat_label)[1])
  fit <- survival::coxph(survival::Surv(time, status) ~ group)
  list(hr = unname(exp(coef(fit))[1]),
       log_hr = unname(coef(fit)[1]),
       se_log_hr = sqrt(diag(vcov(fit)))[1],
       n_events = sum(status))
}

## --- Piecewise-exponential internals ---------------------------------------

## Precompute, per failure period j (boundaries t0_j from cumsum(duration)):
##   theta_j = fail + dropout (all-cause hazard), S_j = P(no event & no dropout
##   by t0_j), and pev_j = P(event by t0_j). Last period extended to Inf.
.pw_precompute <- function(fail_rate) {
  lam   <- fail_rate$fail_rate
  eta   <- fail_rate$dropout_rate
  theta <- lam + eta
  t0    <- c(0, cumsum(fail_rate$duration))
  t0[length(t0)] <- Inf                       # last period extends indefinitely
  J <- length(lam)
  S <- pev <- numeric(J + 1)
  S[1] <- 1; pev[1] <- 0
  for (j in seq_len(J)) {
    d <- t0[j + 1] - t0[j]
    Snext <- if (is.finite(d)) S[j] * exp(-theta[j] * d) else 0
    S[j + 1]   <- Snext
    pev[j + 1] <- pev[j] +
      if (theta[j] > 0) (lam[j] / theta[j]) * (S[j] - Snext) else 0
  }
  list(t0 = t0, lam = lam, eta = eta, theta = theta, S = S, pev = pev, J = J)
}

## P(event by follow-up time tau) for one subject (piecewise exp failure +
## dropout). Vectorised in tau.
.p_event_pw <- function(fail_rate, tau) {
  pw <- .pw_precompute(fail_rate)
  vapply(tau, function(x) {
    if (x <= 0) return(0)
    j <- findInterval(x, pw$t0, rightmost.closed = FALSE)  # period containing x
    j <- min(j, pw$J)
    if (pw$theta[j] > 0)
      pw$pev[j] + (pw$lam[j] / pw$theta[j]) * pw$S[j] *
      (1 - exp(-pw$theta[j] * (x - pw$t0[j])))
    else pw$pev[j]
  }, numeric(1))
}

## Single-arm expected number of events at calendar time 'total_duration'
## (cf. gsDesign2::expected_event). Closed form: the enrollment window is
## split at enrollment-rate change points and at calendar times where the
## follow-up tau = total_duration - u crosses a failure-period boundary;
## within each sub-interval P(event | entry u) = C + D * exp(-theta * tau)
## integrates analytically.
.expected_event_1arm <- function(enroll_rate, fail_rate, total_duration) {
  T <- total_duration
  pw <- .pw_precompute(fail_rate)
  e_bounds <- cumsum(enroll_rate$duration)
  E <- min(e_bounds[length(e_bounds)], T)     # entries after T contribute 0
  if (E <= 0) return(0)
  
  ## breakpoints in entry time u
  br <- c(0, e_bounds[e_bounds < E], E)
  f_break <- T - pw$t0[is.finite(pw$t0) & pw$t0 > 0]   # u where tau hits a period boundary
  br <- sort(unique(c(br, f_break[f_break > 0 & f_break < E])))
  
  total <- 0
  for (k in seq_len(length(br) - 1)) {
    u1 <- br[k]; u2 <- br[k + 1]
    um <- (u1 + u2) / 2
    r  <- enroll_rate$rate[findInterval(um, c(0, e_bounds), rightmost.closed = FALSE)]
    tau_m <- T - um
    j <- min(findInterval(tau_m, pw$t0), pw$J)         # failure period for this slice
    if (pw$theta[j] > 0) {
      a <- (pw$lam[j] / pw$theta[j]) * pw$S[j]
      ## P(event | tau) = pev_j + a * (1 - exp(-theta_j * (tau - t0_j)))
      ## integrate over u in [u1, u2] with tau = T - u  (exponent kept <= 0)
      th <- pw$theta[j]
      total <- total + r * (u2 - u1) * (pw$pev[j] + a) -
        r * a / th * (exp(th * (pw$t0[j] - T + u2)) - exp(th * (pw$t0[j] - T + u1)))
    } else {
      total <- total + r * (u2 - u1) * pw$pev[j]
    }
  }
  total
}

## Expected number of events, BOTH ARMS, under proportional hazards
## (cf. the two-arm assembly in gsDesign2::pw_info: control uses fail_rate,
## experimental uses fail_rate * hr, dropout shared).
##
## NEW interface (piecewise):
##   expected_events(hr = 0.7, enroll_rate = define_enroll_rate(...),
##                   fail_rate = define_fail_rate(...),
##                   total_duration = 30, ratio = 1)
##   ratio = experimental : control randomization ratio; enrollment is split
##   ratio/(1+ratio) vs 1/(1+ratio).
##
## LEGACY interface (kept for backward compatibility; positional args
## unchanged): expected_events(n_treat, n_ctrl, lambda_ctrl, hr,
##   follow_time, enroll_time = 0, dropout_rate = 0)
##   -> single exponential failure rate, uniform accrual over enroll_time
##   (enroll_time = 0 means everyone enters at time 0), analysis at
##   enroll_time + follow_time. dropout_rate = 0 reproduces old results.
expected_events <- function(n_treat = NULL, n_ctrl = NULL, lambda_ctrl = NULL,
                            hr, follow_time = NULL, enroll_time = 0,
                            dropout_rate = 0,
                            enroll_rate = NULL, fail_rate = NULL,
                            total_duration = NULL, ratio = 1) {
  if (!is.null(enroll_rate) && !is.null(fail_rate)) {
    ## ---- new piecewise path ----
    if (is.null(total_duration)) stop("Provide total_duration with enroll_rate/fail_rate.")
    q_e <- ratio / (1 + ratio); q_c <- 1 - q_e
    enroll_c <- transform(enroll_rate, rate = rate * q_c)
    enroll_e <- transform(enroll_rate, rate = rate * q_e)
    fail_c <- fail_rate
    fail_e <- transform(fail_rate, fail_rate = fail_rate * hr)
    return(.expected_event_1arm(enroll_c, fail_c, total_duration) +
             .expected_event_1arm(enroll_e, fail_e, total_duration))
  }
  ## ---- legacy path: build single-period structures ----
  if (is.null(n_treat) || is.null(lambda_ctrl) || is.null(follow_time))
    stop("Provide either (enroll_rate, fail_rate, total_duration) or the legacy arguments.")
  if (is.null(n_ctrl)) n_ctrl <- n_treat
  fr_c <- define_fail_rate(duration = 1e6, fail_rate = lambda_ctrl,
                           dropout_rate = dropout_rate)
  fr_e <- define_fail_rate(duration = 1e6, fail_rate = lambda_ctrl * hr,
                           dropout_rate = dropout_rate)
  total <- enroll_time + follow_time
  if (enroll_time <= 0) {
    ## instantaneous enrollment: no accrual averaging
    return(n_treat * .p_event_pw(fr_e, total) + n_ctrl * .p_event_pw(fr_c, total))
  }
  er_c <- define_enroll_rate(duration = enroll_time, rate = n_ctrl / enroll_time)
  er_e <- define_enroll_rate(duration = enroll_time, rate = n_treat / enroll_time)
  .expected_event_1arm(er_c, fr_c, total) + .expected_event_1arm(er_e, fr_e, total)
}

## --- Random sampling from the piecewise model (for simulation) -------------

## Piecewise-exponential event times via inverse transform of the cumulative
## hazard (rate table: duration, rate; last period extended to Inf).
.rpwexp <- function(n, duration, rate) {
  t0 <- c(0, cumsum(duration)); t0[length(t0)] <- Inf
  J  <- length(rate)
  H0 <- c(0, cumsum(rate * duration))                 # cum hazard at boundaries
  H0[length(H0)] <- Inf
  u <- rexp(n)                                        # target cumulative hazards
  out <- rep(Inf, n)
  for (j in seq_len(J)) {
    idx <- u >= H0[j] & u < H0[j + 1]
    if (any(idx)) {
      if (rate[j] > 0) out[idx] <- t0[j] + (u[idx] - H0[j]) / rate[j]
      ## rate 0 in last period with u beyond H0[j]: stays Inf (correct)
    }
  }
  out
}

## Entry (enrollment) times from the piecewise-constant enrollment intensity.
.renroll <- function(n, enroll_rate) {
  w <- enroll_rate$duration * enroll_rate$rate
  if (sum(w) <= 0) return(rep(0, n))
  t0 <- c(0, cumsum(enroll_rate$duration))
  piece <- sample.int(length(w), n, replace = TRUE, prob = w)
  runif(n, t0[piece], t0[piece + 1])
}

## Phase 3 PoS for a time-to-event endpoint (log-rank; Schoenfeld analytic).
## Three ways to call it:
##   (A) Analytic only          -> hr + n_events (target Phase 3 event count).
##   (B) Legacy design          -> hr, per-arm n, lambda_ctrl (or median_ctrl),
##                                 follow_time, optional enroll_time and
##                                 dropout_rate. Exponential failure, uniform
##                                 accrual.
##   (C) Piecewise design (new) -> hr, enroll_rate, fail_rate, total_duration,
##                                 ratio. Piecewise enrollment + piecewise
##                                 exponential failure/dropout, PH via scalar hr.
## alloc = proportion allocated to treatment (used only in path A; paths B/C
## derive it from the design).
pos_tte <- function(hr,
                    n_events = NULL,
                    n_treat = NULL, n_ctrl = NULL,
                    lambda_ctrl = NULL, median_ctrl = NULL,
                    follow_time = NULL, enroll_time = 0,
                    dropout_rate = 0,
                    enroll_rate = NULL, fail_rate = NULL,
                    total_duration = NULL, ratio = 1,
                    alloc = 0.5,
                    alpha = 0.05,
                    alternative = c("two.sided", "one.sided"),
                    nsim = 5000, seed = NULL) {
  alternative <- match.arg(alternative)
  loghr <- log(hr)
  if (is.null(lambda_ctrl) && !is.null(median_ctrl)) lambda_ctrl <- log(2) / median_ctrl
  
  ## ---- normalize the design to the piecewise representation ----
  piecewise_in <- !is.null(enroll_rate) && !is.null(fail_rate)
  legacy_in    <- !is.null(n_treat) && !is.null(lambda_ctrl) && !is.null(follow_time)
  do_sim <- piecewise_in || legacy_in
  
  design <- NULL
  if (piecewise_in) {
    if (is.null(total_duration)) stop("Provide total_duration with enroll_rate/fail_rate.")
    q_e <- ratio / (1 + ratio)
    N   <- sum(enroll_rate$duration * enroll_rate$rate)
    n_e <- round(N * q_e); n_c <- round(N) - n_e
    design <- list(
      enroll_c = transform(enroll_rate, rate = rate * (1 - q_e)),
      enroll_e = transform(enroll_rate, rate = rate * q_e),
      fail_c   = fail_rate,
      fail_e   = transform(fail_rate, fail_rate = fail_rate * hr),
      n_treat  = n_e, n_ctrl = n_c, total = total_duration)
    alloc <- q_e
  } else if (legacy_in) {
    if (is.null(n_ctrl)) n_ctrl <- n_treat
    total <- enroll_time + follow_time
    er <- if (enroll_time > 0)
      define_enroll_rate(duration = enroll_time, rate = 1 / enroll_time)  # shape only
    else NULL
    design <- list(
      enroll_c = if (!is.null(er)) transform(er, rate = rate * n_ctrl) else NULL,
      enroll_e = if (!is.null(er)) transform(er, rate = rate * n_treat) else NULL,
      fail_c = define_fail_rate(1e6, lambda_ctrl,      dropout_rate),
      fail_e = define_fail_rate(1e6, lambda_ctrl * hr, dropout_rate),
      n_treat = n_treat, n_ctrl = n_ctrl, total = total)
    alloc <- n_treat / (n_treat + n_ctrl)
  }
  
  ## ---- expected events (analytic, closed form) ----
  if (is.null(n_events) && do_sim) {
    ev1 <- function(enr, fr, n) {
      if (is.null(enr)) n * .p_event_pw(fr, design$total)
      else .expected_event_1arm(enr, fr, design$total)
    }
    n_events <- ev1(design$enroll_e, design$fail_e, design$n_treat) +
      ev1(design$enroll_c, design$fail_c, design$n_ctrl)
  }
  if (is.null(n_events)) stop("Provide n_events, or a legacy design, or a piecewise design.")
  
  ## Analytic (Schoenfeld): SE(log HR) = 1 / sqrt(d * alloc * (1 - alloc)).
  se_alt <- 1 / sqrt(n_events * alloc * (1 - alloc))
  pos_analytic <- .normal_power(z_effect = abs(loghr) / se_alt,
                                alpha = alpha, alternative = alternative)
  
  ## ---- Simulation: piecewise-exp survival + dropout, admin censoring, log-rank
  pos_sim <- NA_real_
  if (do_sim && !is.null(nsim) && nsim > 0) {
    if (!has_survival) {
      warning("'survival' package not available -- skipping TTE simulation.")
    } else {
      if (!is.null(seed)) set.seed(seed)
      zc <- .crit(alpha, alternative)
      n_t <- design$n_treat; n_c <- design$n_ctrl
      hits <- 0L
      for (i in seq_len(nsim)) {
        entry <- c(if (!is.null(design$enroll_e)) .renroll(n_t, design$enroll_e) else rep(0, n_t),
                   if (!is.null(design$enroll_c)) .renroll(n_c, design$enroll_c) else rep(0, n_c))
        tt <- c(.rpwexp(n_t, design$fail_e$duration, design$fail_e$fail_rate),
                .rpwexp(n_c, design$fail_c$duration, design$fail_c$fail_rate))
        dd <- c(.rpwexp(n_t, design$fail_e$duration, design$fail_e$dropout_rate),
                .rpwexp(n_c, design$fail_c$duration, design$fail_c$dropout_rate))
        grp <- c(rep("treat", n_t), rep("ctrl", n_c))
        ## drop subjects enrolled after the data cutoff (possible if the
        ## enrollment table extends beyond total_duration)
        keep <- entry < design$total
        entry <- entry[keep]; tt <- tt[keep]; dd <- dd[keep]; grp <- grp[keep]
        cal_end <- design$total - entry              # entry -> analysis cutoff
        cens <- pmin(dd, cal_end)                    # dropout or admin censoring
        obs  <- pmin(tt, cens)
        ev   <- as.integer(tt <= cens)
        if (sum(ev) == 0) next
        sd_fit <- survival::survdiff(survival::Surv(obs, ev) ~ grp)
        chisq  <- sd_fit$chisq
        ## signed z: positive when treatment has fewer events than expected (benefit).
        idx <- grep("treat", names(sd_fit$n))
        o_e <- sd_fit$obs[idx] - sd_fit$exp[idx]
        z <- sign(-o_e) * sqrt(chisq)
        reject <- if (alternative == "two.sided") abs(z) > zc else z > zc
        if (isTRUE(reject)) hits <- hits + 1L
      }
      pos_sim <- hits / nsim
    }
  }
  
  structure(list(
    endpoint = "time-to-event", method = "conditional power",
    inputs = list(hr = hr, n_events = round(n_events, 1), alloc = alloc,
                  alpha = alpha, alternative = alternative),
    log_hr = loghr, se_phase3 = se_alt, expected_events = n_events,
    pos_analytic = pos_analytic, pos_simulated = pos_sim, nsim = nsim
  ), class = "pos_result")
}

## ===========================================================================
## Dispatcher + pretty printing
## ===========================================================================

## Single entry point. endpoint in {"continuous","binary","tte"}; extra args are
## forwarded to the matching pos_*() function.
prob_success <- function(endpoint = c("continuous", "binary", "tte"), ...) {
  endpoint <- match.arg(endpoint)
  switch(endpoint,
         continuous = pos_continuous(...),
         binary     = pos_binary(...),
         tte        = pos_tte(...))
}

print.pos_result <- function(x, ...) {
  cat("Probability of Success  (", x$method, ")\n", sep = "")
  cat("Endpoint         :", x$endpoint, "\n")
  cat("Design/inputs    :",
      paste(names(x$inputs), unlist(x$inputs), sep = "=", collapse = ", "), "\n")
  if (!is.null(x$effect_size_d)) cat("Std. effect (d)  :", round(x$effect_size_d, 4), "\n")
  if (!is.null(x$risk_difference)) cat("Risk difference  :", round(x$risk_difference, 4),
                                       if (!is.null(x$inputs$n_strata) && x$inputs$n_strata > 1)
                                         sprintf(" (weighted over %d strata)", x$inputs$n_strata) else "", "\n")
  if (!is.null(x$expected_events)) cat("Events (P3)      :", round(x$expected_events, 1), "\n")
  cat("SE (Phase 3)     :", round(x$se_phase3, 5), "\n")
  cat("-----------------------------------------------\n")
  cat("PoS (analytic)   :", sprintf("%.4f", x$pos_analytic), "\n")
  if (!is.na(x$pos_simulated))
    cat("PoS (simulated)  :", sprintf("%.4f", x$pos_simulated),
        sprintf("  [%d sims]", x$nsim), "\n")
  invisible(x)
}

###############################################################################
## WORKED EXAMPLES
## Run with:  Rscript "pos functions.R"   (or source() then call run_examples())
###############################################################################

run_examples <- function() {
  cat("\n=================================================================\n")
  cat("EXAMPLE 1  Continuous endpoint\n")
  cat("Phase 2 observed mean difference = 3.0, pooled SD = 8.0.\n")
  cat("Planned Phase 3: 150 patients per arm, two-sided 5%.\n")
  cat("=================================================================\n")
  ex1 <- pos_continuous(delta = 3.0, sd = 8.0, n_treat = 150, n_ctrl = 150,
                        alpha = 0.05, alternative = "two.sided",
                        nsim = 10000, seed = 1)
  print(ex1)
  
  ## Same example starting from raw Phase 2 data instead of summary stats:
  set.seed(42)
  y_c <- rnorm(60, 10, 8); y_t <- rnorm(60, 13, 8)
  s <- ph2_summary_continuous(y_t, y_c)
  cat("\n(Recomputed from raw P2 data: delta =", round(s$delta, 2),
      ", pooled SD =", round(s$sd, 2), ")\n")
  
  cat("\n=================================================================\n")
  cat("EXAMPLE 2a  Binary endpoint -- UNSTRATIFIED Miettinen-Nurminen\n")
  cat("Phase 2 response rates: treatment 0.45, control 0.30.\n")
  cat("Planned Phase 3: 220 per arm, two-sided 5%.\n")
  cat("=================================================================\n")
  ex2 <- pos_binary(p_treat = 0.45, p_ctrl = 0.30, n_treat = 220, n_ctrl = 220,
                    alpha = 0.05, alternative = "two.sided",
                    nsim = 10000, seed = 2)
  print(ex2)
  
  cat("\n=================================================================\n")
  cat("EXAMPLE 2b  Binary endpoint -- STRATIFIED M&N (3 strata)\n")
  cat("Stratifier: disease severity (mild/moderate/severe).\n")
  cat("Per-stratum control rates 0.25/0.30/0.40, treatment 0.40/0.45/0.55;\n")
  cat("stratum prevalence 0.40/0.35/0.25; 220 per arm, two-sided 5%.\n")
  cat("=================================================================\n")
  ex2b <- pos_binary(p_treat = c(0.40, 0.45, 0.55),
                     p_ctrl  = c(0.25, 0.30, 0.40),
                     n_treat = 220, n_ctrl = 220,
                     strata_prop = c(0.40, 0.35, 0.25),
                     margin = 0, weights = "MH",
                     alpha = 0.05, alternative = "two.sided",
                     nsim = 10000, seed = 22)
  print(ex2b)
  cat("Per-stratum design:\n"); print(ex2b$per_stratum, row.names = FALSE)
  
  cat("\n(Non-inferiority variant: margin = -0.10, one-sided 2.5%)\n")
  ex2c <- pos_binary(p_treat = 0.30, p_ctrl = 0.30, n_treat = 300,
                     margin = -0.10, alpha = 0.025, alternative = "one.sided",
                     nsim = 5000, seed = 23)
  print(ex2c)
  
  cat("\n=================================================================\n")
  cat("EXAMPLE 3a  Time-to-event -- legacy interface (+ dropout)\n")
  cat("Phase 2 estimated HR = 0.70 (treatment vs control).\n")
  cat("Planned Phase 3: 350/arm, control median 12 mo, 18 mo accrual,\n")
  cat("12 mo additional follow-up, dropout 0.1%/mo, two-sided 5%.\n")
  cat("=================================================================\n")
  ex3 <- pos_tte(hr = 0.70, n_treat = 350, n_ctrl = 350,
                 median_ctrl = 12, enroll_time = 18, follow_time = 12,
                 dropout_rate = 0.001,
                 alpha = 0.05, alternative = "two.sided",
                 nsim = 3000, seed = 3)
  print(ex3)
  
  cat("\n=================================================================\n")
  cat("EXAMPLE 3b  Time-to-event -- PIECEWISE design (gsDesign2-style)\n")
  cat("Enrollment ramp-up: 10/mo for 6 mo, then 25/mo for 22 mo (N = 610).\n")
  cat("Control hazard: median 9 mo for the first 3 mo, median 18 mo after\n")
  cat("(piecewise exponential); dropout 0.1%/mo; HR = 0.70 (PH); cutoff 40 mo.\n")
  cat("=================================================================\n")
  ex3b <- pos_tte(hr = 0.70,
                  enroll_rate = define_enroll_rate(duration = c(6, 22),
                                                   rate = c(10, 25)),
                  fail_rate = define_fail_rate(duration = c(3, 100),
                                               fail_rate = log(2) / c(9, 18),
                                               dropout_rate = 0.001),
                  total_duration = 40, ratio = 1,
                  alpha = 0.05, alternative = "two.sided",
                  nsim = 3000, seed = 33)
  print(ex3b)
  
  cat("\n(Analytic-only TTE from a target event count, no simulation:)\n")
  ex3c <- pos_tte(hr = 0.70, n_events = 250, alloc = 0.5,
                  alpha = 0.05, alternative = "two.sided", nsim = 0)
  print(ex3c)
  
  cat("\n=================================================================\n")
  cat("EXAMPLE 4  Same call via the prob_success() dispatcher\n")
  cat("=================================================================\n")
  print(prob_success("continuous", delta = 3, sd = 8, n_treat = 150,
                     nsim = 0))
  invisible(NULL)
}

## Auto-run examples when executed as a script (Rscript "pos functions.R").
if (!interactive() && sys.nframe() == 0L) {
  run_examples()
}