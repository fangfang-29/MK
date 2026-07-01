###############################################################
## Date: 2026-07-01
## Author: Fangfang Jiang
## ------------------------------------------------------------
##
## Functions:
##    - pwexp_deltas()          :  Compute per-interval durations Delta_k(t)
##    - .validate_tau_cutoffs() :  Normalize tau_cutoffs (single Inf = single exponential)
##    - .parse_formula()        :  Extract components from regression formula, e.g., Surv(time,event) ~ trt + covs
##    - .normalize_imputation_covariates():  Accept NULL / formula / character imp-cov arg
##    - .expand_factor_covariates(): Expand factor / character columns into 0/1 dummies
##    - .extract_arm()          :  Extract arm-level pwexp estimates from .fit_pwexp() output
##    - .draw_imputed_time()    :  Draw imputed time from trt-disc T* for one recipient
##    - rubins_rules()          :  Pool M imputations via Rubin's rules on the log-HR scale
##    - get_true_hr_analytical():  Analytical on-treatment HR (lam_on_trt / lam_on_ctrl)
###############################################################

###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(MASS)
})

# ============================================================================
# (1) Piecewise-exp geometry
# ============================================================================

# For time t and cutoffs tau_cutoffs = c(tau_1, ..., tau_{K-1}), returns
# Delta_k(t) = time spent in interval k up to t, for k = 1, ..., K.
pwexp_deltas <- function(t, tau_cutoffs) {
  K <- length(tau_cutoffs) + 1
  taus_full <- c(0, tau_cutoffs, Inf)
  deltas <- matrix(NA_real_, nrow = length(t), ncol = K)
  for (k in seq_len(K)) {
    deltas[, k] <- pmax(0, pmin(t, taus_full[k + 1]) - taus_full[k])
  }
  deltas
}

# Normalize tau_cutoffs. Single `Inf` means single-exponential (K = 1);
# otherwise returns a strictly increasing finite vector.
.validate_tau_cutoffs <- function(tau_cutoffs) {
  if (length(tau_cutoffs) == 1 && is.infinite(tau_cutoffs)) return(numeric(0))
  if (length(tau_cutoffs) < 1)                              return(numeric(0))
  if (any(!is.finite(tau_cutoffs))) {
    stop("tau_cutoffs must be finite (single Inf = single-exponential).")
  }
  if (is.unsorted(tau_cutoffs, strictly = TRUE)) {
    stop("tau_cutoffs must be strictly increasing.")
  }
  tau_cutoffs
}


# ============================================================================
# (2) Formula & covariate parsing
# ============================================================================

# Parse a two-sided analysis formula: Surv(time, event) ~ treatment + covariates.
# Returns a list(time_var, event_expr, trt_var, covariates, surv_formula).
.parse_formula <- function(formula, treatment = "trt01p", dat) {

  if (is.null(formula)) {
    stop("`formula` must be specified.")
  }
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("`formula` must be a two-sided formula: ",
         "Surv(time, event) ~ treatment + covariates")
  }

  # LHS must be a Surv() call
  lhs <- formula[[2]]
  if (!is.call(lhs) || !identical(lhs[[1]], as.name("Surv"))) {
    stop("LHS of `formula` must be a Surv() call, e.g. Surv(aval, 1 - cnsr).")
  }
  lhs_args <- as.list(lhs)[-1]
  if (length(lhs_args) < 2) {
    stop("Surv() must have at least two arguments: Surv(time, event).")
  }
  time_var   <- deparse(lhs_args[[1]])
  event_expr <- deparse(lhs_args[[2]])
  surv_str   <- deparse(lhs)

  # RHS must contain at least the treatment term
  rhs_vars <- all.vars(formula[[3]])
  if (length(rhs_vars) == 0) {
    stop("RHS of `formula` is empty. It must contain at least the treatment ",
         "variable '", treatment, "'.")
  }
  if (!(treatment %in% rhs_vars)) {
    stop("Treatment variable '", treatment, "' not found in the formula RHS. ",
         "RHS contains: ", paste(rhs_vars, collapse = ", "))
  }

  covariates <- setdiff(rhs_vars, treatment)

  # Every referenced column must exist in dat
  event_vars   <- all.vars(lhs_args[[2]])
  needed       <- unique(c(time_var, event_vars, treatment, covariates))
  missing_cols <- setdiff(needed, names(dat))
  if (length(missing_cols) > 0) {
    stop("Column(s) referenced in `formula` not found in dat: ",
         paste(missing_cols, collapse = ", "))
  }

  list(
    time_var     = time_var,
    event_expr   = event_expr,
    trt_var      = treatment,
    covariates   = covariates,
    surv_formula = surv_str
  )
}

# Normalize the `imputation_covariates` argument to a character vector of
# column names. Accepts:
#   - NULL           -> fall back to `analysis_covariates`
#   - character()    -> no imputation covariates (empty allowed)
#   - character(n)   -> explicit list
#   - a formula      -> variables from RHS (e.g. ~ age + sex)
.normalize_imputation_covariates <- function(imputation_covariates,
                                             analysis_covariates,
                                             treatment) {
  if (is.null(imputation_covariates)) return(analysis_covariates)

  if (inherits(imputation_covariates, "formula")) {
    v <- all.vars(imputation_covariates)
    return(setdiff(v, treatment))   # ignore treatment if user included it
  }
  if (is.character(imputation_covariates)) {
    return(setdiff(imputation_covariates, treatment))
  }
  stop("`imputation_covariates` must be NULL, a character vector, or a ",
       "one-sided formula (e.g. ~ age + sex).")
}

# For any factor / character covariate in `dat`, build 0/1 dummy columns and
# append them to `dat`. Update the analysis and imputation covariate name
# vectors to reference the dummy columns instead of the original factor.
#
# Returns list(dat, analysis_covariates, imputation_covariates, dummy_map).
# dummy_map is a named list mapping original factor -> character vector of
# dummy column names (useful for diagnostics / summaries).
.expand_factor_covariates <- function(dat, analysis_covariates,
                                      imputation_covariates) {
  all_covs <- union(analysis_covariates, imputation_covariates)
  factor_covs <- all_covs[vapply(all_covs, function(c)
    is.factor(dat[[c]]) || is.character(dat[[c]]), logical(1))]

  if (length(factor_covs) == 0) {
    return(list(
      dat                    = dat,
      analysis_covariates    = analysis_covariates,
      imputation_covariates  = imputation_covariates,
      dummy_map              = list()
    ))
  }

  dummy_map <- list()
  for (fc in factor_covs) {
    dat[[fc]] <- factor(dat[[fc]])
    mm <- model.matrix(as.formula(paste("~", fc)), data = dat)
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    dnms <- colnames(mm)
    for (j in seq_along(dnms)) {
      dat[[dnms[j]]] <- as.numeric(mm[, j])
    }
    dummy_map[[fc]] <- dnms
  }

  expand_names <- function(v) {
    out <- character(0)
    for (nm in v) {
      out <- c(out, if (nm %in% names(dummy_map)) dummy_map[[nm]] else nm)
    }
    out
  }

  list(
    dat                    = dat,
    analysis_covariates    = expand_names(analysis_covariates),
    imputation_covariates  = expand_names(imputation_covariates),
    dummy_map              = dummy_map
  )
}


# ============================================================================
# (3) Model-fit extraction
# ============================================================================

# Extract arm-level piecewise-exp fit output from .fit_pwexp() (in rd_mi.R).
# `fit` has one row per (trt, interval) plus list-columns beta_hat / Sigma_hat.
.extract_arm <- function(fit, trt) {
  fit_trt <- fit |> filter(trt01p == trt) |> arrange(interval)
  list(
    alpha_hat  = fit_trt$alpha_hat,
    lambda_hat = fit_trt$lambda_hat,
    beta_hat   = fit_trt$beta_hat[[1]],
    Sigma_hat  = fit_trt$Sigma_hat[[1]],
    d          = fit_trt$d,
    R          = fit_trt$R
  )
}


# ============================================================================
# (4) Step 2 draw: imputed event time for a single subject
# ============================================================================
# Recipients: P5 (offipday = 0) and P6 (offipday > 0, study dropout).
#
# Draw u ~ Uniform(0, S(offipday)), invert the piecewise-exp cumulative
# hazard to find T*, then re-censor at EOS.
#
# Note on re-censoring: paper Sec 4.3 formula 5 followup writes
# min(EOS, non-CV death); this implementation drops non-CV death from the
# cutoff because the vital-records tracking that assumption relies on is
# typically unavailable for real study dropouts.
#
# Arguments:
#    lam         : piecewise-exp hazard rates (length K), possibly cov-adjusted
#    offipday    : truncation time in off-treatment clock (>= 0)
#    tau_cutoffs : interval cutoffs (finite vector; empty = single exp)
#    eos         : end-of-study in days (rand clock)
#    aval_obs    : observed time on study at dropout (rand clock)
#
# Returns list(t_star, aval_imp, cnsr_imp).
.draw_imputed_time <- function(lam, offipday, tau_cutoffs, eos, aval_obs) {

  tau_cutoffs <- .validate_tau_cutoffs(tau_cutoffs)
  K <- length(lam)
  taus_full <- c(0, tau_cutoffs, Inf)

  # Survival at truncation point offipday
  D_trunc <- pwexp_deltas(offipday, tau_cutoffs)[1, ]
  S_trunc <- exp(-sum(lam * D_trunc))

  if (S_trunc < 1e-15 || !is.finite(S_trunc)) {
    return(list(t_star = Inf, aval_imp = aval_obs, cnsr_imp = 1L))
  }

  # Truncation trick: u ~ Uniform(0, S_trunc)
  u <- runif(1, min = 0, max = S_trunc)

  # Invert S(T*) = u by walking cumulative hazard interval by interval
  cumhaz <- 0
  t_star <- Inf
  for (k in seq_len(K)) {
    tau_lower <- taus_full[k]
    tau_upper <- taus_full[k + 1]

    if (is.infinite(tau_upper)) {
      if (lam[k] > 0) t_star <- tau_lower - (log(u) + cumhaz) / lam[k]
      break
    }
    width_k    <- tau_upper - tau_lower
    cumhaz_new <- cumhaz + lam[k] * width_k
    S_upper    <- exp(-cumhaz_new)

    if (u > S_upper) {
      if (lam[k] > 0) t_star <- tau_lower - (log(u) + cumhaz) / lam[k]
      break
    }
    cumhaz <- cumhaz_new
  }

  # Sanity: T* must be >= offipday (left-truncated)
  if (!is.infinite(t_star) && t_star < offipday - 1e-9) t_star <- offipday

  aval_imp <- aval_obs + (t_star - offipday)

  # Re-censor at EOS only
  if (aval_imp > eos || !is.finite(aval_imp)) {
    aval_imp <- eos
    cnsr_imp <- 1L
  } else {
    cnsr_imp <- 0L
  }

  list(t_star = t_star, aval_imp = aval_imp, cnsr_imp = cnsr_imp)
}


# ============================================================================
# (5) Step 4 pooling: Rubin's rules on the log-HR scale
# ============================================================================
# Inputs: list of length M, each element list(logHR, var_logHR).
# Returns Rubin summary (logHR, HR, SE, CI, df, p_value, W, B, T, M).
rubins_rules <- function(cox_results, alpha = 0.05) {
  M <- length(cox_results)
  if (M < 1) stop("cox_results must contain at least one imputation.")

  logHRs <- sapply(cox_results, `[[`, "logHR")
  varHRs <- sapply(cox_results, `[[`, "var_logHR")

  Q_bar   <- mean(logHRs)             # pooled logHR
  W       <- mean(varHRs)             # within-imp variance
  B       <- if (M > 1) var(logHRs) else 0    # between-imp variance
  T_total <- W + (1 + 1 / M) * B      # total variance
  se_total <- sqrt(T_total)

  # Rubin's large-sample degrees of freedom
  nu <- if (B > 0 && M > 1) (M - 1) * (1 + W / ((1 + 1 / M) * B))^2 else Inf

  t_stat  <- Q_bar / se_total
  p_value <- 2 * pt(-abs(t_stat), df = nu)
  crit    <- qt(1 - alpha / 2, df = nu)
  ci_low  <- Q_bar - crit * se_total
  ci_high <- Q_bar + crit * se_total

  HR_pooled <- exp(Q_bar)
  se_pooled <- HR_pooled * se_total

  list(
    logHR = Q_bar, HR = HR_pooled,
    SE_logHR = se_total, SE_HR = se_pooled,
    CI_logHR = c(ci_low, ci_high), CI_HR = exp(c(ci_low, ci_high)),
    alpha = alpha, alternative = "two.sided",
    test_statistic = t_stat, df = nu, p_value = p_value,
    within_variance = W, between_variance = B, total_variance = T_total,
    M = M
  )
}


# ============================================================================
# (6) Analytical true HR
# ============================================================================
get_true_hr_analytical <- function(lam_on_trt, lam_on_ctrl) {

  if (!is.numeric(lam_on_trt) || length(lam_on_trt) != 1 ||
      !is.numeric(lam_on_ctrl) || length(lam_on_ctrl) != 1) {
    stop("`lam_on_trt` and `lam_on_ctrl` must each be a single numeric value.")
  }
  if (lam_on_ctrl <= 0) {
    stop("`lam_on_ctrl` must be > 0.")
  }
  if (lam_on_trt < 0) {
    stop("`lam_on_trt` must be >= 0.")
  }

  lam_on_trt / lam_on_ctrl
}
