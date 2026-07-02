###############################################################
## Date: 2026-07-01
## Author: Fangfang Jiang
## ------------------------------------------------------------
## Retrieved-Dropout Multiple Imputation for Time-to-Event Data
##
## The main function `rd_mi()` supports:
##   - with or without covariates
##   - three imputation methods (improper / proper_like / proper)
##   - user-configurable schema column names
##
## Depends on helper_functions.R.
##
## Method structure :
##   Step 1: draw M sets of (lambda, beta) per treatment arm
##   Step 2: for each recipient (P5, P6), draw imputed event time T*
##   Step 3: fit Cox PH on each completed dataset
##   Step 4: pool via Rubin's rules on the log-HR scale
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(MASS)
})

# source("helper_functions.R")


# ============================================================================
# INTERNAL Step 0: Fit piecewise-exp model with (optional) covariates
# ============================================================================
# Poisson GLM on person-period expanded data, which is equivalent to the
# piecewise-exponential likelihood.
# When `covariates` is empty, the Poisson GLM reduces to lambda_hat_k = d_k / R_k.
#
# Returns a data.frame with one row per (trt, interval), containing:
#   alpha_hat, lambda_hat, beta_hat (list), Sigma_hat (list), d, R, n_donors,
#   coef_names (list).
.fit_pwexp <- function(dat, tau_cutoffs, covariates = character(0),
                       trt_var = "trt01p") {
  tau_cutoffs <- .validate_tau_cutoffs(tau_cutoffs)

  # Donor set = {P2, P4, P6, P7} with observed off-treatment days > 0
  donors <- dat |> filter(RD == 1, !is.na(offipday), offipday > 0)

  K         <- length(tau_cutoffs) + 1L
  tau_lower <- c(0, tau_cutoffs)
  tau_upper <- c(tau_cutoffs, Inf)
  has_cov   <- length(covariates) > 0

  out <- donors |>
    group_by(.data[[trt_var]]) |>
    group_modify(function(.x, .y) {

      t_off  <- .x$offipday
      # Event during off-treatment = "off_cnsr == 0"
      ev_off <- as.integer(!is.na(.x$off_cnsr) & .x$off_cnsr == 0)
      n_don  <- nrow(.x)

      if (n_don < 5) {
        warning("Donor set too small (n = ", n_don, ") for ", trt_var, " = ",
                .y[[trt_var]], ": imputation model may be unstable.")
      }

      # Time in each interval + interval containing event (NA if censored)
      D <- pwexp_deltas(t_off, tau_cutoffs)
      ev_interval <- rep(NA_integer_, n_don)
      ev_interval[ev_off == 1] <- findInterval(
        t_off[ev_off == 1], tau_cutoffs) + 1L

      # Expand to person-period rows
      rows <- vector("list", n_don * K)
      idx  <- 0L
      for (i in seq_len(n_don)) {
        for (k in seq_len(K)) {
          delta_ik <- D[i, k]
          if (delta_ik <= 0) next
          idx <- idx + 1L
          row <- data.frame(
            interval  = k,
            d_ik      = as.integer(!is.na(ev_interval[i]) && ev_interval[i] == k),
            log_delta = log(delta_ik)
          )
          if (has_cov) row <- cbind(row, .x[i, covariates, drop = FALSE])
          rows[[idx]] <- row
        }
      }
      pp <- do.call(rbind, rows[seq_len(idx)])
      pp$interval_f <- factor(pp$interval, levels = seq_len(K))

      # Build Poisson GLM formula. K = 1 needs a plain intercept because a
      # 1-level factor can't take contrasts.
      cov_str <- if (has_cov) {
        paste("+", paste(sprintf("`%s`", covariates), collapse = " + "))
      } else ""

      glm_formula <- if (K == 1L) {
        as.formula(paste("d_ik ~ 1", cov_str, "+ offset(log_delta)"))
      } else {
        as.formula(paste("d_ik ~ -1 + interval_f", cov_str, "+ offset(log_delta)"))
      }

      glm_fit <- glm(glm_formula, data = pp, family = poisson(link = "log"))
      coefs   <- coef(glm_fit)
      Sigma   <- vcov(glm_fit)

      alpha_idx <- if (K == 1L) which(names(coefs) == "(Intercept)") else
        grep("^interval_f", names(coefs))
      beta_idx  <- setdiff(seq_along(coefs), alpha_idx)

      alpha_hat <- unname(coefs[alpha_idx])
      beta_hat  <- if (has_cov) unname(coefs[beta_idx]) else numeric(0)

      # Raw d and R for the Gamma-conjugate posterior in the no-cov proper method
      d_vec <- tabulate(ev_interval[!is.na(ev_interval)], nbins = K)
      R_vec <- colSums(D)

      data.frame(
        interval   = seq_len(K),
        tau_lower  = tau_lower,
        tau_upper  = tau_upper,
        alpha_hat  = alpha_hat,
        lambda_hat = exp(alpha_hat),
        d          = as.numeric(d_vec),
        R          = as.numeric(R_vec),
        n_donors   = n_don,
        beta_hat   = I(rep(list(beta_hat), K)),
        Sigma_hat  = I(rep(list(Sigma),    K)),
        coef_names = I(rep(list(names(coefs)), K))
      )
    }) |>
    ungroup()

  # Standardize the grouping column name to `trt01p` so downstream helpers
  # (.extract_arm) can filter on it uniformly regardless of the user's
  # actual treatment column.
  if (trt_var != "trt01p") {
    names(out)[names(out) == trt_var] <- "trt01p"
  }
  out
}


# ============================================================================
# INTERNAL Step 1 samplers: draw (lambda, beta) for one arm
# ============================================================================
# All three samplers return a list of length M, each element list(lambda, beta).

# Method 1 (improper): fixed MLE, shared across M
.sample_improper <- function(fit, trt, M, ...) {
  arm    <- .extract_arm(fit, trt)
  params <- list(lambda = arm$lambda_hat, beta = arm$beta_hat)
  replicate(M, params, simplify = FALSE)
}

# Method 2 (proper_like): draw (alpha, beta) jointly from MVN(MLE, Sigma_hat)
.sample_proper_like <- function(fit, trt, M, ...) {
  arm   <- .extract_arm(fit, trt)
  K     <- length(arm$alpha_hat)
  p     <- length(arm$beta_hat)
  mu    <- c(arm$alpha_hat, arm$beta_hat)
  Sigma <- arm$Sigma_hat

  lapply(seq_len(M), function(m) {
    draw <- MASS::mvrnorm(1, mu = mu, Sigma = Sigma)
    list(
      lambda = exp(draw[seq_len(K)]),
      beta   = if (p > 0) draw[K + seq_len(p)] else numeric(0)
    )
  })
}

# Method 3 (proper): MCMC posterior draws.
#   - No cov: direct Gamma conjugate draws (closed form).
#   - With cov: ARMS-within-Gibbs (armspp), matching SAS PROC PHREG bayes
#     which also uses ARMS internally.
#
# `dat` needed for the covariate case to reconstruct the design matrix.
.sample_proper <- function(fit, trt, M, dat = NULL, covariates = character(0),
                           tau_cutoffs = numeric(0),
                           mcmc_control = list()) {

  # Defaults for MCMC control
  ctrl <- utils::modifyList(list(
    nburnin  = 1000, thin = 100, n_sd = 50,
    a0 = 1e-4, b0 = 1e-2, var_beta = 1e3
  ), mcmc_control)

  arm     <- .extract_arm(fit, trt)
  K       <- length(arm$alpha_hat)
  p       <- length(arm$beta_hat)
  has_cov <- p > 0

  # --- No-cov: exact Gamma conjugate draws ---
  if (!has_cov) {
    draws <- lapply(seq_len(M), function(m) {
      lambda_draw <- vapply(seq_len(K), function(k)
        rgamma(1, shape = arm$d[k] + ctrl$a0, rate = arm$R[k] + ctrl$b0),
        numeric(1))
      list(lambda = lambda_draw, beta = numeric(0))
    })
    attr(draws, "sampler") <- "Gamma-conjugate"
    return(draws)
  }

  # --- Cov case: ARMS-within-Gibbs (requires armspp) ---
  if (!requireNamespace("armspp", quietly = TRUE)) {
    stop("Package 'armspp' is required for `proper` with covariates. ",
         "Install with install.packages('armspp').")
  }
  if (is.null(dat)) {
    stop("`dat` must be supplied to .sample_proper() for the covariate case.")
  }

  donors <- dat |> filter(trt01p == trt, RD == 1, !is.na(offipday), offipday > 0)
  D      <- pwexp_deltas(donors$offipday, .validate_tau_cutoffs(tau_cutoffs))
  ev_off <- as.integer(!is.na(donors$off_cnsr) & donors$off_cnsr == 0)
  ev_int <- rep(NA_integer_, nrow(donors))
  ev_int[ev_off == 1] <- findInterval(
    donors$offipday[ev_off == 1], .validate_tau_cutoffs(tau_cutoffs)) + 1L

  X <- model.matrix(
    as.formula(paste("~", paste(sprintf("`%s`", covariates), collapse = " + "))),
    data = donors
  )
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) != p) {
    stop("Design matrix has ", ncol(X), " column(s) but fit estimated ", p,
         " beta coefficient(s). `covariates` may not match .fit_pwexp() input.")
  }
  d_k <- tabulate(ev_int[!is.na(ev_int)], nbins = K)

  # Log-posterior on theta = c(alpha, beta)
  log_post <- function(theta) {
    alpha  <- theta[seq_len(K)]
    beta   <- theta[K + seq_len(p)]
    lam    <- exp(alpha)
    eta    <- as.numeric(X %*% beta)
    wt     <- exp(eta)
    R_beta <- colSums(D * wt)

    pwexp_ll   <- sum(d_k * alpha - lam * R_beta)   # log-likelihood, alpha part
    cov_ll     <- sum(eta[ev_off == 1])              # log-likelihood, beta part
    lam_prior  <- sum((ctrl$a0 - 1) * alpha - ctrl$b0 * lam)  # Gamma(a0,b0) on lambda
    beta_prior <- -0.5 * sum(beta^2) / ctrl$var_beta          # MVN(0, var_beta I)

    pwexp_ll + cov_ll + lam_prior + beta_prior
  }

  # Bounded support for ARMS: MLE +/- n_sd asymptotic SE per coordinate
  se_theta  <- sqrt(diag(arm$Sigma_hat))
  theta_hat <- c(arm$alpha_hat, arm$beta_hat)
  lower_vec <- theta_hat - ctrl$n_sd * se_theta
  upper_vec <- theta_hat + ctrl$n_sd * se_theta

  dim_theta <- K + p
  theta_cur <- theta_hat
  n_iter    <- ctrl$nburnin + M * ctrl$thin
  draws     <- vector("list", M)
  draw_idx  <- 0L
  chain     <- matrix(NA_real_, nrow = n_iter - ctrl$nburnin, ncol = dim_theta)

  for (iter in seq_len(n_iter)) {
    # Systematic-scan Gibbs: update each theta_j via ARMS from its univariate full conditional.
    for (j in seq_len(dim_theta)) {
      log_fc <- function(x) {
        th    <- theta_cur
        th[j] <- x
        log_post(th)
      }
      theta_cur[j] <- armspp::arms(
        n_samples = 1, log_pdf = log_fc,
        lower = lower_vec[j], upper = upper_vec[j],
        previous = theta_cur[j]
      )
    }

    if (iter > ctrl$nburnin) chain[iter - ctrl$nburnin, ] <- theta_cur

    if (iter > ctrl$nburnin && ((iter - ctrl$nburnin) %% ctrl$thin == 0)) {
      draw_idx <- draw_idx + 1L
      draws[[draw_idx]] <- list(
        lambda = exp(theta_cur[seq_len(K)]),
        beta   = theta_cur[K + seq_len(p)]
      )
    }
  }

  param_names <- c(paste0("alpha", seq_len(K)), colnames(X))
  colnames(chain) <- param_names
  attr(draws, "chain")       <- chain
  attr(draws, "nburnin")     <- ctrl$nburnin
  attr(draws, "thin")        <- ctrl$thin
  attr(draws, "param_names") <- param_names
  attr(draws, "sampler")     <- "ARMS-within-Gibbs"

  draws
}


# ============================================================================
# INTERNAL Step 2: impute all recipients (P5, P6) for each of M draws
# ============================================================================
.impute_recipients <- function(dat, draws_active, draws_placebo,
                               tau_cutoffs, trt_var, covariates) {

  tau_cutoffs <- .validate_tau_cutoffs(tau_cutoffs)
  has_cov     <- length(covariates) > 0
  M           <- length(draws_active)

  completed <- vector("list", M)

  for (m in seq_len(M)) {
    dat_imp  <- dat
    params_a <- draws_active[[m]]
    params_p <- draws_placebo[[m]]

    for (i in seq_len(nrow(dat))) {
      if (!(dat$pattern[i] %in% c("P5", "P6"))) next

      params   <- if (dat[[trt_var]][i] == 1) params_a else params_p
      lam_base <- params$lambda
      beta     <- params$beta

      # Subject-specific hazard: lam_i = lam_base * exp(beta' x_i)
      if (has_cov && length(beta) > 0) {
        x_i   <- as.numeric(dat[i, covariates, drop = FALSE])
        lam_i <- lam_base * exp(sum(beta * x_i))
      } else {
        lam_i <- lam_base
      }

      offipday <- if (is.na(dat$offipday[i])) 0 else dat$offipday[i]

      res <- .draw_imputed_time(
        lam         = lam_i,
        offipday    = offipday,
        tau_cutoffs = tau_cutoffs,
        eos         = dat$eos[i],
        aval_obs    = dat$aval[i]
      )
      dat_imp$aval[i] <- res$aval_imp
      dat_imp$cnsr[i] <- res$cnsr_imp
    }
    completed[[m]] <- dat_imp
  }
  completed
}


# ============================================================================
# INTERNAL Step 3: fit Cox PH per completed dataset
# ============================================================================
.fit_cox_per_imp <- function(completed_datasets, analysis_formula, trt_var) {
  lapply(completed_datasets, function(dat) {
    fit <- coxph(analysis_formula, data = dat)
    list(
      logHR     = unname(coef(fit)[trt_var]),
      var_logHR = vcov(fit)[trt_var, trt_var]
    )
  })
}


# ============================================================================
# rd_mi() -- MAIN USER-FACING FUNCTION
# ============================================================================
# Arguments
# ---------
# data
#     data.frame with the endpoint, treatment, covariates, and RD-MI schema
#     columns (see *_col arguments below).
#
# formula
#     Analysis Cox formula: Surv(time, event) ~ treatment + covariates.
#     Default: Surv(aval, 1 - cnsr) ~ trt01p (no covariates).
#
# treatment
#     Name of the treatment column, must appear in the formula RHS.
#     Coded 0 = control / 1 = active treatment.
#
# imputation_covariates
#     Covariates for the piecewise-exponential imputation model. Accepts:
#       - character(0) (default): no covariates in imputation
#       - character(n)          : explicit column names, e.g. c("age", "sex")
#       - a formula             : e.g. ~ age + sex
#       - NULL                  : opt-in to reusing the analysis covariates
#                                 (short-hand for "same as analysis model")
#
# pattern_col, rd_col, offipday_col, off_event_col, eos_col
#     Column names in `data` holding the RD-MI schema:
#       pattern    P1..P8 categorization
#       RD         retrieved-dropout indicator (1 for P2/P4/P6/P7)
#       offipday   observed off-treatment days
#       off_event  off-treatment event indicator (0 = event, 1 = censored,
#                                                 NA / -1 = not applicable)
#       eos        end-of-study day per subject
#
# method
#     One or more of "improper" / "proper_like" / "proper".
#
# M
#     Number of imputations.
#
# tau_cutoffs
#     Piecewise-exp interval boundaries (finite vector; `Inf` = single exp).
#
# mcmc_control
#     Named list overriding defaults for the `proper` MCMC:
#       nburnin = 1000, thin = 100, n_sd = 50,
#       a0 = 1e-4, b0 = 1e-2, var_beta = 1e3
#
# seed
#     RNG seed for reproducibility (applied once at the top).
#
# keep_imputations
#     If TRUE, return the completed datasets in `$imputations`. Off by
#     default because they can be large.
#
# Returns
# -------
# An `rd_mi` object (list with class "rd_mi"):
#   $results     data.frame summary, one row per method
#   $details     list keyed by method: full Rubin's output + lambda/beta draws
#   $imputations list keyed by method: completed datasets (if keep_imputations)
#   $formula     analysis formula (with dummy-expanded factor covariates)
#   $imputation_covariates final imp-covariate names (after dummy expansion)
#   $M, $tau_cutoffs, $method, $call metadata
rd_mi <- function(
    data,
    formula               = Surv(aval, 1 - cnsr) ~ trt01p,
    treatment             = "trt01p",
    imputation_covariates = character(0),   # default: no cov in imputation
    pattern_col           = "pattern",
    rd_col                = "RD",
    offipday_col          = "offipday",
    off_event_col         = "off_cnsr",
    eos_col               = "eos",
    method                = c("improper", "proper_like", "proper"),
    M                     = 25,
    tau_cutoffs           = Inf,
    mcmc_control          = list(
      nburnin  = 1000,   # burn-in iterations (discarded)
      thin     = 100,    # keep 1 draw every `thin` sweeps
      n_sd     = 50,     # ARMS bounded-support width in asymptotic SEs
      a0       = 1e-4,   # Gamma prior shape on lambda_k
      b0       = 1e-2,   # Gamma prior rate on lambda_k
      var_beta = 1e3     # MVN prior variance on beta
    ),
    seed                  = NULL,
    keep_imputations      = FALSE
) {

  call_ <- match.call()
  method <- match.arg(method,
                      choices     = c("improper", "proper_like", "proper"),
                      several.ok  = TRUE)
  if (!is.null(seed)) set.seed(seed)
  tau_cutoffs <- .validate_tau_cutoffs(tau_cutoffs)

  # Check if the time-to-event variables exist
  paper_cols <- c("aval", "cnsr")
  miss_paper <- setdiff(paper_cols, names(data))
  if (length(miss_paper) > 0) {
    stop("`data` must include columns: 'aval' (observed ",
         "time) and 'cnsr' (0 = event, 1 = censored). Missing: ",
         paste(miss_paper, collapse = ", "),
         ". Rename before calling rd_mi().")
  }

  # ------------------------------------------------------------------------
  # 1. Rename user schema columns to internal fixed names
  # ------------------------------------------------------------------------
  schema_map <- c(
    pattern   = pattern_col,
    RD        = rd_col,
    offipday  = offipday_col,
    off_cnsr  = off_event_col,
    eos       = eos_col
  )
  missing_schema <- setdiff(schema_map, names(data))
  if (length(missing_schema) > 0) {
    stop("`data` is missing schema column(s): ",
         paste(missing_schema, collapse = ", "))
  }
  dat <- data
  for (internal_name in names(schema_map)) {
    src <- schema_map[[internal_name]]
    if (src != internal_name && !(internal_name %in% names(dat))) {
      dat[[internal_name]] <- dat[[src]]
    }
  }

  # ------------------------------------------------------------------------
  # 2. Parse analysis formula; extract analysis covariates
  # ------------------------------------------------------------------------
  parsed <- .parse_formula(formula, treatment, dat)
  analysis_covariates <- parsed$covariates

  # ------------------------------------------------------------------------
  # 3. Resolve imputation covariates
  # ------------------------------------------------------------------------
  imp_covariates <- .normalize_imputation_covariates(
    imputation_covariates, analysis_covariates, treatment
  )
  missing_imp <- setdiff(imp_covariates, names(dat))
  if (length(missing_imp) > 0) {
    stop("Imputation covariate column(s) not found in data: ",
         paste(missing_imp, collapse = ", "))
  }

  # ------------------------------------------------------------------------
  # 4. Expand factor / character covariates into dummy columns
  # ------------------------------------------------------------------------
  exp_result          <- .expand_factor_covariates(dat, analysis_covariates, imp_covariates)
  dat                 <- exp_result$dat
  analysis_covariates <- exp_result$analysis_covariates
  imp_covariates      <- exp_result$imputation_covariates
  dummy_map           <- exp_result$dummy_map

  # Rebuild analysis formula with dummy-expanded covariate names
  analysis_cov_str <- if (length(analysis_covariates) > 0) {
    paste("+", paste(sprintf("`%s`", analysis_covariates), collapse = " + "))
  } else ""
  analysis_fml <- as.formula(
    paste(parsed$surv_formula, "~", parsed$trt_var, analysis_cov_str)
  )

  # ------------------------------------------------------------------------
  # 5. Rename treatment column to internal `trt01p` if needed
  # ------------------------------------------------------------------------
  # Several internal steps (.sample_proper, .fit_pwexp, .extract_arm) assume
  # the treatment column is named `trt01p`. Rename here so downstream code
  # doesn't have to thread the user's actual name through every helper.
  if (parsed$trt_var != "trt01p") {
    if ("trt01p" %in% names(dat)) {
      stop("`data` already has a column named 'trt01p' distinct from the ",
           "treatment column '", parsed$trt_var, "'. Rename one of them ",
           "before calling rd_mi().")
    }
    dat$trt01p       <- dat[[parsed$trt_var]]
    parsed$trt_var   <- "trt01p"
  }

  # ------------------------------------------------------------------------
  # 6. Validate treatment coding
  # ------------------------------------------------------------------------
  trt_vals <- unique(stats::na.omit(dat$trt01p))
  if (!all(c(0, 1) %in% trt_vals)) {
    stop("Treatment column must contain both 0 (control) and 1 (treatment).")
  }

  # Rebuild analysis formula against the (possibly renamed) treatment column
  analysis_fml <- as.formula(
    paste(parsed$surv_formula, "~", parsed$trt_var, analysis_cov_str)
  )

  # ------------------------------------------------------------------------
  # 7. Fit piecewise-exp model once (per arm, via Poisson GLM)
  # ------------------------------------------------------------------------
  fit <- .fit_pwexp(dat, tau_cutoffs,
                    covariates = imp_covariates,
                    trt_var    = parsed$trt_var)

  # ------------------------------------------------------------------------
  # 8. Loop over requested methods: Step 1..4
  # ------------------------------------------------------------------------
  sampler_map <- list(
    improper    = .sample_improper,
    proper_like = .sample_proper_like,
    proper      = .sample_proper
  )

  results_rows <- list()
  details      <- list()
  imputations  <- list()

  for (mth in method) {
    sampler <- sampler_map[[mth]]

    # Step 1
    draws_active <- sampler(fit, trt = 1, M = M,
                            dat = dat, covariates = imp_covariates,
                            tau_cutoffs = tau_cutoffs,
                            mcmc_control = mcmc_control)
    draws_placebo <- sampler(fit, trt = 0, M = M,
                             dat = dat, covariates = imp_covariates,
                             tau_cutoffs = tau_cutoffs,
                             mcmc_control = mcmc_control)

    # Step 2
    completed <- .impute_recipients(dat, draws_active, draws_placebo,
                                    tau_cutoffs   = tau_cutoffs,
                                    trt_var       = parsed$trt_var,
                                    covariates    = imp_covariates)

    # Step 3
    cox_res <- .fit_cox_per_imp(completed, analysis_fml, parsed$trt_var)

    # Step 4
    pooled <- rubins_rules(cox_res)

    # Assemble diagnostics
    lam_active  <- do.call(rbind, lapply(draws_active,  `[[`, "lambda"))
    lam_placebo <- do.call(rbind, lapply(draws_placebo, `[[`, "lambda"))
    colnames(lam_active)  <- paste0("lam", seq_len(ncol(lam_active)))
    colnames(lam_placebo) <- paste0("lam", seq_len(ncol(lam_placebo)))

    if (length(imp_covariates) > 0 && length(draws_active[[1]]$beta) > 0) {
      beta_active  <- do.call(rbind, lapply(draws_active,  `[[`, "beta"))
      beta_placebo <- do.call(rbind, lapply(draws_placebo, `[[`, "beta"))
      colnames(beta_active)  <- imp_covariates
      colnames(beta_placebo) <- imp_covariates
    } else {
      beta_active  <- matrix(nrow = M, ncol = 0)
      beta_placebo <- matrix(nrow = M, ncol = 0)
    }

    details[[mth]] <- c(pooled, list(
      lambda_draws_active  = lam_active,
      lambda_draws_placebo = lam_placebo,
      beta_draws_active    = beta_active,
      beta_draws_placebo   = beta_placebo,
      # Full post-burn-in MCMC chain (present only for ARMS-in-Gibbs;
      # NULL for Gamma-conjugate / MVN / fixed samplers). Consumed by
      # diagnose_rd_mi().
      chain_active         = attr(draws_active,  "chain"),
      chain_placebo        = attr(draws_placebo, "chain"),
      mcmc_nburnin         = attr(draws_active,  "nburnin"),
      mcmc_thin            = attr(draws_active,  "thin"),
      sampler_active       = attr(draws_active,  "sampler"),
      sampler_placebo      = attr(draws_placebo, "sampler")
    ))

    if (keep_imputations) imputations[[mth]] <- completed

    results_rows[[mth]] <- data.frame(
      method    = mth,
      logHR     = pooled$logHR,
      HR        = pooled$HR,
      SE_logHR  = pooled$SE_logHR,
      df        = pooled$df,
      p_value   = pooled$p_value,
      CI_HR_low = pooled$CI_HR[1],
      CI_HR_hi  = pooled$CI_HR[2],
      W         = pooled$within_variance,
      B         = pooled$between_variance,
      total_var = pooled$total_variance,
      M         = pooled$M,
      row.names = NULL,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, results_rows)
  rownames(results) <- NULL

  structure(
    list(
      results               = results,
      details               = details,
      imputations           = if (keep_imputations) imputations else NULL,
      formula               = analysis_fml,
      imputation_covariates = imp_covariates,
      dummy_map             = dummy_map,
      M                     = M,
      tau_cutoffs           = tau_cutoffs,
      method                = method,
      call                  = call_
    ),
    class = "rd_mi"
  )
}


# ============================================================================
# diagnose_rd_mi() -- MCMC diagnostics for the `proper` method
# ============================================================================
# For the `proper` method with covariates (ARMS-within-Gibbs), reports ESS,
# Geweke z-statistic, posterior summary, and optionally trace/density plots.
# For the no-covariate proper method (Gamma-conjugate closed-form) the draws
# are iid, so MCMC diagnostics don't apply and only a posterior summary is
# returned.
#
# Arguments:
#   res    : an `rd_mi` object returned by rd_mi()
#   method : method name (default "proper" -- typically the only one to check)
#   arm    : which arm(s) to diagnose: "active", "placebo", or both
#   plot   : if TRUE, opens base-R plot() windows with trace/density panels
#
# Returns: (invisibly) a list keyed by arm, each element containing
#   $sampler, $ess, $geweke, $summary (or $note for the Gamma-conjugate case).
diagnose_rd_mi <- function(res,
                           method = "proper",
                           arm    = c("active", "placebo"),
                           plot   = FALSE) {

  if (!inherits(res, "rd_mi")) {
    stop("`res` must be an rd_mi object (from rd_mi()).")
  }
  if (!method %in% names(res$details)) {
    stop("Method '", method, "' was not run in this rd_mi() call. ",
         "Available: ", paste(names(res$details), collapse = ", "))
  }
  arm <- match.arg(arm, choices = c("active", "placebo"), several.ok = TRUE)

  d <- res$details[[method]]
  out <- list()

  for (a in arm) {
    chain_key   <- paste0("chain_",   a)
    sampler_key <- paste0("sampler_", a)
    lam_key     <- paste0("lambda_draws_", a)
    beta_key    <- paste0("beta_draws_",   a)

    chain   <- d[[chain_key]]
    sampler <- d[[sampler_key]]

    # ---- Gamma-conjugate case: iid closed-form draws, no MCMC diagnostics ----
    if (is.null(chain)) {
      lam_mat  <- d[[lam_key]]
      beta_mat <- d[[beta_key]]
      draws_mat <- cbind(lam_mat, beta_mat)
      cat(sprintf("\n[%s / %s] Sampler: %s -- iid closed-form draws\n",
                  method, a, sampler))
      cat("MCMC diagnostics are not applicable. Posterior summary:\n")
      summ <- data.frame(
        mean = colMeans(draws_mat),
        sd   = apply(draws_mat, 2, sd),
        q025 = apply(draws_mat, 2, quantile, 0.025),
        q500 = apply(draws_mat, 2, median),
        q975 = apply(draws_mat, 2, quantile, 0.975)
      )
      print(round(summ, 4))
      out[[a]] <- list(sampler = sampler,
                       note    = "iid closed-form draws; MCMC diagnostics N/A",
                       summary = summ)
      next
    }

    # ---- ARMS-within-Gibbs case: real MCMC diagnostics via coda ----
    if (!requireNamespace("coda", quietly = TRUE)) {
      stop("Package 'coda' is required for MCMC diagnostics. ",
           "Install with install.packages('coda').")
    }

    thin <- d$mcmc_thin %||% 1
    if (is.null(thin)) thin <- 1

    mcmc_obj <- coda::mcmc(chain, start = 1, thin = thin)
    ess      <- coda::effectiveSize(mcmc_obj)
    geweke   <- coda::geweke.diag(mcmc_obj)
    summ     <- summary(mcmc_obj)

    cat(sprintf("\n[%s / %s] Sampler: %s\n", method, a, sampler))
    cat(sprintf("Chain length (post burn-in): %d, thinning: %d\n",
                nrow(chain), thin))
    cat("\nEffective sample size (ESS):\n"); print(round(ess, 1))
    cat("\nGeweke z-statistic (should be roughly |z| < 2 for stationarity):\n")
    print(round(geweke$z, 2))
    cat("\nPosterior summary:\n"); print(summ)

    if (plot) {
      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mfrow = c(2, 2))
      for (nm in colnames(chain)) {
        plot(coda::mcmc(chain[, nm, drop = FALSE]),
             main = paste(a, "-", nm),
             auto.layout = FALSE)
      }
    }

    out[[a]] <- list(
      sampler  = sampler,
      chain    = mcmc_obj,
      ess      = ess,
      geweke   = geweke,
      summary  = summ
    )
  }

  invisible(out)
}

# Small null-coalesce operator used above
`%||%` <- function(a, b) if (is.null(a)) b else a


# ============================================================================
# print.rd_mi() -- pretty display of the results table
# ============================================================================
print.rd_mi <- function(x, digits = 4, ...) {
  cat("\nRetrieved-Dropout Multiple Imputation\n")
  cat(strrep("-", 66), "\n", sep = "")

  cat("Analysis formula      : ", deparse(x$formula), "\n", sep = "")
  if (length(x$imputation_covariates) > 0) {
    cat("Imputation covariates : ",
        paste(x$imputation_covariates, collapse = ", "), "\n", sep = "")
  } else {
    cat("Imputation covariates : (none)\n")
  }
  cat("M (imputations)       : ", x$M, "\n", sep = "")
  cat("tau_cutoffs           : ",
      if (length(x$tau_cutoffs) == 0) "single exponential"
      else paste(x$tau_cutoffs, collapse = ", "), "\n", sep = "")
  cat(strrep("-", 66), "\n", sep = "")

  r <- x$results
  ci_lvl <- 100 * (1 - x$details[[1]]$alpha)
  print_df <- data.frame(
    Method    = r$method,
    HR        = round(r$HR, digits),
    SE_logHR  = round(r$SE_logHR, digits),
    CI_low    = round(r$CI_HR_low, digits),
    CI_high   = round(r$CI_HR_hi, digits),
    p_value   = signif(r$p_value, digits),
    check.names = FALSE
  )
  names(print_df)[4:5] <- c(sprintf("%.0f%% CI_low", ci_lvl),
                            sprintf("%.0f%% CI_hi",  ci_lvl))
  print(print_df, row.names = FALSE)
  invisible(x)
}



# ============================================================================
# USAGE EXAMPLES (commented)
# ============================================================================
# source("data_generator_unified.R")   # simulation_scenario(), generate_covariates()
# source("helper_new.R")                # dedup helpers required by rd_mi()
# source("rd_mi.R")                     # rd_mi(), diagnose_rd_mi()
#
# # -- 1) No-covariate scenario (scenario 2) --
# ex  <- simulation_scenario(n_per_group = 5000, scenario = 2, seed = 1)
# dat <- ex$ADaM_df |> transmute(
#   trt01p   = TRT01PN,
#   aval     = AVAL,
#   cnsr     = CNSR,
#   pattern  = pattern,
#   RD       = RD,
#   offipday = offipday,
#   off_cnsr = off_status,
#   eos      = as.numeric(EOSDT - RANDDT)
# )
# res <- rd_mi(dat, M = 25, seed = 42)
# print(res)
#
# # -- 2) Covariates in analysis AND imputation (same covariate set) --
# beta   <- c(X11 = 0.4, X21 = 0.3, X22 = 0.7, X3B = -0.2, X3C = 0.15)
# cov_df <- generate_covariates(N = 2 * 5000)
# ex  <- simulation_scenario(n_per_group = 5000, cov_df = cov_df,
#                            beta_on = beta, seed = 1)
# dat <- ex$ADaM_df |> transmute(
#   trt01p = TRT01PN, aval = AVAL, cnsr = CNSR, pattern = pattern,
#   RD = RD, offipday = offipday, off_cnsr = off_status,
#   eos = as.numeric(EOSDT - RANDDT), X1, X2, X3
# )
# res <- rd_mi(
#   dat,
#   formula               = Surv(aval, 1 - cnsr) ~ trt01p + X1 + X2 + X3,
#   imputation_covariates = ~ X1 + X2 + X3,   # explicit: same covs as analysis
#   M                     = 25,
#   seed                  = 42
# )
# print(res)
#
# # -- 3) Different imputation vs analysis covariate sets --
# # Analysis: adjust for X1, X2, X3. Imputation model uses only X1, X2.
# res <- rd_mi(
#   dat,
#   formula               = Surv(aval, 1 - cnsr) ~ trt01p + X1 + X2 + X3,
#   imputation_covariates = ~ X1 + X2,
#   M                     = 25,
#   seed                  = 42
# )
#
# # -- 4) MCMC diagnostics for the `proper` method with covariates --
# diagnose_rd_mi(res, method = "proper", arm = "active",  plot = FALSE)
# diagnose_rd_mi(res, method = "proper", arm = "placebo", plot = TRUE)
#
# # -- 5) User-defined schema column names --
# # Time and event columns MUST be named `aval` (time) and `cnsr` (0 = event,
# # 1 = censored) . Rename if your dataset uses different names.
# # Other schema fields (pattern, RD, offipday, off_cnsr, eos) can be
# # overridden via the *_col arguments.
# my_data <- ex$ADaM_df |>
#   transmute(
#     arm = TRT01PN,
#     time = AVAL,
#     censor = CNSR,
#     trt_pattern = pattern,
#     is_retrieved = RD,
#     days_off_trt = offipday,
#     off_trt_censor = off_status,
#     study_end_day = as.numeric(EOSDT - RANDDT),
#     baseline_age = rnorm(nrow(ex$ADaM_df), mean = 50, sd = 10)
#   )
# my_dat <- my_data |>
#   rename(aval = time, cnsr = censor)        # required rename to paper convention
# res <- rd_mi(
#   my_dat,
#   formula       = Surv(aval, 1 - cnsr) ~ arm + baseline_age,
#   treatment     = "arm",
#   pattern_col   = "trt_pattern",
#   rd_col        = "is_retrieved",
#   offipday_col  = "days_off_trt",
#   off_event_col = "off_trt_censor",
#   eos_col       = "study_end_day"
# )
# print(res)
