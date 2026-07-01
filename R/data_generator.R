###############################################################
## Date: 2026-07-01
## Author: Fangfang Jiang
## ------------------------------------------------------------
## Simulation Data Generator
## ------------------------------------------------------------
## Two modes, driven by the `scenario` argument:
##   (A) scenario in 1..6  -> paper-fixed lambdas, p_disc, and
##       has_off_period (strict paper reproduction).
##       Covariates are NOT allowed in this mode.
##   (B) scenario = NULL   -> custom mode.
##       User supplies lambdas via arguments (defaults provided).
##       Covariates optional: supply a pre-built `cov_df`
##       (via generate_covariates()) together with `beta_on`
##       (and optionally `beta_off`).
##
## Covariate structure (custom mode only): baseline proportional
## hazards on the on- and off-treatment event hazards only:
##   lam_on_i  = lam_on_base[group_i]  * exp(beta_on'  x_i)
##   lam_off_i = lam_off_base[group_i] * exp(beta_off' x_i)
## Treatment discontinuation, study dropout, and non-CV death
## remain covariate-free. Randomization is stratified by
## interaction(cov_df) when covariates are supplied.
##
## Notation:
##    - lam_on         : on-treatment hazard
##    - lam_off        : off-treatment hazard
##    - lam_disc       : treatment discontinuation rate
##    - lam_cens       : study-dropout rate (S5-6 only)
##    - lam_noncv      : non-CV death hazard (Inf = disabled)
##    - EOS            : administrative censoring, default = 100
##    - p_disc         : Pr(withdrawal at trt disc | disc'd)
##
## Output (list):
##    - input     : per-arm summary of the lambdas used
##    - df        : plain data.frame (with cov columns if used)
##    - ADaM_df   : ADaM-style data.frame (with cov columns)
##    - cov_df    : covariate frame (NULL if not used)
##    - beta_on   : on-treatment beta (NULL if not used)
##    - beta_off  : off-treatment beta (NULL if not used)
##    - scenario  : the scenario id passed in (NULL for custom)
##    - has_cov   : logical
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

# ============================================================================
# 1. generate_covariates()
# ============================================================================
# Convenience helper for building a cov_df for the custom mode.
# Arguments:
#   N     : number of subjects (must equal 2 * n_per_group when passed to
#           simulation_scenario)
#   specs : optional named list. Each element is list(levels=..., probs=...).
#           Defaults to three covariates: X1 (2-level), X2 (3-level ordinal),
#           X3 (3-level nominal).
# Returns a data.frame of factor columns.
generate_covariates <- function(N, specs = NULL) {
  if (is.null(specs)) {
    specs <- list(
      X1 = list(levels = c("0", "1"),          probs = c(0.50, 0.50)),
      X2 = list(levels = c("0", "1", "2"),     probs = c(0.40, 0.35, 0.25)),
      X3 = list(levels = c("A", "B", "C"),     probs = c(0.50, 0.30, 0.20))
    )
  }
  cov_df <- lapply(names(specs), function(nm) {
    s <- specs[[nm]]
    factor(sample(s$levels, size = N, replace = TRUE, prob = s$probs),
           levels = s$levels)
  })
  names(cov_df) <- names(specs)
  as.data.frame(cov_df, stringsAsFactors = TRUE)
}


# ============================================================================
# Helper: build linear predictor eta_i = beta' x_i
# ============================================================================
# Reference level of each factor is absorbed into the baseline hazard.
# `beta` is a named numeric vector; names must match the dummy column names
# produced by model.matrix(~ ., data = cov_df) (minus the "(Intercept)").
.build_linear_predictor <- function(cov_df, beta) {
  mm <- model.matrix(~ ., data = cov_df)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]

  b <- setNames(rep(0, ncol(mm)), colnames(mm))
  common <- intersect(names(beta), colnames(mm))
  if (length(common) == 0) {
    warning("No beta names match design-matrix columns; eta will be all zero.\n",
            "Design columns: ", paste(colnames(mm), collapse = ", "))
  }
  b[common] <- beta[common]
  as.numeric(mm %*% b)
}


# ============================================================================
# 2. simulation_scenario()
# ============================================================================
# Arguments:
#   n_per_group    : subjects per arm
#   EOS            : administrative censoring time
#   scenario       : NULL (custom) OR one of 1..6 (paper-fixed).
#                    Paper Table 1 encoding:
#                      1,3,5 -> HR = 1        ; 2,4,6 -> HR ~ 0.8
#                      1,2,5,6 -> p_disc = (0.20, 0.60)
#                      3,4     -> p_disc = (0.50, 0.90)
#                      5,6 have off-treatment observation window; others don't
#   lam_*, p_disc_*, has_off_period : custom-mode arguments; ignored (with a
#                    message) when scenario is set.
#   cov_df         : optional data.frame of factor covariates (custom mode
#                    only). Must have 2 * n_per_group rows.
#   beta_on        : named numeric vector of true betas for on-treatment
#                    hazard (required when cov_df is supplied).
#   beta_off       : named numeric vector for off-treatment hazard.
#                    NULL -> reuse beta_on (common choice).
#   seed           : RNG seed for reproducibility.
simulation_scenario <- function(
    n_per_group    = 5000,
    EOS            = 100,
    scenario       = NULL,
    # ---- custom-mode lambdas (defaults ~ scenario 2 style with off-period) ----
    lam_on_ctrl    = 0.010, lam_off_ctrl  = 0.020,
    lam_on_trt     = 0.008, lam_off_trt   = 0.016,
    lam_disc_ctrl  = 0.005, lam_disc_trt  = 0.005,
    lam_cens_ctrl  = 0.02,  lam_cens_trt  = 0.06,
    p_disc_ctrl    = 0.20,  p_disc_trt    = 0.60,
    has_off_period = TRUE,
    lam_noncv_ctrl = 0,   lam_noncv_trt = 0,
    # ---- covariates (custom mode only) ----
    cov_df   = NULL,
    beta_on  = NULL,
    beta_off = NULL,
    seed     = NULL
) {

  if (!is.null(seed)) set.seed(seed)
  N <- n_per_group * 2

  # -------------------------------------------------------------------
  # Mode dispatch: scenario 1..6 (strict, no cov) vs. custom (cov optional)
  # -------------------------------------------------------------------
  use_cov <- !is.null(cov_df)

  if (!is.null(scenario)) {
    if (!scenario %in% 1:6) {
      stop("`scenario` must be one of 1..6, or NULL for custom mode.")
    }
    if (use_cov || !is.null(beta_on) || !is.null(beta_off)) {
      stop("Covariates are not supported in scenario 1..6 mode.")
    }

    # Paper Table 1 -- lambdas and p_disc per scenario
    if (scenario %in% c(1, 3, 5)) {
      lam_on_ctrl <- 0.01; lam_off_ctrl <- 0.02
      lam_on_trt  <- 0.01; lam_off_trt  <- 0.02       # HR = 1
    } else {                                           # 2, 4, 6
      lam_on_ctrl <- 0.01; lam_off_ctrl <- 0.02
      lam_on_trt  <- 0.008; lam_off_trt <- 0.016      # HR ~ 0.8
    }
    if (scenario %in% c(1, 2, 5, 6)) {
      p_disc_ctrl <- 0.20; p_disc_trt <- 0.60
    } else {                                           # 3, 4
      p_disc_ctrl <- 0.50; p_disc_trt <- 0.90
    }
    lam_disc_ctrl <- 0.005; lam_disc_trt <- 0.005
    has_off_period <- scenario %in% c(5, 6)
    if (has_off_period) {
      lam_cens_ctrl <- 0.02; lam_cens_trt <- 0.06
    } else {
      lam_cens_ctrl <- NA;   lam_cens_trt <- NA
    }
  } else {
    # Custom mode: validate covariate arguments if cov_df is supplied.
    if (use_cov) {
      if (nrow(cov_df) != N) {
        stop("nrow(cov_df) (", nrow(cov_df), ") must equal 2 * n_per_group (",
             N, "). Build with generate_covariates(N = 2 * n_per_group).")
      }
      if (is.null(beta_on)) {
        stop("When cov_df is supplied, `beta_on` must also be provided ",
             "(named numeric vector matching model.matrix(~ ., cov_df) columns).")
      }
      if (is.null(beta_off)) beta_off <- beta_on
    }
  }

  # -------------------------------------------------------------------
  # Treatment assignment: stratified 1:1 within cov strata, else block
  # -------------------------------------------------------------------
  if (use_cov) {
    strata_id <- interaction(cov_df, drop = TRUE)
    group     <- integer(N)
    for (s in levels(strata_id)) {
      idx  <- which(strata_id == s)
      n_s  <- length(idx)
      n_t  <- floor(n_s / 2)
      # Random permutation so the leftover (if n_s odd) goes to a random arm.
      group[idx] <- sample(c(rep(1L, n_t), rep(0L, n_s - n_t)))
    }
  } else {
    # Simple 1:1 randomization: exactly n_per_group per arm, order shuffled.
    group <- sample(rep(c(0L, 1L), each = n_per_group))
  }

  # -------------------------------------------------------------------
  # Linear predictors from covariates (0 when not used)
  # -------------------------------------------------------------------
  if (use_cov) {
    eta_on  <- .build_linear_predictor(cov_df, beta_on)
    eta_off <- .build_linear_predictor(cov_df, beta_off)
  } else {
    eta_on  <- rep(0, N)
    eta_off <- rep(0, N)
  }

  # -------------------------------------------------------------------
  # Latent times
  # -------------------------------------------------------------------
  # On-treatment event time (subject-specific rate)
  lam_on_base <- ifelse(group == 1, lam_on_trt, lam_on_ctrl)
  lam_on_i    <- lam_on_base * exp(eta_on)
  T_on        <- rexp(N, rate = lam_on_i)

  # Treatment-discontinuation time (covariate-free)
  lam_disc   <- ifelse(group == 1, lam_disc_trt, lam_disc_ctrl)
  T_trt_disc <- rexp(N, rate = lam_disc)

  # Off-treatment event time (measured from T_trt_disc; subject-specific rate)
  lam_off_base <- ifelse(group == 1, lam_off_trt, lam_off_ctrl)
  lam_off_i    <- lam_off_base * exp(eta_off)
  T_off_latent <- rexp(N, rate = lam_off_i)

  # True event time: on-trt event if it precedes disc, else t_disc + off-latent
  T_event_true <- ifelse(T_on < T_trt_disc, T_on, T_trt_disc + T_off_latent)

  # Non-CV death (Inf rate = disabled)
  lam_noncv     <- ifelse(group == 1, lam_noncv_trt, lam_noncv_ctrl)
  T_noncv_death <- rep(Inf, N)
  has_noncv     <- is.finite(lam_noncv) & lam_noncv > 0
  if (any(has_noncv)) {
    T_noncv_death[has_noncv] <- rexp(sum(has_noncv), rate = lam_noncv[has_noncv])
  }

  # Study termination time via the two-stage mechanism:
  #   Stage 1: will_studydisc = TRUE subjects exit at T_trt_disc (P5)
  #   Stage 2: remaining discontinuers may exit at T_study_term if the draw
  #            falls strictly after T_trt_disc (P6); otherwise Inf (they
  #            continue until event or EOS -> P2/P4).
  disc_trt       <- T_trt_disc < pmin(T_event_true, EOS)
  p_disc         <- ifelse(group == 1, p_disc_trt, p_disc_ctrl)
  will_studydisc <- rep(FALSE, N)
  will_studydisc[disc_trt] <- runif(sum(disc_trt)) < p_disc[disc_trt]

  if (has_off_period) {
    lam_cens     <- ifelse(group == 1, lam_cens_trt, lam_cens_ctrl)
    T_study_term <- rexp(N, rate = lam_cens)
    T_study_disc <- ifelse(will_studydisc, T_trt_disc,
                           ifelse(disc_trt & T_study_term > T_trt_disc,
                                  T_study_term, Inf))
  } else {
    # Scenarios 1-4 style: study dropout coincides with treatment disc (P5 only)
    T_study_disc <- ifelse(will_studydisc, T_trt_disc, Inf)
  }

  # -------------------------------------------------------------------
  # Determine observed pattern per subject
  # -------------------------------------------------------------------
  pattern    <- character(N)
  status     <- rep(NA_integer_, N)   # 0 = event, 1 = censored
  time       <- rep(NA_real_, N)
  off_time   <- rep(NA_real_, N)      # observed off-treatment days
  off_status <- rep(NA_integer_, N)   # 0 = event, 1 = censored, -1 = N/A
  reason     <- character(N)

  for (i in seq_len(N)) {
    t_disc  <- T_trt_disc[i]; t_sdisc <- T_study_disc[i]
    t_ev    <- T_event_true[i]; t_noncv <- T_noncv_death[i]

    cands     <- c(EV = t_ev, NCVD = t_noncv, SDISC = t_sdisc, EOS = EOS)
    reason[i] <- names(which.min(cands))
    t_end     <- min(cands)
    status[i] <- if (reason[i] == "EV") 0L else 1L

    if (reason[i] == "EV") {
      # P3: event before treatment disc; P4: event after treatment disc
      pattern[i]  <- ifelse(t_ev <= t_disc, "P3", "P4")
      time[i]     <- t_end
      off_time[i] <- ifelse(pattern[i] == "P3", NA, t_ev - t_disc)
    } else if (reason[i] == "NCVD") {
      # P8: non-CV death before trt disc; P7: after trt disc
      pattern[i]  <- ifelse(t_noncv <= t_disc, "P8", "P7")
      time[i]     <- t_end
      off_time[i] <- ifelse(pattern[i] == "P8", NA, t_noncv - t_disc)
    } else if (reason[i] == "EOS") {
      # P1: never discontinued; P2: retrieved dropout who reached EOS
      pattern[i]  <- ifelse(t_disc < EOS, "P2", "P1")
      time[i]     <- EOS
      off_time[i] <- ifelse(pattern[i] == "P1", NA, t_end - t_disc)
    } else {  # SDISC
      # P5: study dropout at trt disc; P6: study dropout later
      pattern[i]  <- ifelse(t_sdisc <= t_disc, "P5", "P6")
      time[i]     <- t_end
      off_time[i] <- ifelse(pattern[i] == "P5", NA, t_end - t_disc)
    }

    if (pattern[i] %in% c("P2", "P4", "P6", "P7")) {
      off_status[i] <- if (reason[i] == "EV") 0L else 1L
    } else {
      off_status[i] <- -1L
    }
  }

  RD <- as.integer(pattern %in% c("P2", "P4", "P6", "P7"))

  # -------------------------------------------------------------------
  # Assemble outputs
  # -------------------------------------------------------------------
  sim_data <- data.frame(
    id       = seq_len(N),
    group    = group,
    pattern  = pattern,
    randdt   = 0,
    adt      = 0 + pmin(T_event_true, T_noncv_death,
                        T_trt_disc + replace(off_time, is.na(off_time), 0), EOS),
    EOS      = EOS,
    T_event_true  = T_event_true,
    T_noncv_death = T_noncv_death,
    T_trt_disc    = T_trt_disc,
    T_study_disc  = T_study_disc,
    time     = time,
    status   = status,
    RD       = RD,
    reason   = reason,
    off_time = off_time,
    off_status = off_status
  )
  if (use_cov) sim_data <- cbind(sim_data, cov_df)

  RANDDT <- as.Date("2026-01-01")
  PACUDT <- as.Date("2026-12-31")
  sim_data_ADaM <- data.frame(
    USUBJID = seq_len(N),
    TRT01PN = group,
    TRT01P  = as.character(group),
    AVAL    = time,
    CNSR    = status,
    RANDDT  = RANDDT,
    PACUDT  = PACUDT,
    ADT     = RANDDT + pmin(T_event_true, T_noncv_death,
                            T_trt_disc + replace(off_time, is.na(off_time), 0), EOS),
    EOSDT   = RANDDT + EOS,
    TRTEDT  = RANDDT + T_trt_disc,
    pattern = pattern,
    RD      = RD,
    reason  = reason,
    offipday   = off_time,
    off_status = off_status,
    T_event_true  = T_event_true,
    T_noncv_death = T_noncv_death,
    T_trt_disc    = T_trt_disc,
    T_study_disc  = T_study_disc
  )
  if (use_cov) sim_data_ADaM <- cbind(sim_data_ADaM, cov_df)

  input <- data.frame(
    group         = c("Control", "Treatment"),
    EOS           = EOS,
    lambda_on     = c(lam_on_ctrl, lam_on_trt),
    lambda_off    = c(lam_off_ctrl, lam_off_trt),
    lambda_disc   = c(lam_disc_ctrl, lam_disc_trt),
    lambda_cens   = c(lam_cens_ctrl, lam_cens_trt),
    lambda_noncv  = c(lam_noncv_ctrl, lam_noncv_trt),
    noncv_enabled = is.finite(c(lam_noncv_ctrl, lam_noncv_trt)) & (c(lam_noncv_ctrl, lam_noncv_trt) > 0),
    p_disc        = c(p_disc_ctrl, p_disc_trt)
  )

  list(
    input    = input,
    df       = sim_data,
    ADaM_df  = sim_data_ADaM,
    cov_df   = if (use_cov) cov_df else NULL,
    beta_on  = if (use_cov) beta_on else NULL,
    beta_off = if (use_cov) beta_off else NULL,
    scenario = scenario,
    has_cov  = use_cov
  )
}


# ============================================================================
# 3. diagnose_simulation()
# ============================================================================
# Verbatim from data_generator.R (cov-agnostic).
diagnose_simulation <- function(df) {

  overall_nc <- mean(df$pattern %in% c("P5", "P6"))

  nc_by_group <- df |>
    group_by(group) |>
    summarise(
      nc_pct = mean(pattern %in% c("P5", "P6")),
      .groups = "drop"
    )

  missing_by_group <- df |>
    filter(T_trt_disc < pmin(T_event_true, EOS)) |>
    group_by(group) |>
    summarise(
      missing_pct = mean(pattern %in% c("P5", "P6")),
      .groups = "drop"
    )

  HR <- exp(coef(coxph(Surv(time, status == 0) ~ group, data = df)))

  pattern_by_group <- df |>
    count(group, pattern, name = "n") |>
    group_by(group) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    mutate(group = as.character(group))

  pattern_overall <- df |>
    count(pattern, name = "n") |>
    mutate(group = "Overall", prop = n / sum(n)) |>
    dplyr::select(group, pattern, n, prop)

  pattern_tbl <- dplyr::bind_rows(pattern_by_group, pattern_overall) |>
    arrange(group, pattern) |>
    dplyr::select(group, pattern, prop) |>
    tidyr::pivot_wider(names_from = pattern, values_from = prop)

  retrieved_tbl <- df |>
    mutate(group = as.character(group), retrieved = RD) |>
    group_by(group) |>
    summarise(N = sum(retrieved), prop = mean(retrieved), .groups = "drop") |>
    bind_rows(
      df |>
        mutate(retrieved = RD) |>
        summarise(group = "Overall", N = sum(retrieved),
                  prop = mean(retrieved))
    ) |>
    mutate(value = sprintf("%d (%.1f%%)", N, prop * 100)) |>
    dplyr::select(group, value)

  cat(sprintf("Test 1-1 (Overall non-completers): %.3f\n\n", overall_nc))
  cat("Test 1-2 (Non-completers by group):\n"); print(nc_by_group); cat("\n")
  cat("Test 2 (Missing among discontinued):\n"); print(missing_by_group); cat("\n")
  cat("Test 3 (Estimated HR):\n"); print(HR); cat("\n")
  cat("Test 4-1 (Distribution of patterns):\n"); print(pattern_tbl); cat("\n")
  cat("Test 4-2 (Retrieved dropouts, RD = 1; patterns P2/4/6/7):\n")
  print(retrieved_tbl); cat("\n")

  invisible(list(
    overall_nc = overall_nc,
    nc_by_group = nc_by_group,
    missing_by_group = missing_by_group,
    HR = HR,
    pattern_tbl = pattern_tbl,
    retrieved_tbl = retrieved_tbl
  ))
}


# ============================================================================
# 4. Usage examples (commented)
# ============================================================================
# library(dplyr); library(tidyverse); library(survival)
#
# ## (A) Strict paper reproduction (no covariates) -- scenario 1
# ex <- simulation_scenario(n_per_group = 5000, scenario = 1, seed = 0)
# print(ex$input)
# diagnose_simulation(ex$df)
# head(ex$ADaM_df)
#
# ## (B) Custom mode, no covariates (defaults ~ scenario 2 with off-period)
# ex <- simulation_scenario(n_per_group = 5000, seed = 0)
# diagnose_simulation(ex$df)
#
# ## (C) Custom mode with covariates
# cov_df <- generate_covariates(N = 2 * 5000)   # default 3-cov spec
# beta   <- c(X11 = 0.40, X21 = 0.30, X22 = 0.70, X3B = -0.20, X3C = 0.15)
# ex <- simulation_scenario(
#   n_per_group = 5000,
#   cov_df      = cov_df,
#   beta_on     = beta,          # beta_off defaults to beta_on
#   seed        = 0
# )
#
# ## (D) Custom mode with covariates AND a scenario-4 style p_disc
# ex <- simulation_scenario(
#   n_per_group = 5000,
#   p_disc_ctrl = 0.50, p_disc_trt = 0.90,
#   has_off_period = FALSE,
#   cov_df = cov_df, beta_on = beta,
#   seed   = 0
# )
