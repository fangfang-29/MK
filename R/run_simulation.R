###############################################################
## Date: 2026-07-01
## Author: Fangfang Jiang
## ------------------------------------------------------------
## Simulation study runner for RD-MI
##
## For each replication:
##   1. Generate a dataset via simulation_scenario()
##   2. Fit an unadjusted / adjusted Cox PH baseline (no imputation)
##   3. Run rd_mi() with the requested methods
##   4. Store HR, SE(logHR), and 95% CI coverage vs. reference HR
##
## Reference HR
## ------------
## Represents the treatment-policy estimand the simulation is trying to
## recover. Two ways to obtain it:
##   (a) analytical: lam_on_trt / lam_on_ctrl (only exact when no off-trt
##       contribution and no covariates -- otherwise a proxy);
##   (b) empirical:  fit Cox on a huge N reference dataset with the same
##       covariate structure but treatment discontinuation and study
##       dropout disabled -- gives the marginal treatment effect that a
##       "perfect" trial would recover.
## Users can also pass `true_hr` explicitly and skip the reference call.
##
## Depends on: data_generator.R, helper_functions.R, rd_mi.R
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

# Assumes these are sourced in the calling script:
# source("data_generator_unified.R")
# source("helper_new.R")
# source("rd_mi.R")


# ============================================================================
# .adam_to_dat(): convert simulation_scenario()$ADaM_df to rd_mi()-ready dat
# ============================================================================
.adam_to_dat <- function(adam_df, cov_names = character(0)) {
  base <- adam_df |>
    transmute(
      trt01p   = TRT01PN,
      aval     = AVAL,
      cnsr     = CNSR,
      pattern  = pattern,
      RD       = RD,
      offipday = offipday,
      off_cnsr = off_status,
      eos      = as.numeric(EOSDT - RANDDT)
    )
  if (length(cov_names) > 0) {
    base <- cbind(base, adam_df[, cov_names, drop = FALSE])
  }
  base
}


# ============================================================================
# .compute_reference_hr(): empirical HR from a huge no-dropout dataset
# ============================================================================
# Regenerates data with treatment discontinuation and study dropout disabled.
# Every subject stays on treatment until event or EOS. Fits Cox with the analysis
# formula and returns the trt coefficient's HR.
.compute_reference_hr <- function(gen_args, formula, treatment, N_ref = 1e5,
                                  cov_names = character(0), seed = NULL) {
  # Override discontinuation / dropout rates to zero
  ref_gen_args <- utils::modifyList(gen_args, list(
    lam_disc_ctrl  = 0,
    lam_disc_trt   = 0,
    p_disc_ctrl    = 0,
    p_disc_trt     = 0,
    has_off_period = FALSE
  ))
  # scenario 1..6 hard-codes these values, so we force custom mode by
  # dropping scenario
  ref_gen_args$scenario <- NULL

  ref_gen_args$n_per_group <- N_ref
  ref_gen_args$seed        <- seed

  sim <- do.call(simulation_scenario, ref_gen_args)
  dat <- .adam_to_dat(sim$ADaM_df, cov_names = cov_names)

  fit  <- coxph(formula, data = dat)
  trt_row <- which(names(coef(fit)) == treatment)
  exp(unname(coef(fit)[trt_row]))
}


# ============================================================================
# run_simulation() -- main entry point
# ============================================================================
# Arguments
# ---------
# n_sims               : replications
# n_per_group          : per-arm sample size
# gen_args             : named list of arguments forwarded to simulation_scenario().
#                        Example: list(scenario = 6)
#                                 list(scenario = NULL,
#                                      lam_on_ctrl = 0.01, lam_on_trt = 0.008,
#                                      cov_df = <df>, beta_on = <named vec>)
# formula              : analysis Cox formula, e.g. Surv(aval, 1 - cnsr) ~ trt01p
# treatment            : treatment column name (default "trt01p")
# imputation_covariates: passed to rd_mi() (default character(0) -- no cov)
# method               : one or more of "improper", "proper_like", "proper"
# M                    : imputations per replication (default 25)
# tau_cutoffs          : piecewise-exp cutoffs (default Inf)
# mcmc_control         : optional MCMC settings for `proper` method
# true_hr              : reference HR. If NULL, computed empirically via
#                        .compute_reference_hr() using `ref_hr_N`.
# ref_hr_N             : per-arm N for the reference HR calculation (default 1e5)
# seed                 : base RNG seed. Replication i uses seed + i - 1.
# progress             : print iteration messages every 5 sims
#
# Returns list with:
#   $summary   data.frame: method, true_HR, Mean_HR, Bias_pct, Emp_SE, Coverage_pct
#   $hr        matrix (n_sims x n_methods) of pooled HRs
#   $se_logHR  matrix (n_sims x n_methods) of pooled SE(logHR)
#   $coverage  matrix (n_sims x n_methods) of 95% CI coverage indicators
#   $true_hr   scalar
#   $call      match.call() for reproducibility
run_simulation <- function(
    n_sims                = 100,
    n_per_group           = 5000,
    gen_args              = list(scenario = 6),
    formula               = Surv(aval, 1 - cnsr) ~ trt01p,
    treatment             = "trt01p",
    imputation_covariates = character(0),
    method                = c("improper", "proper_like", "proper"),
    M                     = 25,
    tau_cutoffs           = Inf,
    mcmc_control          = list(),
    true_hr               = NULL,
    ref_hr_N              = 1e5,
    seed                  = 2026,
    progress              = TRUE
) {
  call_ <- match.call()
  method <- match.arg(method,
                      choices = c("improper", "proper_like", "proper"),
                      several.ok = TRUE)

  # Extract covariate column names from formula (needed for .adam_to_dat)
  rhs_vars <- all.vars(formula[[3]])
  cov_names <- setdiff(rhs_vars, treatment)

  # -----------------------------------------------------------------------
  # 1. Reference HR
  # -----------------------------------------------------------------------
  if (is.null(true_hr)) {
    if (progress) message("Computing reference HR from N=", ref_hr_N,
                          " no-dropout dataset ...")
    true_hr <- .compute_reference_hr(
      gen_args  = gen_args,
      formula   = formula,
      treatment = treatment,
      N_ref     = ref_hr_N,
      cov_names = cov_names,
      seed      = seed + 999999L
    )
    if (progress) message(sprintf("Reference HR = %.4f", true_hr))
  }

  # -----------------------------------------------------------------------
  # 2. Storage
  # -----------------------------------------------------------------------
  method_labels <- c(
    Cox         = "Cox",
    improper    = "RD-MI Improper",
    proper_like = "RD-MI Proper-like",
    proper      = "RD-MI Proper"
  )
  cols <- unname(c("Cox", method_labels[method]))

  hr_store  <- matrix(NA_real_, nrow = n_sims, ncol = length(cols),
                      dimnames = list(NULL, cols))
  se_store  <- hr_store
  cov_store <- hr_store

  # -----------------------------------------------------------------------
  # 3. Replicate
  # -----------------------------------------------------------------------
  for (i in seq_len(n_sims)) {
    if (progress && i %% 10 == 0) message("iteration ", i)
    sim_seed <- seed + i - 1L

    # Generate
    this_gen_args <- utils::modifyList(gen_args, list(
      n_per_group = n_per_group,
      seed        = sim_seed
    ))
    sim <- do.call(simulation_scenario, this_gen_args)
    dat <- .adam_to_dat(sim$ADaM_df, cov_names = cov_names)

    # ---- Method 1: Cox baseline (no imputation) ----
    cox_fit  <- coxph(formula, data = dat)
    trt_row  <- which(names(coef(cox_fit)) == treatment)
    loghr    <- unname(coef(cox_fit)[trt_row])
    se_loghr <- sqrt(vcov(cox_fit)[trt_row, trt_row])
    hr       <- exp(loghr)
    ci       <- exp(c(loghr - 1.96 * se_loghr, loghr + 1.96 * se_loghr))

    hr_store[i, "Cox"]  <- hr
    se_store[i, "Cox"]  <- se_loghr
    cov_store[i, "Cox"] <- as.integer(true_hr >= ci[1] && true_hr <= ci[2])

    # ---- Methods 2-4: RD-MI (all requested methods in one call) ----
    mi_res <- rd_mi(
      data                  = dat,
      formula               = formula,
      treatment             = treatment,
      imputation_covariates = imputation_covariates,
      method                = method,
      M                     = M,
      tau_cutoffs           = tau_cutoffs,
      mcmc_control          = mcmc_control,
      seed                  = sim_seed
    )
    for (mth in method) {
      col <- method_labels[[mth]]
      row <- mi_res$results[mi_res$results$method == mth, ]
      hr_store[i, col]  <- row$HR
      se_store[i, col]  <- row$SE_logHR
      cov_store[i, col] <- as.integer(
        true_hr >= row$CI_HR_low && true_hr <= row$CI_HR_hi
      )
    }
  }

  # -----------------------------------------------------------------------
  # 4. Aggregate
  # -----------------------------------------------------------------------
  mean_hr  <- colMeans(hr_store, na.rm = TRUE)
  bias_pct <- 100 * (mean_hr - true_hr) / true_hr
  emp_se   <- apply(hr_store, 2, sd, na.rm = TRUE)
  coverage <- 100 * colMeans(cov_store, na.rm = TRUE)

  summary_df <- data.frame(
    method                     = cols,
    true_HR                    = true_hr,
    Mean_HR                    = as.numeric(mean_hr),
    `Bias_%`                   = as.numeric(bias_pct),
    Empirical_SE               = as.numeric(emp_se),
    `Coverage_95%`             = as.numeric(coverage),
    check.names = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    summary   = summary_df,
    hr        = hr_store,
    se_logHR  = se_store,
    coverage  = cov_store,
    true_hr   = true_hr,
    n_sims    = n_sims,
    call      = call_
  )
}


# ============================================================================
# USAGE EXAMPLES
# ============================================================================
# source("data_generator_unified.R")
# source("helper_new.R")
# source("rd_mi.R")
# source("run_simulation.R")
#
#
# -----------------------------------------------------------------------
# (a) Scenario 1-6 (no covariates)
# -----------------------------------------------------------------------
# Loop over all 6 scenarios and print summary tables.
results_by_scenario <- list()
for (sc in 1:6) {
  message("\n===== Scenario ", sc, " =====")
  res <- run_simulation(
    n_sims       = 100,
    n_per_group  = 5000,
    gen_args     = list(scenario = sc),
    formula      = Surv(aval, 1 - cnsr) ~ trt01p,
    method       = c("improper", "proper_like", "proper"),
    M            = 25,
    tau_cutoffs  = Inf,
    # Analytical reference for paper scenarios: skip the empirical call
    # by supplying true_hr directly. Scenario 5 has HR = 1, scenario 6 ~= 0.8.
    true_hr      = get_true_hr_analytical(
                      lam_on_trt  = if (sc %in% c(2,4,6)) 0.008 else 0.01,
                      lam_on_ctrl = 0.01
                   ),
    seed         = 2026
  )
  print(res$summary)
  results_by_scenario[[paste0("S", sc)]] <- res
}


# -----------------------------------------------------------------------
# (b) Covariate case using scenario 5-6 style lambdas
# -----------------------------------------------------------------------
# Scenario 6 has HR ~= 0.8 with an off-treatment observation window.
# Reuse those lambdas in CUSTOM mode + attach covariates + provide beta.
# The reference HR here is EMPIRICAL (from a big no-dropout run) because
# covariates + off-trt period distort the analytical on-trt HR.

n_per_group <- 5000
cov_df      <- generate_covariates(N = 2 * n_per_group)   # default 3-cov spec
beta        <- c(X11 = 0.4, X21 = 0.3, X22 = 0.7,
                 X3B = -0.2, X3C = 0.15)

res_cov <- run_simulation(
  n_sims       = 100,
  n_per_group  = n_per_group,
  gen_args     = list(
    # Scenario 6 lambdas + p_disc + off period
    scenario       = NULL,
    lam_on_ctrl    = 0.01,   lam_off_ctrl   = 0.02,
    lam_on_trt     = 0.008,  lam_off_trt    = 0.016,
    lam_disc_ctrl  = 0.005,  lam_disc_trt   = 0.005,
    lam_cens_ctrl  = 0.02,   lam_cens_trt   = 0.06,
    p_disc_ctrl    = 0.20,   p_disc_trt     = 0.60,
    has_off_period = TRUE,
    # Covariates
    cov_df         = cov_df,
    beta_on        = beta
  ),
  formula               = Surv(aval, 1 - cnsr) ~ trt01p + X1 + X2 + X3,
  imputation_covariates = ~ X1 + X2 + X3,   # same covariates in imp model
  method                = c("improper", "proper_like", "proper"),
  M                     = 25,
  tau_cutoffs           = Inf,
  # true_hr = NULL -> empirical reference from a big no-dropout run
  ref_hr_N              = 1e5,
  seed                  = 2026
)
print(res_cov$summary)


# -----------------------------------------------------------------------
# (c) OPTIONAL: MCMC diagnostics from ONE replication of the covariate case
# -----------------------------------------------------------------------
# For the proper method with covariates, inspect ESS / Geweke / trace.
sim   <- simulation_scenario(
           n_per_group = n_per_group,
           scenario = NULL,
           lam_on_ctrl = 0.01, lam_off_ctrl = 0.02,
           lam_on_trt = 0.008, lam_off_trt = 0.016,
           lam_cens_ctrl = 0.02, lam_cens_trt = 0.06,
           p_disc_ctrl = 0.20, p_disc_trt = 0.60,
           has_off_period = TRUE,
           cov_df = cov_df, beta_on = beta,
           seed = 2026)
dat   <- .adam_to_dat(sim$ADaM_df, cov_names = c("X1","X2","X3"))
one   <- rd_mi(
           data = dat,
           formula = Surv(aval, 1 - cnsr) ~ trt01p + X1 + X2 + X3,
           imputation_covariates = ~ X1 + X2 + X3,
           method = "proper",
           M = 25,
           seed = 2026
         )
diagnose_rd_mi(one, method = "proper", arm = "active", plot = TRUE)
