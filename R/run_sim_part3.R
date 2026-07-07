###############################################################
## Date: 2026-07-07
## Author: Fangfang Jiang
## ------------------------------------------------------------
## Simulation driver for slides Part 3 (Slides 14--16).
##
## DGP:  paper Table 1 Scenario 6 lambdas + custom covariates.
##       (has_off_period = TRUE so RD-MI has a donor pool.)
##
## Two analysis conditions on the SAME simulated data:
##   (A) Unadjusted: Cox ~ trt01p ;  imputation ignores X
##   (B) Adjusted:   Cox ~ trt01p + X1 + X2 + X3 ;  imputation uses X1..X3
##
## Outputs:
##   res_unadj / res_adj         -- list from run_simulation()
##   results_slides_part3.rds    -- both objects, for reproducibility
##   results_slides_part3.csv    -- tidy long table for slide 15 fill-in
##
## Depends: data_generator_unified.R, helper_new.R, rd_mi.R, run_simulation.R
###############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})

# ---- Sources (adjust paths if this script is moved) ------------------------
source("data_generator_unified.R")
source("helper_new.R")
source("rd_mi.R")
source("run_simulation.R")

# ---- Config ---------------------------------------------------------------
n_per_group <- 5000
n_sims      <- 100
M           <- 25
seed        <- 2026

# ---- Covariate frame (fixed across reps; resampled inside run_simulation
# ----- when n_per_group differs; here sizes match so it's used verbatim) ---
set.seed(seed)
cov_df <- generate_covariates(N = 2 * n_per_group)   # default: X1, X2, X3
beta   <- c(X11 = 0.4, X21 = 0.3, X22 = 0.7,
            X3B = -0.2, X3C = 0.15)

# ---- Shared DGP: paper Scenario 6 lambdas + covariates ---------------------
gen_args_cov <- list(
  scenario       = NULL,      # custom mode (scenario 1..6 doesn't accept cov)
  lam_on_ctrl    = 0.010,     lam_off_ctrl   = 0.020,
  lam_on_trt     = 0.008,     lam_off_trt    = 0.016,     # HR ~ 0.8
  lam_disc_ctrl  = 0.005,     lam_disc_trt   = 0.005,
  lam_cens_ctrl  = 0.020,     lam_cens_trt   = 0.060,     # retrieval on
  p_disc_ctrl    = 0.20,      p_disc_trt     = 0.60,      # MNAR mechanism
  has_off_period = TRUE,
  cov_df         = cov_df,
  beta_on        = beta                                   # beta_off <- beta_on
)

# ---- (A) Unadjusted analysis ----------------------------------------------
cat("\n===== A. Unadjusted analysis =====\n")
res_unadj <- run_simulation(
  n_sims                = n_sims,
  n_per_group           = n_per_group,
  gen_args              = gen_args_cov,
  formula               = Surv(aval, 1 - cnsr) ~ trt01p,     # no cov in Cox
  imputation_covariates = character(0),                       # no cov in imp
  method                = c("improper", "proper_like", "proper"),
  M                     = M,
  tau_cutoffs           = Inf,
  true_hr               = 0.80,
  seed                  = seed
)
print(res_unadj$summary)

# ---- (B) Adjusted analysis ------------------------------------------------
cat("\n===== B. Adjusted analysis (X in both imp & Cox) =====\n")
res_adj <- run_simulation(
  n_sims                = n_sims,
  n_per_group           = n_per_group,
  gen_args              = gen_args_cov,
  formula               = Surv(aval, 1 - cnsr) ~ trt01p + X1 + X2 + X3,
  imputation_covariates = ~ X1 + X2 + X3,
  method                = c("improper", "proper_like", "proper"),
  M                     = M,
  tau_cutoffs           = Inf,
  true_hr               = 0.80,
  seed                  = seed
)
print(res_adj$summary)

# ---- Tidy combined table for Slide 15 -------------------------------------
tab_A <- res_unadj$summary; tab_A$Condition <- "A. Unadjusted"
tab_B <- res_adj$summary;   tab_B$Condition <- "B. Adjusted"

results_combined <- rbind(tab_A, tab_B)[, c("Condition", "method", "true_HR",
                                             "Mean_HR", "Bias_%",
                                             "Empirical_SE", "Coverage_95%")]
cat("\n===== Combined table (for Slide 15) =====\n")
print(results_combined, row.names = FALSE, digits = 4)

# ---- Persist for slide fill-in and forest plot ----------------------------
saveRDS(list(unadj = res_unadj, adj = res_adj, combined = results_combined),
        file = "results_slides_part3.rds")
write.csv(results_combined, "results_slides_part3.csv", row.names = FALSE)

cat("\nDone. Saved: results_slides_part3.rds, results_slides_part3.csv\n")

# ============================================================================
# Optional: quick forest-plot preview in R (matches Slide 16 layout)
# ============================================================================
# library(ggplot2)
# df <- results_combined
# df$lo <- df$Mean_HR - 1.96 * df$Empirical_SE
# df$hi <- df$Mean_HR + 1.96 * df$Empirical_SE
# df$method <- factor(df$method,
#                     levels = c("Cox", "RD-MI Improper",
#                                "RD-MI Proper-like", "RD-MI Proper"))
# ggplot(df, aes(x = Mean_HR, y = method, colour = Condition)) +
#   geom_vline(xintercept = 0.80, linetype = "dashed") +
#   geom_pointrange(aes(xmin = lo, xmax = hi),
#                   position = position_dodge(width = 0.5)) +
#   scale_x_log10() +
#   labs(x = "HR (log scale)", y = NULL,
#        title = "Simulation forest plot: unadjusted vs adjusted") +
#   theme_minimal(base_size = 12)
