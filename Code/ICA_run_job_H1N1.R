# P_INFL_OO6                    ################################################
# Script: Batch runner - H1N1   ################################################
# Author: Trang Tu              ################################################


## =====================================================================
## Runs the nested bootstrap + sensitivity analysis on one dataset and
## writes the pooled draws to a timestamped .csv.
## =====================================================================

# =====================================================================
# 0. Setup ============================================================

## Load packages and functions
if (!require("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)
pacman::p_load(
  here,
  dplyr
)

source(here("Code", "Helper_functions.R"))
source(here("Code","ICA.BinCont_full_loop.R"))
source(here("Code", "SPF.BinCont.R"))

## Load data

analysis_h1n1 <- import(here("Raw_data", "analysis_h1n1.sas7bdat"))

# Convert character empty cells -> NA
analysis_h1n1 <- analysis_h1n1 |> clean_missing()

# Recode the treatment indicator
analysis_h1n1 <- analysis_h1n1 %>%
  mutate(
    ARM = factor(ARM, levels = c(1, 2), labels = c("NIV", "QIV"))
  ) %>%
  arrange(log_hi)


pi_1. <- mean(analysis_h1n1$no_h1n1[analysis_h1n1$ARM == "NIV"])
pi_0. <- 1 - pi_1.
pi_.1  <- mean(analysis_h1n1$no_h1n1[analysis_h1n1$ARM == "QIV"])
pi_.0 = 1 - pi_.1

up_g_pi_10 <- min(pi_.0, pi_1.)

# =====================================================================
# 1. H1N1 Titer, no assumption ========================================

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_main <- run_nested_ICA(
   Dataset    = analysis_h1n1,
   surr_name  = "log_hi", 
   true_name = "no_h1n1",
   treat_name = "treat",    
   cov_name  = "center_log_hi_base",
   Center.Cov=FALSE,
   G_pi_10 = c(0, up_g_pi_10),
   Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
   Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
   inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
   x_eval = 0,
   Monotonicity = FALSE, Independence = FALSE,
   HAA = FALSE, Cond_ind = FALSE,
   Nboot = 300, 
   M = 1000, 
   global_seed = 12345,
   out_dir = "Results", 
   file_prefix = "h1n1_ICA_main", 
   save_csv = TRUE
   )

## No bootstrap version --------------------------------------------------------

h1n1_ICA_main_noboot <- ICA.BinCont.Cov(
  Dataset    = analysis_h1n1,
  Surr = log_hi,
  True = no_h1n1,             
  Treat = treat,
  Cov = center_log_hi_base,
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  x_eval = 0,
  Monotonicity = FALSE, Independence = FALSE,
  HAA = FALSE, Cond_ind = FALSE,
  M = 1000, 
  Seed = 2026,
  Show.Details = TRUE
)


# =====================================================================
# 2. H1N1 Titer, monotonicity  ========================================

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_mono <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = TRUE, Independence = FALSE,
  HAA = FALSE, Cond_ind = FALSE,
  Nboot = 300, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_mono", 
  save_csv = TRUE
)


# =====================================================================
# 3. H1N1 Titer, independence= ========================================
# Run in remote computer

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_indep <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,   
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = FALSE, Independence = TRUE,
  HAA = FALSE, Cond_ind = FALSE,
  Nboot = 300, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_indep", 
  save_csv = TRUE
)

# =====================================================================
# 4. H1N1 Titer, HVCA =================================================
# Run in remote computer

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_HVCA <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = FALSE, Independence = FALSE,
  HAA = TRUE, Cond_ind = FALSE,
  Nboot = 300, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_HVCA", 
  save_csv = TRUE
)

# =====================================================================
# 5. H1N1 Titer, Conditional independence =============================

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_Cond_ind <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = FALSE, Independence = FALSE,
  HAA = FALSE, Cond_ind = TRUE,
  Nboot = 300, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_Cond_ind", 
  save_csv = TRUE
)

# =====================================================================
# 6. H1N1 Titer, full identifiability with monotonicity ===============

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_full_id_mono <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = TRUE, Independence = FALSE,
  HAA = TRUE, Cond_ind = TRUE,
  Nboot = 300, 
  M = 1, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_full_id_mono", 
  save_csv = TRUE
)

# =====================================================================
# 7. H1N1 Titer, full identifiability with independence ===============

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_full_id_indep <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = FALSE, Independence = TRUE,
  HAA = TRUE, Cond_ind = TRUE,
  Nboot = 300, 
  M = 1, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_full_id_indep", 
  save_csv = TRUE
)

# =====================================================================
# 8. H1N1 Titer, no covariate, no assumption ==========================

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_no_covar <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = NULL,
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(1.6, 2.5, 3.0, 4.5, 1),
  Theta.S_1  = c(3.0, 4.0, 4.5, 5.5, 1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = FALSE, Independence = FALSE,
  HAA = FALSE, Cond_ind = FALSE,
  Nboot = 300, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_no_covar", 
  save_csv = TRUE
)

# No bootstrap version
h1n1_ICA_no_covar_noboot <- ICA.BinCont.v2(
  Dataset = analysis_h1n1,
  Surr = log_hi,
  True = no_h1n1,             
  Treat = treat,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(1.6, 2.5, 3.0, 4.5, 1),
  Theta.S_1  = c(3.0, 4.0, 4.5, 5.5, 1),
  M = 1000,
  Seed = 2026,
  Show.Details = TRUE
)


h1n1_ICA_no_covar_noboot_df <- fit_to_df(h1n1_ICA_no_covar_noboot, boot_id = 1)


# =====================================================================
# 9. Test SPF calculation of model 8 ==========================

# ====================================================================================
# 10. H1N1 Titer, no assumption, ICA eval at pre-vacc HI = ln(5) (1st quantile) ======

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_main_hi_q1 <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = log(5) - mean(analysis_h1n1$log_hi_base),
  Monotonicity = FALSE, Independence = FALSE,
  HAA = FALSE, Cond_ind = FALSE,
  Nboot = 300, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_main_hi_q1", 
  save_csv = TRUE
)




# ====================================================================================
# 12. H1N1 Titer, no assumption, ICA eval at pre-vacc HI = ln(20) (3rd quantile) ======

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_main_hi_q3 <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = log(20) - mean(analysis_h1n1$log_hi_base),
  Monotonicity = FALSE, Independence = FALSE,
  HAA = FALSE, Cond_ind = FALSE,
  Nboot = 300, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_main_hi_q3", 
  save_csv = TRUE
)

# =====================================================================
# 11. H1N1 Titer, no assumption, 200 bootstrap ========================

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_main <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = FALSE, Independence = FALSE,
  HAA = FALSE, Cond_ind = FALSE,
  Nboot = 200, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_main_200boot", 
  save_csv = TRUE
)

# =====================================================================
# 13. H1N1 Titer, no assumption, 500 bootstrap ========================

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h1n1_ICA_main <- run_nested_ICA(
  Dataset    = analysis_h1n1,
  surr_name  = "log_hi", 
  true_name = "no_h1n1",
  treat_name = "treat",    
  cov_name  = "center_log_hi_base",
  Center.Cov=FALSE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0  = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),   # 4 intercepts + slope + sigma
  Theta.S_1  = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
  x_eval = 0,
  Monotonicity = FALSE, Independence = FALSE,
  HAA = FALSE, Cond_ind = FALSE,
  Nboot = 500, 
  M = 1000, 
  global_seed = 12345,
  out_dir = "Results", 
  file_prefix = "h1n1_ICA_main_500boot", 
  save_csv = TRUE
)

