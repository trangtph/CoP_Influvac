# P_INFL_OO6                    ################################################
# Script: Batch runner - H3N2   ################################################
# Author: Trang Tu              ################################################


## =====================================================================
## Runs the nested bootstrap + sensitivity analysis on one dataset and
## writes the pooled draws to a timestamped .csv.
## =====================================================================

# =====================================================================
# 0. Setup ============================================================
# =====================================================================

## Load packages and functions
library(pacman)
pacman::p_load(
  here,
  dplyr,
  Surrogate,
  tictoc
)

source(here("Code", "Helper_functions.R"))
source(here("Code","ICA.BinCont_full_loop.R"))
source(here("Code", "ICA.BinCont.Cov_inner.R"))

## Load data
analysis_h3n2 <- import(here("Raw_data", "analysis_h3n2.sas7bdat"))

# Convert character empty cells -> NA
analysis_h3n2 <- analysis_h3n2 |> clean_missing()

# Recode the treatment indicator

analysis_h3n2 <- analysis_h3n2 %>%
  mutate(
    ARM = factor(ARM, levels = c(1, 2), labels = c("NIV", "QIV"))         
  ) %>%
  arrange(log_hi)

pi_1. <- mean(analysis_h3n2$no_h3n2[analysis_h3n2$ARM == "NIV"])
pi_0. <- 1 - pi_1.
pi_.1  <- mean(analysis_h3n2$no_h3n2[analysis_h3n2$ARM == "QIV"])
pi_.0 = 1 - pi_.1

up_g_pi_10 <- min(pi_.0, pi_1.)

# =====================================================================
# 1. Compare models with and without covariate ========================
# =====================================================================

# Model without covariate --------------------------------------------

# Set parameters
h3_params <- data.frame(
  parameter = c(
    "mu0_00", "mu0_01", "mu0_10", "mu0_11",
    "mu1_00", "mu1_01", "mu1_10", "mu1_11"
  ),
  group = c(
    "Control", "Control", "Control", "Control",
    "QIV", "QIV", "QIV", "QIV"
  ),
  stratum = c(
    "00", "01", "10", "11",
    "00", "01", "10", "11"
  ),
  explanation = c(
    "Mean HI in control, infected under both",
    "Mean HI in control, infected under control only",
    "Mean HI in control, infected under vaccination only",
    "Mean HI in control, protected under both",
    "Mean HI in QIV, infected under both",
    "Mean HI in QIV, infected under control only",
    "Mean HI in QIV, infected under vaccination only",
    "Mean HI in QIV, protected under both"
  ),
  value = c(1.6, 2.5, 3.0, 4.5, 4.0, 5.5 , 5.5, 7)
)

mu0_h3_params <- h3_params$value[h3_params$group=="Control"]
mu1_h3_params <- h3_params$value[h3_params$group=="QIV"]

# Fit models

sigma_S0 <- 1
sigma_S1 <- 1

tic("ICA h3_sd1_no_boot")
ica_h3_sd1_no_boot <- ICA.BinCont.v2(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(mu0_h3_params, sigma_S0),
  Theta.S_1 = c(mu1_h3_params, sigma_S1),
  M = 1000,
  Seed = 2026,
  Show.Details = TRUE
)
toc()

sigma_S0 <- 1.5
sigma_S1 <- 1.5


tic("ICA h3_sd15_no_boot")
ica_h3_sd15_no_covar_no_boot <- ICA.BinCont.v2(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(mu0_h3_params, sigma_S0),
  Theta.S_1 = c(mu1_h3_params, sigma_S1),
  M = 1000,
  Seed = 2026,
  Show.Details = TRUE
)
toc()

# Abort model without covariate because of very poor convergence, model took too long to run

# =====================================================================
# 1. H3N2 Titer, no assumption ========================================
# =====================================================================

## Pre-vaccination titer as covariate

## Sigma = 0.8 ---------------------------------------------

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 0.8
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 5.5
beta_sv_S1 <- 0
sigma_S1 <- 0.8
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_s08_noboot <- ICA.BinCont.Cov(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = center_log_hi_base,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=FALSE,
  Show.Details = TRUE
)


h3n2_ICA_s08_noboot_df <- fit_to_df(h3n2_ICA_s08_noboot, boot_id = 1)
h3n2_ICA_s08_noboot_summary <- summarise_ICA(h3n2_ICA_s08_noboot_df)
summary(h3n2_ICA_s08_noboot_df$dev_S0)
summary(h3n2_ICA_s08_noboot_df$dev_S1)

summary(h3n2_ICA_s08_noboot_df$mean_Y_S0)
summary(h3n2_ICA_s08_noboot_df$mean_Y_S1)
summary(h3n2_ICA_s08_noboot_df$var_Y_S0)
summary(h3n2_ICA_s08_noboot_df$var_Y_S1)

# Dev 1196 (med 1310) - 1347

## Sigma = 1 ---------------------------------------------

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 5.5
beta_sv_S1 <- 0
sigma_S1 <- 1
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_s1_noboot <- ICA.BinCont.Cov(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = center_log_hi_base,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=FALSE,
  Show.Details = TRUE
)


h3n2_ICA_s1_noboot_df <- fit_to_df(h3n2_ICA_s1_noboot, boot_id = 1)
h3n2_ICA_s1_noboot_df <- h3n2_ICA_s1_noboot_df %>%
  mutate(
    AIC_S0 = dev_S0 + 12,
    AIC_S1 = dev_S1 + 12
  )
export(h3n2_ICA_s1_noboot_df, here("Results", "h3n2_ICA_s1_noboot_2026_07_26.csv"))

h3n2_ICA_s1_noboot_summary <- summarise_ICA(h3n2_ICA_s1_noboot_df)
summary(h3n2_ICA_s1_noboot_df$dev_S0)
summary(h3n2_ICA_s1_noboot_df$dev_S1)

summary(h3n2_ICA_s1_noboot_df$mean_Y_S0)
summary(h3n2_ICA_s1_noboot_df$mean_Y_S1)
summary(h3n2_ICA_s1_noboot_df$var_Y_S0)
summary(h3n2_ICA_s1_noboot_df$var_Y_S1)
summary(h3n2_ICA_s1_noboot_df$AIC_S0) #Med 1320
summary(h3n2_ICA_s1_noboot_df$AIC_S1) # Med 3487

# Dev 1242.2 (med 1308) - 3476

## Sigma = 1.5 -------------------------------------------------------

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 5.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_noboot <- ICA.BinCont.Cov(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = center_log_hi_base,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=FALSE,
  Show.Details = TRUE
)


h3n2_ICA_noboot_df <- fit_to_df(h3n2_ICA_noboot, boot_id = 1)
h3n2_ICA_noboot_df <- h3n2_ICA_noboot_df %>%
  mutate(
    AIC_S0 = dev_S0 + 12,
    AIC_S1 = dev_S1 + 12
         )
export(h3n2_ICA_noboot_df, here("Results", "h3n2_ICA_s15_noboot_2026_07_26.csv"))
h3n2_ICA_noboot_summary <- summarise_ICA(h3n2_ICA_noboot_df)
summary(h3n2_ICA_noboot_df$dev_S0)
summary(h3n2_ICA_noboot_df$mean_Y_S0)
summary(h3n2_ICA_noboot_df$mean_Y_S1)
summary(h3n2_ICA_noboot_df$var_Y_S0)
summary(h3n2_ICA_noboot_df$var_Y_S1)
summary(h3n2_ICA_noboot_df$AIC_S0) #Med 1321
summary(h3n2_ICA_noboot_df$AIC_S1) # Med 3488

#Dev 1206 (median 1309) - 3476



## Sigma = 1 & 2 -------------------------------------------------

## Configuration
intercept_sv_S0 <- 2.5
beta_sv_S0 <- 0
sigma_S0 <- 1
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 5.5
beta_sv_S1 <- 0
sigma_S1 <- 2
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_s1_2_noboot <- ICA.BinCont.Cov(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = center_log_hi_base,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=FALSE,
  Show.Details = TRUE
)


h3n2_ICA_s1_2_noboot_df <- fit_to_df(h3n2_ICA_s1_2_noboot, boot_id = 1)

h3n2_ICA_s1_2_noboot_df <- h3n2_ICA_s1_2_noboot_df %>%
  mutate(
    AIC_S0 = dev_S0 + 12,
    AIC_S1 = dev_S1 + 12
  )
export(h3n2_ICA_s1_2_noboot_df, here("Results", "h3n2_ICA_s1_2_noboot_2026_07_26.csv"))

h3n2_ICA_s1_2_noboot_summary <- summarise_ICA(h3n2_ICA_s1_2_noboot_df)
summary(h3n2_ICA_s1_2_noboot_df$dev_S0)
summary(h3n2_ICA_s1_2_noboot_df$mean_Y_S0)
summary(h3n2_ICA_s1_2_noboot_df$mean_Y_S1)
summary(h3n2_ICA_s1_2_noboot_df$var_Y_S0)
summary(h3n2_ICA_s1_2_noboot_df$var_Y_S1)
summary(h3n2_ICA_s1_2_noboot_df$AIC_S0) #Med 1320
summary(h3n2_ICA_s1_2_noboot_df$AIC_S1) # Med 3487

# Dev 1242.1 (median 1308) - 3475.7

## Sigma = 1.5, intercep s1 = 4.5 -------------------------------------------------

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

h3n2_ICA_is1_45_noboot <- ICA.BinCont.Cov(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = center_log_hi_base,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=FALSE,
  Show.Details = TRUE
)


h3n2_ICA_is1_45_noboot_df <- fit_to_df(h3n2_ICA_is1_45_noboot, boot_id = 1)
h3n2_ICA_is1_45_noboot_summary <- summarise_ICA(h3n2_ICA_is1_45_noboot_df)
summary(h3n2_ICA_is1_45_noboot$dev_S0)
summary(h3n2_ICA_is1_45_noboot$mean_Y_S0)
summary(h3n2_ICA_is1_45_noboot$mean_Y_S1)
summary(h3n2_ICA_is1_45_noboot$var_Y_S0)
summary(h3n2_ICA_is1_45_noboot$var_Y_S1)

## Age as covariate -----

### SD = 1.5 ----
## Configuration
intercept_sv_S0 <- 2.4
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 5.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_age_noboot <- ICA.BinCont.Cov(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = AGE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=TRUE,
  Show.Details = TRUE
)


h3n2_ICA_age_noboot_df <- fit_to_df(h3n2_ICA_age_noboot, boot_id = 1)
h3n2_ICA_age_noboot_df <- h3n2_ICA_age_noboot_df %>%
  mutate(AIC_S0 = dev_S0 + 2*6,
         AIC_S1 = dev_S1 + 1*6)
#export(h3n2_ICA_age_noboot_df, here("Results", "h3n2_ICA_age_noboot_2026_07-17.csv"))
#h3n2_ICA_age_noboot_df <- import(here("Results", "h3n2_ICA_age_noboot_2026_07-17.csv"))
h3n2_ICA_age_noboot_summary <- summarise_ICA(h3n2_ICA_age_noboot_df)

summary(h3n2_ICA_age_noboot_df$dev_S0)
summary(h3n2_ICA_age_noboot_df$mean_Y_S0)
summary(h3n2_ICA_age_noboot_df$mean_Y_S1)
summary(h3n2_ICA_age_noboot_df$var_Y_S0)
summary(h3n2_ICA_age_noboot_df$var_Y_S1)

summary(h3n2_ICA_age_noboot_df$AIC_S0) #Med 3359
summary(h3n2_ICA_age_noboot_df$AIC_S1) #Med 3857


# Dev 3346, 3844

### SD = 2 ----
intercept_sv_S0 <- 2.4
beta_sv_S0 <- 0
sigma_S0 <- 2
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 5.5
beta_sv_S1 <- 0
sigma_S1 <- 2
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_age_s2_noboot <- ICA.BinCont.Cov(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = AGE,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, sigma_S0),
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=TRUE,
  Show.Details = TRUE
)


h3n2_ICA_age_s2_noboot_df <- fit_to_df(h3n2_ICA_age_s2_noboot, boot_id = 1)
#export(h3n2_ICA_age_s2_noboot_df, here("Results", "h3n2_ICA_age_s2_noboot_2026_07-22.csv"))
#h3n2_ICA_age_s2_noboot_df <- import(here("Results", "h3n2_ICA_age_s2_noboot_2026_07-22.csv"))
h3n2_ICA_age_s2_noboot_summary <- summarise_ICA(h3n2_ICA_age_s2_noboot_df)
# Dev 3345 (med 3347) - 3844
summary(h3n2_ICA_age_s2_noboot_df$AIC_S0) #Med 3359
summary(h3n2_ICA_age_s2_noboot_df$AIC_S1) # Med 3857


## Categorical age as covariate -----
### Sigma 1.5 -----
## Configuration
intercept_sv_S0 <- 1.7
beta_sv_S0 <- 0
sigma_S0 <- 1.5
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.7
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_agecat_noboot <- ICA.BinCont.CovCat(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = age_grp,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, 0,0,0, sigma_S0), # 3 coefficents for 3 dummies
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, 0,0,0, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=FALSE,
  Show.Details = TRUE
)


h3n2_ICA_agecat_noboot_df <- fit_to_df2(h3n2_ICA_agecat_noboot, boot_id = 1)
#export(h3n2_ICA_agecat_noboot_df, here("Results", "h3n2_ICA_agecat_noboot_2026_07-15.csv"))
#h3n2_ICA_agecat_noboot_df <- import(here("Results", "h3n2_ICA_agecat_noboot_2026_07-15.csv"))
h3n2_ICA_agecat_noboot_summary <- summarise_ICA(h3n2_ICA_agecat_noboot_df)
# Dev 2120, 3847
summary(h3n2_ICA_agecat_noboot_df$AIC_S0) #Median 2136
summary(h3n2_ICA_agecat_noboot_df$AIC_S1) #Median 3863



### Sigma 1.0 -----
## Configuration
intercept_sv_S0 <- 1.7
beta_sv_S0 <- 0
sigma_S0 <- 1.0
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 4.7
beta_sv_S1 <- 0
sigma_S1 <- 1.0
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_agecat_s1_noboot <- ICA.BinCont.CovCat(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov = age_grp,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, 0,0,0, sigma_S0), # 3 coefficents for 3 dummies
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, 0,0,0, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=FALSE,
  Show.Details = TRUE
)


h3n2_ICA_agecat_s1_noboot_df <- fit_to_df2(h3n2_ICA_agecat_s1_noboot, boot_id = 1)
export(h3n2_ICA_agecat_s1_noboot_df, here("Results", "h3n2_ICA_agecat_s1_noboot_2026_07-18.csv"))
h3n2_ICA_agecat_s1_noboot_summary <- summarise_ICA(h3n2_ICA_agecat_s1_noboot_df)
summary(h3n2_ICA_agecat_s1_noboot_df$AIC_S0) #Median 2137
summary(h3n2_ICA_agecat_s1_noboot_df$AIC_S1) #Median 3863

# dev 2122 - 3847

## Both age and HI as covariates -----
## Configuration
intercept_sv_S0 <- 2.4
beta_sv_S0 <- 0
sigma_S0 <- 0.7
b1_SV_S0 <- intercept_sv_S0 - 1.5
b2_SV_S0 <- intercept_sv_S0 - 0.5
b3_SV_S0 <- intercept_sv_S0
b4_SV_S0 <- intercept_sv_S0 + 0.5

intercept_sv_S1 <- 5.5
beta_sv_S1 <- 0
sigma_S1 <- 1.5
b1_SV_S1 <- intercept_sv_S1 - 1.5
b2_SV_S1 <- intercept_sv_S1 - 0.5
b3_SV_S1 <- intercept_sv_S1
b4_SV_S1 <- intercept_sv_S1 + 0.5

## Execution

h3n2_ICA_ageHI_noboot <- ICA.BinCont.Cov2(
  Dataset = analysis_h3n2,
  Surr = log_hi,
  True = no_h3n2,             
  Treat = treat,
  Cov1 = AGE,
  Cov2 = log_hi_base,
  G_pi_10 = c(0, up_g_pi_10),
  Theta.S_0 = c(b1_SV_S0, b2_SV_S0, b3_SV_S0, b4_SV_S0, beta_sv_S0, beta_sv_S0, sigma_S0), # 2 slopes
  Theta.S_1 = c(b1_SV_S1, b2_SV_S1, b3_SV_S1, b4_SV_S1, beta_sv_S1, beta_sv_S1, sigma_S1),
  M = 1000,
  Seed = 2026,
  Center.Cov=TRUE,
  Show.Details = TRUE
)

h3n2_ICA_ageHI_noboot_df <- fit_to_df(h3n2_ICA_ageHI_noboot, boot_id = 1)
export(h3n2_ICA_ageHI_noboot_df, here("Results", "h3n2_ICA_ageHI_noboot_df_2026_07_23.csv"))
h3n2_ICA_ageHI_noboot_summary <- summarise_ICA(h3n2_ICA_ageHI_noboot_df)
summary(h3n2_ICA_ageHI_noboot_df$AIC_S0) #Median 1330.8
summary(h3n2_ICA_ageHI_noboot_df$AIC_S1) #Median 3477
# Dev 1114, 3462
