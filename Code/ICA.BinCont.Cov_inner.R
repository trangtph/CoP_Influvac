# P_INFL_OO6                    ################################################
# Script: ICA_Bincont_inner     ################################################
# Author: Trang Tu              ################################################


## =====================================================================
## Extensions of function ICA.BinCont() from package `Surrogate`
## Alonso Abad et al. (Stat Med 2024;43:1083-1102)
## The following extensions were made:
## - Optimization of the original function 
## - One covariate in the component means mu_ij
## - Two covaraites
## - Categorical covariates
##
## All functions do not include bootstrap options.
## To perform bootstrap, a separte wrapper function was written in script `ICA.BinCont_full_loop.R` 
## =====================================================================

## =====================================================================
## ICA.BinCont with NO covariate in the component means mu_ij ==========
## NO-BOOTSTRAP INNER VERSION (for use inside an external bootstrap loop)
## ---------------------------------------------------------------------
##
## Model (homoscedastic case, Diff.Sigma = FALSE):
##   Control arm:       S0_k ~ sum_ij pi_ij * N( mu0_ij , sigma00 )
##   Experimental arm:  S1_k ~ sum_ij pi_ij * N( mu1_ij , sigma11 )
##
##
## STARTING VALUES (Theta.S_0, Theta.S_1):
##   Diff.Sigma = FALSE : 5 values = c(mu_00, mu_01, mu_10, mu_11, sigma)
##   Diff.Sigma = TRUE  : 8 values = c(mu_00, mu_01, mu_10, mu_11,
##                                     sigma_00, sigma_01, sigma_10, sigma_11)
##

## DIFFERENCES in workflow FROM ICA.BinCont.R
## (A) BOOTSTRAP REMOVED. The internal BS argument and both BS branches are
##     gone. This function ALWAYS operates on the dataset it is given. Draw
##     bootstrap samples in an OUTER loop and call this function on each
##     resampled dataset.
##     Because the data are fixed within a call, the per-arm surrogate and
##     covariate vectors and the marginal probabilities are computed ONCE,
##     before the sensitivity loop.
##
## (B) SEED SCHEME. Two kinds of randomness are now handled separately:
##
##     1. SCIENTIFIC randomness (the runif draws of pi_10 and the four rho's).
##        Governed by one stream. The stream is set once and then flows through
##        the loop. Control it via `Seed`:
##          * Seed = <integer> : reproducible standalone run (set.seed(Seed)).
##          * Seed = NULL       : do not reseed. Use the RNG stream already in
##                                effect, e.g. the per-task L'Ecuyer-CMRG stream
##                                that future/%dofuture% (seed = TRUE) assigns to
##                                each bootstrap worker in . 
##                                Use this in the outer loop so each bootstrap gets an 
##                                independent,reproducible stream that this
##                                function does not overwrite.
##
##     2. MONTE CARLO INTEGRATION noise (the rnorm samples for I_00, I_01, I_10,
##        I_11). The original reset set.seed(1..4) every iteration so this noise
##        was identical across iterations (a common-random-numbers device). The
##        identical pattern is now drawn once, before the loop, as four fixed
##        N(0,1) vectors z_00, z_01, z_10, z_11 (seeds 1..4), and reused inside
##        the loop via the location-scale identity
##            rnorm(n, mu, sd) == mu + sd * z,   z ~ N(0,1).
##        The changing per-iteration mu and sd are applied to these fixed
##        patterns, so the integral estimates are bit-for-bit the same as the
##        original while no RNG is touched inside the loop. The four patterns
##        are generated without disturbing the scientific stream (the stream
##        state is saved before and restored after).
## =====================================================================

ICA.BinCont.v2 <- function(Dataset, Surr, True, Treat,
                           G_pi_10 = c(0,1),
                           G_rho_01_00=c(-1,1),
                           G_rho_01_01=c(-1,1),
                           G_rho_01_10=c(-1,1),
                           G_rho_01_11=c(-1,1),
                           Theta.S_0, # 4 means + sigma(s)
                           Theta.S_1,
                           M=1000,        # number of sensitivity-analysis draws
                           Seed=123,      # integer for standalone use; NULL inside an outer parallel loop
                           N.MC=10000,    # Monte Carlo size for the integrals I_ij
                           Monotonicity=FALSE,
                           Independence=FALSE,
                           HAA=FALSE,
                           Cond_ind=FALSE,
                           Plots=TRUE, Save.Plots="No",
                           Show.Details=FALSE){
  # 1. PREPARATION --------------------------------------------------------------

  Diff.Sigma = FALSE # All 4 SD values in each group are the same
  
  ## Extract variables
  Surr  <- Dataset[,paste(substitute(Surr))]
  True  <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  
  ## Check data format
  data_no_miss.original <- data.frame(na.exclude(cbind(Surr, True, Treat)))
  if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (max(na.exclude(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
  if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  
  ## CHECK LENGTH OF STARTING VALUES 
  need <- if (Diff.Sigma) 8 else 5
  if (length(Theta.S_0)!=need | length(Theta.S_1)!=need){
    stop(paste0("\nWithout a covariate, Theta.S_0 and Theta.S_1 must have ",
                need, " values (4 means",
                if (Diff.Sigma) " + 4 sigmas)." else " + 1 sigma)."))
  }
  
  ## SEED SETUP 
  # Scientific stream: set once here, then left to flow through the loop.
  if (!is.null(Seed)) set.seed(Seed)
  # Make sure a random seed already exist, if not create one
  if (!exists(".Random.seed", envir = .GlobalEnv)) stats::runif(1)
  # Store the current state of the random seed to use later
  .sci_state <- get(".Random.seed", envir = .GlobalEnv)
  
  # Fixed N(0,1) integration patterns (reproduce the original seeds 1..4),
  # drawn once and held constant for every iteration. The scientific stream
  # saved above is restored immediately after, so these draws do not perturb it.
  set.seed(1L); z_10 <- rnorm(N.MC)
  set.seed(2L); z_01 <- rnorm(N.MC)
  set.seed(3L); z_00 <- rnorm(N.MC)
  set.seed(4L); z_11 <- rnorm(N.MC)
  
  assign(".Random.seed", .sci_state, envir = .GlobalEnv)
  
  ## Functions to calculate mixture deviances -----------------------------
  # The mixing weights pi_ij are passed in as arguments (p00..p11).
  if (Diff.Sigma==TRUE){
    # Unequal sigma's: theta = c(mu00,mu01,mu10,mu11, s00,s01,s10,s11)
    mixt.deviance <- function(theta, data, p00, p01, p10, p11){
      # Specify the density function of normal mixture models
      pdf <- p00*dnorm(data, theta[1], theta[5]) +
        p01*dnorm(data, theta[2], theta[6]) +
        p10*dnorm(data, theta[3], theta[7]) +
        p11*dnorm(data, theta[4], theta[8])
      # Deviance = -2 log-likelihood
      -2*sum(log(pdf))
    }
  } else {
    ## Equal sigma's: theta = c(mu00,mu01,mu10,mu11, sigma)
    mixt.deviance <- function(theta, data, p00, p01, p10, p11){
      sigma_all <- theta[5]
      pdf <- p00*dnorm(data, theta[1], sigma_all) +
        p01*dnorm(data, theta[2], sigma_all) +
        p10*dnorm(data, theta[3], sigma_all) +
        p11*dnorm(data, theta[4], sigma_all)
      -2*sum(log(pdf))
    }
  }
  
  ## PRE-ALLOCATE RESULT ACCUMULATORS
  R2_H_all <- R2_C_all <- numeric(M)
  pi_00_all <- pi_01_all <- pi_10_all <- pi_11_all <- numeric(M)
  G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- numeric(M)
  pi_Delta_T_min1_all <- pi_Delta_T_0_all <- pi_Delta_T_1_all <- numeric(M)
  pi_0_00_e_all <- pi_0_01_e_all <- pi_0_10_e_all <- pi_0_11_e_all <- numeric(M)
  mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- numeric(M)
  sigma_00_00_all <- sigma_00_01_all <- sigma_00_10_all <- sigma_00_11_all <- numeric(M)
  pi_1_00_e_all <- pi_1_01_e_all <- pi_1_10_e_all <- pi_1_11_e_all <- numeric(M)
  mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- numeric(M)
  sigma_11_00_all <- sigma_11_01_all <- sigma_11_10_all <- sigma_11_11_all <- numeric(M)
  mean_Y_S0_all <- mean_Y_S1_all <- numeric(M)
  var_Y_S0_all <- var_Y_S1_all <- numeric(M)
  dev_S0_all <- dev_S1_all <- numeric(M)
  AIC_S0_all <- AIC_S1_all <- AIC_tot_all <- numeric(M)   # AIC per arm and combined
  code_nlm_0_all <- code_nlm_1_all <- numeric(M)
  aantal <- 0
  
  ## Per-arm vectors and marginal probabilities ------------------------
  dnm <- data_no_miss.original
  data_conttreat <- dnm[dnm$Treat == -1, , drop=FALSE]   # control arm
  data_exptreat  <- dnm[dnm$Treat ==  1, , drop=FALSE]   # experimental arm
  S_0 <- data_conttreat$Surr
  S_1 <- data_exptreat$Surr
  
  # Estimate the marginal probabilities of T
  pi_punt1 <- mean(data_exptreat$True)  # E(T|Z=1)
  pi_punt0 <- 1 - pi_punt1
  pi_1punt <- mean(data_conttreat$True) # E(T|Z=0)
  pi_0punt <- 1 - pi_1punt
  if (Monotonicity==TRUE){
    if (pi_punt0 > pi_0punt | pi_1punt > pi_punt1) {stop("\nThe monotonicity assumption is not compatible with the data")}
    G_pi_10 = c(0,0)
  }
  
  # 2. SENSITIVITY ANALYSIS ----------------------------------------------------

  
  while (aantal < M){
    
    ## 2.1. Sample the unidentifiable parameters -------------------------------
    ## (no per-iteration set.seed; the stream set up above advances naturally)
    if(HAA==FALSE){
      G_rho_01_00_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
      G_rho_01_01_hier <- runif(1, min=G_rho_01_01[1], max=G_rho_01_01[2])
      G_rho_01_10_hier <- runif(1, min=G_rho_01_10[1], max=G_rho_01_10[2])
      G_rho_01_11_hier <- runif(1, min=G_rho_01_11[1], max=G_rho_01_11[2])
    }
    if(HAA==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
    }
    if(Cond_ind==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- 0
    }
    if(Show.Details==TRUE){
      flush.console()
      cat("\n\nG_rho:", cbind(G_rho_01_00_hier, G_rho_01_01_hier, G_rho_01_10_hier, G_rho_01_11_hier), "\n")
    }
    
    # Sample pi_10 and derive the remaining pi_ij compatible with the margins
    if(Independence==FALSE){
      pi_10_hier <- runif(1, min=G_pi_10[1], max=G_pi_10[2])
      pi_00_hier <- pi_punt0 - pi_10_hier
      pi_01_hier <- pi_0punt - pi_00_hier
      pi_11_hier <- pi_1punt - pi_10_hier
    }
    if(Independence==TRUE){
      pi_10_hier <- pi_1punt*pi_punt0
      pi_00_hier <- pi_0punt*pi_punt0
      pi_01_hier <- pi_0punt*pi_punt1
      pi_11_hier <- pi_1punt*pi_punt1
    }
    
    if ((pi_00_hier >= 0 & pi_01_hier >= 0 & pi_10_hier >= 0 & pi_11_hier >= 0) &
        (pi_00_hier <= 1 & pi_01_hier <= 1 & pi_10_hier <= 1 & pi_11_hier <= 1)){
      
      if(Show.Details==TRUE){
        flush.console()
        cat("\n\nPi's used here:", cbind(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), "\n")
      }
      
      ## 2.2. Calculate pi_delta_T -----------------------------------------------
      # individual causal effect of the treatment on the true endpoint: delta_T = T1 - T0
      pi_Delta_T_min1 <- pi_10_hier # Delta_T = -1
      pi_Delta_T_0    <- pi_00_hier + pi_11_hier # Delta_T = 0
      pi_Delta_T_1    <- pi_01_hier # Delta_T = 1
      
      ## 2.3. Estimate Gaussian mixture models of Y_S| Y_T -----------------------
      
      # Fit model: maximize likelihood = minimizing deviance
      # Iterative process: start with an initial guess of parameter values,
      # then iteratively updat them to minimize deviance, using Newton-Raphson type algorithm
      
      options(warn = -999)
      mix1 <- nlm(f = mixt.deviance, p = Theta.S_0, data = S_0,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000 #the maximum number of iterations to be performed before the program is terminated
      )
      mix2 <- nlm(f = mixt.deviance, p = Theta.S_1, data = S_1,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000)
      options(warn = 1)
      
      if (exists("mix1")==TRUE & exists("mix2")==TRUE & is.null(mix1)==FALSE & is.null(mix2)==FALSE &
          mix1$code<=3 & mix2$code<=3){
        
        ## 2.4. Extract fitted parameters (means, residual variances) & deviance ----
        # Control arm (S0)
        mu_0_00 <- mix1$estimate[1]; mu_0_01 <- mix1$estimate[2]
        mu_0_10 <- mix1$estimate[3]; mu_0_11 <- mix1$estimate[4]
        if (Diff.Sigma==FALSE){ sigma_00_00 <- sigma_00_01 <- sigma_00_10 <- sigma_00_11 <- mix1$estimate[5]**2 }
        if (Diff.Sigma==TRUE){
          sigma_00_00 <- mix1$estimate[5]**2; sigma_00_01 <- mix1$estimate[6]**2
          sigma_00_10 <- mix1$estimate[7]**2; sigma_00_11 <- mix1$estimate[8]**2 }
        pi_0_00_e <- pi_00_hier; pi_0_01_e <- pi_01_hier; pi_0_10_e <- pi_10_hier; pi_0_11_e <- pi_11_hier
        code_nlm_0 <- mix1$code
        
        # Experimental arm (S1)
        mu_1_00 <- mix2$estimate[1]; mu_1_01 <- mix2$estimate[2]
        mu_1_10 <- mix2$estimate[3]; mu_1_11 <- mix2$estimate[4]
        if (Diff.Sigma==FALSE){ sigma_11_00 <- sigma_11_01 <- sigma_11_10 <- sigma_11_11 <- mix2$estimate[5]**2 }
        if (Diff.Sigma==TRUE){
          sigma_11_00 <- mix2$estimate[5]**2; sigma_11_01 <- mix2$estimate[6]**2
          sigma_11_10 <- mix2$estimate[7]**2; sigma_11_11 <- mix2$estimate[8]**2 }
        pi_1_00_e <- pi_00_hier; pi_1_01_e <- pi_01_hier; pi_1_10_e <- pi_10_hier; pi_1_11_e <- pi_11_hier
        code_nlm_1 <- mix2$code
        
        # Conditional deviance of each fitted mixture (for GOF check)
        dev_S0 <- mix1$minimum
        dev_S1 <- mix2$minimum
        
        # AIC = deviance + 2*k. k = number of parameters ESTIMATED per arm
        # (4 intercepts + slope + sigma(s) = length of Theta.S, i.e. `need`).
        # The mixing weights pi_ij are fixed inputs of the sensitivity draw,
        # not estimated by MLE, so they are NOT counted in k. Each arm is a
        # separate fit; AIC_tot is the joint criterion for the two independent
        # arms (deviances add, parameters add).
        k_arm   <- need
        AIC_S0  <- dev_S0 + 2*k_arm
        AIC_S1  <- dev_S1 + 2*k_arm
        AIC_tot <- (dev_S0 + dev_S1) + 2*(2*k_arm)
        
        if (Plots == TRUE){
          ## Fitted marginal density overlaid on the marginal histogram of S.
          sd_0 <- sqrt(sigma_00_00); sd_1 <- sqrt(sigma_11_00)
          x_S_0 <- sort(S_0)
          S0_M1 <- pi_00_hier * dnorm(x_S_0, mu_0_00, sd_0)
          S0_M2 <- pi_01_hier * dnorm(x_S_0, mu_0_01, sd_0)
          S0_M3 <- pi_10_hier * dnorm(x_S_0, mu_0_10, sd_0)
          S0_M4 <- pi_11_hier * dnorm(x_S_0, mu_0_11, sd_0)
          Sum_M_S_0 <- S0_M1 + S0_M2 + S0_M3 + S0_M4
          x_S_1 <- sort(S_1)
          S1_M1 <- pi_00_hier * dnorm(x_S_1, mu_1_00, sd_1)
          S1_M2 <- pi_01_hier * dnorm(x_S_1, mu_1_01, sd_1)
          S1_M3 <- pi_10_hier * dnorm(x_S_1, mu_1_10, sd_1)
          S1_M4 <- pi_11_hier * dnorm(x_S_1, mu_1_11, sd_1)
          Sum_M_S_1 <- S1_M1 + S1_M2 + S1_M3 + S1_M4
          
          plot_one <- function(){
            par(mfrow=c(1,2))
            max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
            hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]),
                 main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_0, Sum_M_S_0, col="red", lwd=4, lty=1)
            lines(x_S_0, S0_M1, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M2, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M3, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M4, col="blue", lwd=2, lty=3)
            max_val <- max(max(density(S_1)$y), max(Sum_M_S_1))
            hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]),
                 main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_1, Sum_M_S_1, col="red", lwd=4, lty=1)
            lines(x_S_1, S1_M1, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M2, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M3, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M4, col="blue", lwd=2, lty=3)
            par(mfrow=c(1,1))
          }
          if (Save.Plots!="No"){
            options(warn=-1); dir.create(Save.Plots); options(warn=0)
            pdf(paste(Save.Plots, "/Run", (aantal+1), ".pdf", sep=""),width=6,height=4,paper='special')
            plot_one(); dev.off()
          }
          plot_one()
        }
        
        ## 2.5 Covariance-matrix checks (PD and condition number) --------------
        # If checks do not pass, the results are discared (so we always have M valid results)
        mat_a <- matrix(c(sigma_00_00,
                          (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                          (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                          sigma_11_00), nrow = 2, byrow = TRUE)
        eigen_mat_a <- min(eigen((mat_a))$values)
        mat_b <- matrix(c(sigma_00_01,
                          (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                          (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                          sigma_11_01), nrow = 2, byrow = TRUE)
        eigen_mat_b <- min(eigen((mat_b))$values)
        mat_c <- matrix(c(sigma_00_10,
                          (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                          (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                          sigma_11_10), nrow = 2, byrow = TRUE)
        eigen_mat_c <- min(eigen((mat_c))$values)
        mat_d <- matrix(c(sigma_00_11,
                          (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                          (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                          sigma_11_11), nrow = 2, byrow = TRUE)
        eigen_mat_d <- min(eigen((mat_d))$values)
        
        CN_mat_a <- max(svd(mat_a)$d)/min(svd(mat_a)$d)
        CN_mat_b <- max(svd(mat_b)$d)/min(svd(mat_b)$d)
        CN_mat_c <- max(svd(mat_c)$d)/min(svd(mat_c)$d)
        CN_mat_d <- max(svd(mat_d)$d)/min(svd(mat_d)$d)
        
        if (eigen_mat_a > 0 & eigen_mat_b > 0 & eigen_mat_c > 0 & eigen_mat_d > 0 &
            CN_mat_a < 50 & CN_mat_b < 50 & CN_mat_c < 50 & CN_mat_d < 50){
          
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nMin Eigenvalues:", cbind(eigen_mat_a, eigen_mat_b, eigen_mat_c, eigen_mat_d), "\n")
            cat("\n\nCondition number:", cbind(CN_mat_a, CN_mat_b, CN_mat_c, CN_mat_d), "\n")
            cat("\n\nChecks passed...")
          }
          
          ## 2.6 Marginal moments of S0, S1 (for GOF check) --------------------
          mean_int_S0 <- (pi_0_00_e*mu_0_00) + (pi_0_01_e*mu_0_01) + (pi_0_10_e*mu_0_10) + (pi_0_11_e*mu_0_11) # Equation 2 in paper 
          mean_int_S1 <- (pi_1_00_e*mu_1_00) + (pi_1_01_e*mu_1_01) + (pi_1_10_e*mu_1_10) + (pi_1_11_e*mu_1_11)
          
          mean_Y_S0 <- mean_int_S0   # E(S | Z=0)
          mean_Y_S1 <- mean_int_S1   # E(S | Z=1)
          
          var_Y_S0 <- ((pi_0_00_e*mu_0_00^2)+(pi_0_01_e*mu_0_01^2)+(pi_0_10_e*mu_0_10^2)+(pi_0_11_e*mu_0_11^2)) - mean_int_S0^2 +
            ((pi_0_00_e*sigma_00_00)+(pi_0_01_e*sigma_00_01)+(pi_0_10_e*sigma_00_10)+(pi_0_11_e*sigma_00_11)) # Equation 3 in paper
          
          var_Y_S1 <- ((pi_1_00_e*mu_1_00^2)+(pi_1_01_e*mu_1_01^2)+(pi_1_10_e*mu_1_10^2)+(pi_1_11_e*mu_1_11^2)) - mean_int_S1^2 +
            ((pi_1_00_e*sigma_11_00)+(pi_1_01_e*sigma_11_01)+(pi_1_10_e*sigma_11_10)+(pi_1_11_e*sigma_11_11))
          
          ## 2.7 Calculate ICA via Monte Carlo ---------------------------------
          ## Treatment effect on the surrogate within each T-stratum. (Equation 7 in paper)
          ## mu*_ij = mu1_ij - mu0_ij
          mu_s00 <- mu_1_00 - mu_0_00
          mu_s01 <- mu_1_01 - mu_0_01
          mu_s10 <- mu_1_10 - mu_0_10
          mu_s11 <- mu_1_11 - mu_0_11
          
          # Weights of pi_00 and pi_11
          omega_1 <- pi_00_hier / (pi_00_hier + pi_11_hier); omega_2 <- 1 - omega_1
          ## Equation 8: Variance of Delta S in each subgroup ----
          sigma_s00 <- sigma_00_00 + sigma_11_00 - (2 * (sqrt(sigma_00_00 * sigma_11_00) * G_rho_01_00_hier))
          sigma_s01 <- sigma_00_01 + sigma_11_01 - (2 * (sqrt(sigma_00_01 * sigma_11_01) * G_rho_01_01_hier))
          sigma_s10 <- sigma_00_10 + sigma_11_10 - (2 * (sqrt(sigma_00_10 * sigma_11_10) * G_rho_01_10_hier))
          sigma_s11 <- sigma_00_11 + sigma_11_11 - (2 * (sqrt(sigma_00_11 * sigma_11_11) * G_rho_01_11_hier))
          
          f_Delta_S <- function(val){
            (pi_00_hier * dnorm(val, mu_s00, sd = sqrt(sigma_s00))) +
              (pi_01_hier * dnorm(val, mu_s01, sd = sqrt(sigma_s01))) +
              (pi_10_hier * dnorm(val, mu_s10, sd = sqrt(sigma_s10))) +
              (pi_11_hier * dnorm(val, mu_s11, sd = sqrt(sigma_s11)))
          }
          
          ## Sample S1, S2, S3, S4 from normal distribution
          ## S = mu + sd * z  is identical to rnorm(n, mu, sd) under the same seed.
          ## z values were sampled from N(0,1) on line 127-130
          S1mc <- mu_s10 + sqrt(sigma_s10) * z_10
          I_10 <- mean(log(dnorm(S1mc, mu_s10, sd = sqrt(sigma_s10)) / f_Delta_S(S1mc)))
          
          S2mc <- mu_s01 + sqrt(sigma_s01) * z_01
          I_01 <- mean(log(dnorm(S2mc, mu_s01, sd = sqrt(sigma_s01)) / f_Delta_S(S2mc)))
          
          S3mc <- mu_s00 + sqrt(sigma_s00) * z_00
          I_00 <- mean(log(((omega_1 * dnorm(S3mc, mu_s00, sd = sqrt(sigma_s00))) +
                              (omega_2 * dnorm(S3mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S3mc)))
          
          S4mc <- mu_s11 + sqrt(sigma_s11) * z_11
          I_11 <- mean(log(((omega_1 * dnorm(S4mc, mu_s00, sd = sqrt(sigma_s00))) +
                              (omega_2 * dnorm(S4mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S4mc)))
          
          # Calculate ICA
          I_Delta_T_Delta_S <- (pi_00_hier * I_00) + (pi_01_hier * I_01) + (pi_10_hier * I_10) + (pi_11_hier * I_11)
          H_Delta_T <- -((pi_Delta_T_min1 * log(pi_Delta_T_min1)) +
                           (pi_Delta_T_0 * log(pi_Delta_T_0)) +
                           (pi_Delta_T_1 * log(pi_Delta_T_1)))
          if(Monotonicity==TRUE){
            H_Delta_T <- -((pi_Delta_T_0 * log(pi_Delta_T_0)) +
                             (pi_Delta_T_1 * log(pi_Delta_T_1)))
          }
          R2_H <- I_Delta_T_Delta_S / H_Delta_T
          R2_C <- 1 - exp(-2*I_Delta_T_Delta_S)
          
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nR2H:", R2_H, "  R2C:", R2_C, "\n")
          }
          
          aantal <- aantal + 1
          if(Show.Details==TRUE){
            cat("\n", (aantal/M)*100, "% done. \n", sep="")
          }
          
          # save results (indexed assignment into pre-allocated vectors)
          R2_H_all[aantal] <- R2_H
          R2_C_all[aantal] <- R2_C
          G_rho_01_00_all[aantal] <- G_rho_01_00_hier
          G_rho_01_01_all[aantal] <- G_rho_01_01_hier
          G_rho_01_10_all[aantal] <- G_rho_01_10_hier
          G_rho_01_11_all[aantal] <- G_rho_01_11_hier
          pi_00_all[aantal] <- pi_00_hier
          pi_01_all[aantal] <- pi_01_hier
          pi_10_all[aantal] <- pi_10_hier
          pi_11_all[aantal] <- pi_11_hier
          pi_Delta_T_min1_all[aantal] <- pi_Delta_T_min1
          pi_Delta_T_0_all[aantal] <- pi_Delta_T_0
          pi_Delta_T_1_all[aantal] <- pi_Delta_T_1
          pi_0_00_e_all[aantal] <- pi_0_00_e
          pi_0_01_e_all[aantal] <- pi_0_01_e
          pi_0_10_e_all[aantal] <- pi_0_10_e
          pi_0_11_e_all[aantal] <- pi_0_11_e
          mu_0_00_all[aantal] <- mu_0_00
          mu_0_01_all[aantal] <- mu_0_01
          mu_0_10_all[aantal] <- mu_0_10
          mu_0_11_all[aantal] <- mu_0_11
          sigma_00_00_all[aantal] <- sigma_00_00
          sigma_00_01_all[aantal] <- sigma_00_01
          sigma_00_10_all[aantal] <- sigma_00_10
          sigma_00_11_all[aantal] <- sigma_00_11
          pi_1_00_e_all[aantal] <- pi_1_00_e
          pi_1_01_e_all[aantal] <- pi_1_01_e
          pi_1_10_e_all[aantal] <- pi_1_10_e
          pi_1_11_e_all[aantal] <- pi_1_11_e
          mu_1_00_all[aantal] <- mu_1_00
          mu_1_01_all[aantal] <- mu_1_01
          mu_1_10_all[aantal] <- mu_1_10
          mu_1_11_all[aantal] <- mu_1_11
          sigma_11_00_all[aantal] <- sigma_11_00
          sigma_11_01_all[aantal] <- sigma_11_01
          sigma_11_10_all[aantal] <- sigma_11_10
          sigma_11_11_all[aantal] <- sigma_11_11
          mean_Y_S0_all[aantal] <- mean_Y_S0
          mean_Y_S1_all[aantal] <- mean_Y_S1
          var_Y_S0_all[aantal] <- var_Y_S0
          var_Y_S1_all[aantal] <- var_Y_S1
          dev_S0_all[aantal] <- dev_S0
          dev_S1_all[aantal] <- dev_S1
          AIC_S0_all[aantal] <- AIC_S0
          AIC_S1_all[aantal] <- AIC_S1
          AIC_tot_all[aantal] <- AIC_tot
          code_nlm_0_all[aantal] <- code_nlm_0
          code_nlm_1_all[aantal] <- code_nlm_1
          
        } # min eigen and CN requirements
      } # mixture distribution available
    } # pi's requirements
  } # end while loop
  
  # 3. OUTPUT RESULTS ----------------------------------------------------
  fit <- list(R2_H=R2_H_all, R2_C=R2_C_all,
              pi_00=as.numeric(pi_00_all), pi_01=as.numeric(pi_01_all),
              pi_10=as.numeric(pi_10_all), pi_11=as.numeric(pi_11_all),
              G_rho_01_00=as.numeric(G_rho_01_00_all),
              G_rho_01_01=as.numeric(G_rho_01_01_all),
              G_rho_01_10=as.numeric(G_rho_01_10_all),
              G_rho_01_11=as.numeric(G_rho_01_11_all),
              pi_Delta_T_min1=as.numeric(pi_Delta_T_min1_all),
              pi_Delta_T_0=as.numeric(pi_Delta_T_0_all),
              pi_Delta_T_1=as.numeric(pi_Delta_T_1_all),
              pi_0_00=as.numeric(pi_0_00_e_all), pi_0_01=as.numeric(pi_0_01_e_all),
              pi_0_10=as.numeric(pi_0_10_e_all), pi_0_11=as.numeric(pi_0_11_e_all),
              mu_0_00=as.numeric(mu_0_00_all), mu_0_01=as.numeric(mu_0_01_all),
              mu_0_10=as.numeric(mu_0_10_all), mu_0_11=as.numeric(mu_0_11_all),
              sigma2_00_00=as.numeric(sigma_00_00_all), sigma2_00_01=as.numeric(sigma_00_01_all),
              sigma2_00_10=as.numeric(sigma_00_10_all), sigma2_00_11=as.numeric(sigma_00_11_all),
              pi_1_00=as.numeric(pi_1_00_e_all), pi_1_01=as.numeric(pi_1_01_e_all),
              pi_1_10=as.numeric(pi_1_10_e_all), pi_1_11=as.numeric(pi_1_11_e_all),
              mu_1_00=as.numeric(mu_1_00_all), mu_1_01=as.numeric(mu_1_01_all),
              mu_1_10=as.numeric(mu_1_10_all), mu_1_11=as.numeric(mu_1_11_all),
              sigma2_11_00=as.numeric(sigma_11_00_all), sigma2_11_01=as.numeric(sigma_11_01_all),
              sigma2_11_10=as.numeric(sigma_11_10_all), sigma2_11_11=as.numeric(sigma_11_11_all),
              mean_Y_S0=as.numeric(mean_Y_S0_all), mean_Y_S1=as.numeric(mean_Y_S1_all),
              var_Y_S0=as.numeric(var_Y_S0_all), var_Y_S1=as.numeric(var_Y_S1_all),
              dev_S0=as.numeric(dev_S0_all), dev_S1=as.numeric(dev_S1_all),
              AIC_S0=as.numeric(AIC_S0_all), AIC_S1=as.numeric(AIC_S1_all),
              AIC_tot=as.numeric(AIC_tot_all),       # joint AIC of the two arms
              npar=need,                             # estimated parameters PER arm (k)
              code_nlm_0=as.numeric(code_nlm_0_all), code_nlm_1=as.numeric(code_nlm_1_all),
              mean.S0=mean(S_0, na.rm=TRUE), var.S0=var(S_0, na.rm=TRUE),
              mean.S1=mean(S_1, na.rm=TRUE), var.S1=var(S_1, na.rm=TRUE),
              Call=match.call())
  
  class(fit) <- "ICA.BinCont"
  fit
} # end of function

## =====================================================================
## ICA.BinCont with ONE covariate in the component means mu_ij
## NO-BOOTSTRAP INNER VERSION (for use inside an external bootstrap loop)
## ---------------------------------------------------------------------
## Extension of Alonso Abad et al. (Stat Med 2024;43:1083-1102) ICA.BinCont().
##
## Model for distribution of f(Y_S |X) (homoscedastic case):
##   Control arm:       S0_k ~ sum_ij pi_ij * N( beta0_ij + gamma0 * X_k , sigma00 )
##   Experimental arm:  S1_k ~ sum_ij pi_ij * N( beta1_ij + gamma1 * X_k , sigma11 )
##   where X_k = (pre-vaccination titer)_k, centered at its sample mean.
##
## This is a finite mixture of LINEAR REGRESSIONS with KNOWN mixing weights pi_ij. 
## The intercepts beta_ij replace the old means mu_ij; gamma is the covariate slope.
##
## Consequence for the causal contrast Delta S = S1 - S0 & ICA calculation:
##   mu*_ij(x) = (beta1_ij - beta0_ij) + (gamma1 - gamma0) * x
##   sigma*_ij = sigma00 + sigma11 - 2*rho_ij*sqrt(sigma00*sigma11)   (unchanged form)
##
##   * If gamma0 == gamma1 (common slope across arms), the covariate CANCELS
##     in Delta S. The ICA is then a single COVARIATE-ADJUSTED value: the
##     association after removing baseline-titer variation. The only effect
##     of the covariate is to shrink the residual sigma00, sigma11.
##   * If gamma0 != gamma1, mu*_ij depends on x, so the ICA is a FUNCTION of
##     the baseline titer. Evaluate it at representative x (argument x_eval,
##     on the centered scale) or loop over a grid.
##
##
## STARTING VALUES (Theta.S_0, Theta.S_1):
##   Diff.Sigma = FALSE : 6 values = c(beta_00, beta_01, beta_10, beta_11,
##                                     gamma, sigma)
##   Diff.Sigma = TRUE  : 9 values = c(beta_00, beta_01, beta_10, beta_11,
##                                     gamma, sigma_00, sigma_01,
##                                     sigma_10, sigma_11)
##
## ARGUMENTS:
##   Cov         : unquoted column name of the covariate in Dataset.
##   Center.Cov  : center the covariate at its sample mean (default TRUE).
##   x_eval      : covariate value (centered scale) at which strata means of
##                 Delta S are evaluated for the ICA (default 0 = sample mean).
##   N.MC        : Monte Carlo sample size for the integrals I_ij (default 10000).
## =====================================================================

ICA.BinCont.Cov <- function(Dataset, Surr, True, Treat, Cov,
                        G_pi_10 = c(0,1),
                        G_rho_01_00=c(-1,1),
                        G_rho_01_01=c(-1,1),
                        G_rho_01_10=c(-1,1),
                        G_rho_01_11=c(-1,1),
                        Theta.S_0, # 4 intercepts + slope + sigma(s)
                        Theta.S_1,
                        M=1000,        # number of sensitivity-analysis draws
                        Seed=123,      # integer for standalone use; NULL inside an outer parallel loop
                        N.MC=10000,    # Monte Carlo size for the integrals I_ij
                        Center.Cov=TRUE,
                        x_eval=0,
                        Monotonicity=FALSE,
                        Independence=FALSE,
                        HAA=FALSE,
                        Cond_ind=FALSE,
                        Plots=TRUE, Save.Plots="No",
                        Show.Details=FALSE){
# 1. PREPARATION --------------------------------------------------------------

  Diff.Sigma = FALSE # All 4 SD values in each group are the same

  ## EXTRACT VARIABLE 
  Surr  <- Dataset[,paste(substitute(Surr))]
  True  <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Cov   <- Dataset[,paste(substitute(Cov))]

  ## CHECK DATA FORMAT
  data_no_miss.original <- data.frame(na.exclude(cbind(Surr, True, Treat, Cov)))
  if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (max(na.exclude(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
  if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}

  ## CHECK LENGTH OF STARTING VALUES
  need <- if (Diff.Sigma) 9 else 6
  if (length(Theta.S_0)!=need | length(Theta.S_1)!=need){
    stop(paste0("\nWith a covariate, Theta.S_0 and Theta.S_1 must have ",
                need, " values (4 intercepts + slope",
                if (Diff.Sigma) " + 4 sigmas)." else " + 1 sigma)."))
  }

  ## SEED SETUP 
  # Scientific stream: set once here, then left to flow through the loop.
  if (!is.null(Seed)) set.seed(Seed)
  # Make sure a random seed already exist, if not create one
  if (!exists(".Random.seed", envir = .GlobalEnv)) stats::runif(1)
  # Store the current state of the random seed to use later
  .sci_state <- get(".Random.seed", envir = .GlobalEnv)

  # Fixed N(0,1) integration patterns (reproduce the original seeds 1..4),
  # drawn once and held constant for every iteration. The scientific stream
  # saved above is restored immediately after, so these draws do not perturb it.
  set.seed(1L); z_10 <- rnorm(N.MC)
  set.seed(2L); z_01 <- rnorm(N.MC)
  set.seed(3L); z_00 <- rnorm(N.MC)
  set.seed(4L); z_11 <- rnorm(N.MC)
  
  assign(".Random.seed", .sci_state, envir = .GlobalEnv)

  ## COVARIATE CENTERING CONSTANT
  cov_center <- if (Center.Cov) mean(data_no_miss.original$Cov) else 0

  ## BUILD FUNCTIONS TO CALCULATE MIXTURE DEVIANCE -----------------------------
  # The mixing weights pi_ij are passed in as arguments (p00..p11).
  if (Diff.Sigma==TRUE){
  # Unequal sigma's: theta = c(b00,b01,b10,b11, gamma, s00,s01,s10,s11)
    mixt.deviance <- function(theta, data, covar, p00, p01, p10, p11){
      gc <- theta[5]*covar          # gamma * covariate (shared by all 4 means)
      # Specify the density function of normal mixture models
      pdf <- p00*dnorm(data, theta[1]+gc, theta[6]) +
             p01*dnorm(data, theta[2]+gc, theta[7]) +
             p10*dnorm(data, theta[3]+gc, theta[8]) +
             p11*dnorm(data, theta[4]+gc, theta[9])
      # Deviance = -2 log-likelihood
      -2*sum(log(pdf))
    }
  } else {
    ## Equal sigma's: theta = c(b00,b01,b10,b11, gamma, sigma)
    mixt.deviance <- function(theta, data, covar, p00, p01, p10, p11){
      sigma_all <- theta[6]
      gc <- theta[5]*covar          # gamma * covariate (shared by all 4 means)
      pdf <- p00*dnorm(data, theta[1]+gc, sigma_all) +
             p01*dnorm(data, theta[2]+gc, sigma_all) +
             p10*dnorm(data, theta[3]+gc, sigma_all) +
             p11*dnorm(data, theta[4]+gc, sigma_all)
      -2*sum(log(pdf))
    }
  }

  ## PRE-ALLOCATE RESULT ACCUMULATORS -------------------------------------------
  R2_H_all <- R2_C_all <- numeric(M)
  pi_00_all <- pi_01_all <- pi_10_all <- pi_11_all <- numeric(M)
  G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- numeric(M)
  pi_Delta_T_min1_all <- pi_Delta_T_0_all <- pi_Delta_T_1_all <- numeric(M)
  pi_0_00_e_all <- pi_0_01_e_all <- pi_0_10_e_all <- pi_0_11_e_all <- numeric(M)
  mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- numeric(M)
  sigma_00_00_all <- sigma_00_01_all <- sigma_00_10_all <- sigma_00_11_all <- numeric(M)
  pi_1_00_e_all <- pi_1_01_e_all <- pi_1_10_e_all <- pi_1_11_e_all <- numeric(M)
  mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- numeric(M)
  sigma_11_00_all <- sigma_11_01_all <- sigma_11_10_all <- sigma_11_11_all <- numeric(M)
  gamma_0_all <- gamma_1_all <- numeric(M)                              # covariate slopes
  mean_Y_S0_all <- mean_Y_S1_all <- numeric(M)
  var_Y_S0_all <- var_Y_S1_all <- numeric(M)
  dev_S0_all <- dev_S1_all <- numeric(M)
  AIC_S0_all <- AIC_S1_all <- AIC_tot_all <- numeric(M)   # AIC per arm and combined
  code_nlm_0_all <- code_nlm_1_all <- numeric(M)
  aantal <- 0

  ## COMPUTE PER-ARM VECTORS AND MARGINAL PROBABILITIES ------------------------
  dnm <- data_no_miss.original
  data_conttreat <- dnm[dnm$Treat == -1, , drop=FALSE]   # control arm
  data_exptreat  <- dnm[dnm$Treat ==  1, , drop=FALSE]   # experimental arm
  S_0 <- data_conttreat$Surr
  S_1 <- data_exptreat$Surr
  X_0 <- data_conttreat$Cov - cov_center # centered covariate, control
  X_1 <- data_exptreat$Cov  - cov_center # centered covariate, experimental

  # Estimate the marginal probabilities of T
  pi_punt1 <- mean(data_exptreat$True)  # E(T|Z=1)
  pi_punt0 <- 1 - pi_punt1
  pi_1punt <- mean(data_conttreat$True) # E(T|Z=0)
  pi_0punt <- 1 - pi_1punt
  if (Monotonicity==TRUE){
    if (pi_punt0 > pi_0punt | pi_1punt > pi_punt1) {stop("\nThe monotonicity assumption is not compatible with the data")}
    G_pi_10 = c(0,0)
  }

  # 2. SENSITIVITY ANALYSIS ----------------------------------------------------

  
  while (aantal < M){

    ## 2.1. Sample the unidentifiable parameters -------------------------------
    ## (no per-iteration set.seed; the stream set up above advances naturally)
    if(HAA==FALSE){
      G_rho_01_00_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
      G_rho_01_01_hier <- runif(1, min=G_rho_01_01[1], max=G_rho_01_01[2])
      G_rho_01_10_hier <- runif(1, min=G_rho_01_10[1], max=G_rho_01_10[2])
      G_rho_01_11_hier <- runif(1, min=G_rho_01_11[1], max=G_rho_01_11[2])
    }
    if(HAA==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
    }
    if(Cond_ind==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- 0
    }
    if(Show.Details==TRUE){
      flush.console()
      cat("\n\nG_rho:", cbind(G_rho_01_00_hier, G_rho_01_01_hier, G_rho_01_10_hier, G_rho_01_11_hier), "\n")
    }

    # Sample pi_10 and derive the remaining pi_ij compatible with the margins
    if(Independence==FALSE){
      pi_10_hier <- runif(1, min=G_pi_10[1], max=G_pi_10[2])
      pi_00_hier <- pi_punt0 - pi_10_hier
      pi_01_hier <- pi_0punt - pi_00_hier
      pi_11_hier <- pi_1punt - pi_10_hier
    }
    if(Independence==TRUE){
      pi_10_hier <- pi_1punt*pi_punt0
      pi_00_hier <- pi_0punt*pi_punt0
      pi_01_hier <- pi_0punt*pi_punt1
      pi_11_hier <- pi_1punt*pi_punt1
    }

    if ((pi_00_hier >= 0 & pi_01_hier >= 0 & pi_10_hier >= 0 & pi_11_hier >= 0) &
        (pi_00_hier <= 1 & pi_01_hier <= 1 & pi_10_hier <= 1 & pi_11_hier <= 1)){

      if(Show.Details==TRUE){
        flush.console()
        cat("\n\nPi's used here:", cbind(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), "\n")
      }
      
      ## 2.2. Calculate pi_delta_T -----------------------------------------------
      # individual causal effect of the treatment on the true endpoint: delta_T = T1 - T0
      pi_Delta_T_min1 <- pi_10_hier # Delta_T = -1
      pi_Delta_T_0    <- pi_00_hier + pi_11_hier # Delta_T = 0
      pi_Delta_T_1    <- pi_01_hier # Delta_T = 1

      ## 2.3. Estimate Gaussian mixture models of Y_S| Y_T -----------------------
      
      # Fit model: maximize likelihood = minimizing deviance
      # Iterative process: start with an initial guess of parameter values, 
      # then iteratively updat them to minimize deviance, using Newton-Raphson type algorithm

      options(warn = -999)
      mix1 <- nlm(f = mixt.deviance, p = Theta.S_0, data = S_0, covar = X_0,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000 #the maximum number of iterations to be performed before the program is terminated
                  )
      mix2 <- nlm(f = mixt.deviance, p = Theta.S_1, data = S_1, covar = X_1,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000)
      options(warn = 1)

      if (exists("mix1")==TRUE & exists("mix2")==TRUE & is.null(mix1)==FALSE & is.null(mix2)==FALSE &
          mix1$code<=3 & mix2$code<=3){

        ## 2.4. Extract fitted parameters (intercepts, slope, residual variances) & deviance ----
        # Control arm (S0)
        mu_0_00 <- mix1$estimate[1]; mu_0_01 <- mix1$estimate[2]
        mu_0_10 <- mix1$estimate[3]; mu_0_11 <- mix1$estimate[4]
        gamma_0 <- mix1$estimate[5]
        if (Diff.Sigma==FALSE){ sigma_00_00 <- sigma_00_01 <- sigma_00_10 <- sigma_00_11 <- mix1$estimate[6]**2 }
        if (Diff.Sigma==TRUE){
          sigma_00_00 <- mix1$estimate[6]**2; sigma_00_01 <- mix1$estimate[7]**2
          sigma_00_10 <- mix1$estimate[8]**2; sigma_00_11 <- mix1$estimate[9]**2 }
        pi_0_00_e <- pi_00_hier; pi_0_01_e <- pi_01_hier; pi_0_10_e <- pi_10_hier; pi_0_11_e <- pi_11_hier
        code_nlm_0 <- mix1$code

        # Experimental arm (S1)
        mu_1_00 <- mix2$estimate[1]; mu_1_01 <- mix2$estimate[2]
        mu_1_10 <- mix2$estimate[3]; mu_1_11 <- mix2$estimate[4]
        gamma_1 <- mix2$estimate[5]
        if (Diff.Sigma==FALSE){ sigma_11_00 <- sigma_11_01 <- sigma_11_10 <- sigma_11_11 <- mix2$estimate[6]**2 }
        if (Diff.Sigma==TRUE){
          sigma_11_00 <- mix2$estimate[6]**2; sigma_11_01 <- mix2$estimate[7]**2
          sigma_11_10 <- mix2$estimate[8]**2; sigma_11_11 <- mix2$estimate[9]**2 }
        pi_1_00_e <- pi_00_hier; pi_1_01_e <- pi_01_hier; pi_1_10_e <- pi_10_hier; pi_1_11_e <- pi_11_hier
        code_nlm_1 <- mix2$code

        # Conditional deviance of each fitted mixture-of-regressions (for GOF check)
        dev_S0 <- mix1$minimum
        dev_S1 <- mix2$minimum
        
        # AIC = deviance + 2*k. k = number of parameters ESTIMATED per arm
        # (4 intercepts + slope + sigma(s) = length of Theta.S, i.e. `need`).
        # AIC_tot is the joint criterion for the two independent
        # arms (deviances add, parameters add).
        k_arm   <- need
        AIC_S0  <- dev_S0 + 2*k_arm
        AIC_S1  <- dev_S1 + 2*k_arm
        AIC_tot <- (dev_S0 + dev_S1) + 2*(2*k_arm)

        if (Plots == TRUE){
          ## Marginal fitted density, averaged over the observed covariate values,
          ## overlaid on the marginal histogram of S (O(n^2), so only when plotting).
          sd_0 <- sqrt(sigma_00_00); sd_1 <- sqrt(sigma_11_00)
          comp_dens <- function(grid, b, gamma, covar, sdv){
            sapply(grid, function(s) mean(dnorm(s, b + gamma*covar, sdv)))
          }
          x_S_0 <- sort(S_0)
          S0_M1 <- pi_00_hier * comp_dens(x_S_0, mu_0_00, gamma_0, X_0, sd_0)
          S0_M2 <- pi_01_hier * comp_dens(x_S_0, mu_0_01, gamma_0, X_0, sd_0)
          S0_M3 <- pi_10_hier * comp_dens(x_S_0, mu_0_10, gamma_0, X_0, sd_0)
          S0_M4 <- pi_11_hier * comp_dens(x_S_0, mu_0_11, gamma_0, X_0, sd_0)
          Sum_M_S_0 <- S0_M1 + S0_M2 + S0_M3 + S0_M4
          x_S_1 <- sort(S_1)
          S1_M1 <- pi_00_hier * comp_dens(x_S_1, mu_1_00, gamma_1, X_1, sd_1)
          S1_M2 <- pi_01_hier * comp_dens(x_S_1, mu_1_01, gamma_1, X_1, sd_1)
          S1_M3 <- pi_10_hier * comp_dens(x_S_1, mu_1_10, gamma_1, X_1, sd_1)
          S1_M4 <- pi_11_hier * comp_dens(x_S_1, mu_1_11, gamma_1, X_1, sd_1)
          Sum_M_S_1 <- S1_M1 + S1_M2 + S1_M3 + S1_M4

          plot_one <- function(){
            par(mfrow=c(1,2))
            max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
            hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]),
                 main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_0, Sum_M_S_0, col="red", lwd=4, lty=1)
            lines(x_S_0, S0_M1, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M2, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M3, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M4, col="blue", lwd=2, lty=3)
            max_val <- max(max(density(S_1)$y), max(Sum_M_S_1))
            hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]),
                 main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_1, Sum_M_S_1, col="red", lwd=4, lty=1)
            lines(x_S_1, S1_M1, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M2, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M3, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M4, col="blue", lwd=2, lty=3)
            par(mfrow=c(1,1))
          }
          if (Save.Plots!="No"){
            options(warn=-1); dir.create(Save.Plots); options(warn=0)
            pdf(paste(Save.Plots, "/Run", (aantal+1), ".pdf", sep=""),width=6,height=4,paper='special')
            plot_one(); dev.off()
          }
          plot_one()
        }

        ## 2.5 Covariance-matrix checks (PD and condition number) --------------
        # If checks do not pass, the results are discared (so we always have M valid results)
        mat_a <- matrix(c(sigma_00_00,
                         (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                         (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                          sigma_11_00), nrow = 2, byrow = TRUE)
        eigen_mat_a <- min(eigen((mat_a))$values)
        mat_b <- matrix(c(sigma_00_01,
                         (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                         (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                          sigma_11_01), nrow = 2, byrow = TRUE)
        eigen_mat_b <- min(eigen((mat_b))$values)
        mat_c <- matrix(c(sigma_00_10,
                         (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                         (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                          sigma_11_10), nrow = 2, byrow = TRUE)
        eigen_mat_c <- min(eigen((mat_c))$values)
        mat_d <- matrix(c(sigma_00_11,
                         (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                         (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                          sigma_11_11), nrow = 2, byrow = TRUE)
        eigen_mat_d <- min(eigen((mat_d))$values)

        CN_mat_a <- max(svd(mat_a)$d)/min(svd(mat_a)$d)
        CN_mat_b <- max(svd(mat_b)$d)/min(svd(mat_b)$d)
        CN_mat_c <- max(svd(mat_c)$d)/min(svd(mat_c)$d)
        CN_mat_d <- max(svd(mat_d)$d)/min(svd(mat_d)$d)

        if (eigen_mat_a > 0 & eigen_mat_b > 0 & eigen_mat_c > 0 & eigen_mat_d > 0 &
            CN_mat_a < 50 & CN_mat_b < 50 & CN_mat_c < 50 & CN_mat_d < 50){

          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nMin Eigenvalues:", cbind(eigen_mat_a, eigen_mat_b, eigen_mat_c, eigen_mat_d), "\n")
            cat("\n\nCondition number:", cbind(CN_mat_a, CN_mat_b, CN_mat_c, CN_mat_d), "\n")
            cat("\n\nChecks passed...")
          }

          ## 2.6 Marginal moments of S0, S1 (for GOF check) --------------------
          m0 <- mean(X_0); v0 <- var(X_0)   # control-arm covariate mean and variance
          m1 <- mean(X_1); v1 <- var(X_1)   # experimental-arm covariate mean and variance

          mean_int_S0 <- (pi_0_00_e*mu_0_00) + (pi_0_01_e*mu_0_01) + (pi_0_10_e*mu_0_10) + (pi_0_11_e*mu_0_11)
          mean_int_S1 <- (pi_1_00_e*mu_1_00) + (pi_1_01_e*mu_1_01) + (pi_1_10_e*mu_1_10) + (pi_1_11_e*mu_1_11)

          mean_Y_S0 <- mean_int_S0 + gamma_0 * m0   # E(S | Z=0)
          mean_Y_S1 <- mean_int_S1 + gamma_1 * m1   # E(S | Z=1)

          var_Y_S0 <- ((pi_0_00_e*mu_0_00^2)+(pi_0_01_e*mu_0_01^2)+(pi_0_10_e*mu_0_10^2)+(pi_0_11_e*mu_0_11^2)) - mean_int_S0^2 +
                      ((pi_0_00_e*sigma_00_00)+(pi_0_01_e*sigma_00_01)+(pi_0_10_e*sigma_00_10)+(pi_0_11_e*sigma_00_11)) +
                      (gamma_0^2 * v0)

          var_Y_S1 <- ((pi_1_00_e*mu_1_00^2)+(pi_1_01_e*mu_1_01^2)+(pi_1_10_e*mu_1_10^2)+(pi_1_11_e*mu_1_11^2)) - mean_int_S1^2 +
                      ((pi_1_00_e*sigma_11_00)+(pi_1_01_e*sigma_11_01)+(pi_1_10_e*sigma_11_10)+(pi_1_11_e*sigma_11_11)) +
                      (gamma_1^2 * v1)

          ## 2.7 Calculate ICA via Monte Carlo ---------------------------------
          ## Treatment effect on the surrogate within each T-stratum, at covariate x_eval.
          ## mu*_ij(x) = (beta1_ij - beta0_ij) + (gamma1 - gamma0) * x_eval
          xs <- x_eval
          mu_s00 <- (mu_1_00 - mu_0_00) + (gamma_1 - gamma_0)*xs
          mu_s01 <- (mu_1_01 - mu_0_01) + (gamma_1 - gamma_0)*xs
          mu_s10 <- (mu_1_10 - mu_0_10) + (gamma_1 - gamma_0)*xs
          mu_s11 <- (mu_1_11 - mu_0_11) + (gamma_1 - gamma_0)*xs
          
          # Weights of pi_00 and pi_11
          omega_1 <- pi_00_hier / (pi_00_hier + pi_11_hier); omega_2 <- 1 - omega_1

          sigma_s00 <- sigma_00_00 + sigma_11_00 - (2 * (sqrt(sigma_00_00 * sigma_11_00) * G_rho_01_00_hier))
          sigma_s01 <- sigma_00_01 + sigma_11_01 - (2 * (sqrt(sigma_00_01 * sigma_11_01) * G_rho_01_01_hier))
          sigma_s10 <- sigma_00_10 + sigma_11_10 - (2 * (sqrt(sigma_00_10 * sigma_11_10) * G_rho_01_10_hier))
          sigma_s11 <- sigma_00_11 + sigma_11_11 - (2 * (sqrt(sigma_00_11 * sigma_11_11) * G_rho_01_11_hier))

          f_Delta_S <- function(val){
            (pi_00_hier * dnorm(val, mu_s00, sd = sqrt(sigma_s00))) +
            (pi_01_hier * dnorm(val, mu_s01, sd = sqrt(sigma_s01))) +
            (pi_10_hier * dnorm(val, mu_s10, sd = sqrt(sigma_s10))) +
            (pi_11_hier * dnorm(val, mu_s11, sd = sqrt(sigma_s11)))
          }

          ## Sample S1, S2, S3, S4 from normal distribution
          ## S = mu + sd * z  is identical to rnorm(n, mu, sd) under the same seed.
          ## z values were sampled from N(0,1) on line 127-130
          S1mc <- mu_s10 + sqrt(sigma_s10) * z_10
          I_10 <- mean(log(dnorm(S1mc, mu_s10, sd = sqrt(sigma_s10)) / f_Delta_S(S1mc)))

          S2mc <- mu_s01 + sqrt(sigma_s01) * z_01
          I_01 <- mean(log(dnorm(S2mc, mu_s01, sd = sqrt(sigma_s01)) / f_Delta_S(S2mc)))

          S3mc <- mu_s00 + sqrt(sigma_s00) * z_00
          I_00 <- mean(log(((omega_1 * dnorm(S3mc, mu_s00, sd = sqrt(sigma_s00))) +
                            (omega_2 * dnorm(S3mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S3mc)))

          S4mc <- mu_s11 + sqrt(sigma_s11) * z_11
          I_11 <- mean(log(((omega_1 * dnorm(S4mc, mu_s00, sd = sqrt(sigma_s00))) +
                            (omega_2 * dnorm(S4mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S4mc)))
          
          # Calculate ICA
          I_Delta_T_Delta_S <- (pi_00_hier * I_00) + (pi_01_hier * I_01) + (pi_10_hier * I_10) + (pi_11_hier * I_11)
          H_Delta_T <- -((pi_Delta_T_min1 * log(pi_Delta_T_min1)) +
                         (pi_Delta_T_0 * log(pi_Delta_T_0)) +
                         (pi_Delta_T_1 * log(pi_Delta_T_1)))
          if(Monotonicity==TRUE){
            H_Delta_T <- -((pi_Delta_T_0 * log(pi_Delta_T_0)) +
                           (pi_Delta_T_1 * log(pi_Delta_T_1)))
          }
          R2_H <- I_Delta_T_Delta_S / H_Delta_T
          R2_C <- 1 - exp(-2*I_Delta_T_Delta_S)

          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nR2H:", R2_H, "  R2C:", R2_C, "\n")
          }

          aantal <- aantal + 1
          if(Show.Details==TRUE){
            cat("\n", (aantal/M)*100, "% done. \n", sep="")
          }

          # save results (indexed assignment into pre-allocated vectors)
          R2_H_all[aantal] <- R2_H
          R2_C_all[aantal] <- R2_C
          G_rho_01_00_all[aantal] <- G_rho_01_00_hier
          G_rho_01_01_all[aantal] <- G_rho_01_01_hier
          G_rho_01_10_all[aantal] <- G_rho_01_10_hier
          G_rho_01_11_all[aantal] <- G_rho_01_11_hier
          pi_00_all[aantal] <- pi_00_hier
          pi_01_all[aantal] <- pi_01_hier
          pi_10_all[aantal] <- pi_10_hier
          pi_11_all[aantal] <- pi_11_hier
          pi_Delta_T_min1_all[aantal] <- pi_Delta_T_min1
          pi_Delta_T_0_all[aantal] <- pi_Delta_T_0
          pi_Delta_T_1_all[aantal] <- pi_Delta_T_1
          pi_0_00_e_all[aantal] <- pi_0_00_e
          pi_0_01_e_all[aantal] <- pi_0_01_e
          pi_0_10_e_all[aantal] <- pi_0_10_e
          pi_0_11_e_all[aantal] <- pi_0_11_e
          mu_0_00_all[aantal] <- mu_0_00
          mu_0_01_all[aantal] <- mu_0_01
          mu_0_10_all[aantal] <- mu_0_10
          mu_0_11_all[aantal] <- mu_0_11
          sigma_00_00_all[aantal] <- sigma_00_00
          sigma_00_01_all[aantal] <- sigma_00_01
          sigma_00_10_all[aantal] <- sigma_00_10
          sigma_00_11_all[aantal] <- sigma_00_11
          pi_1_00_e_all[aantal] <- pi_1_00_e
          pi_1_01_e_all[aantal] <- pi_1_01_e
          pi_1_10_e_all[aantal] <- pi_1_10_e
          pi_1_11_e_all[aantal] <- pi_1_11_e
          mu_1_00_all[aantal] <- mu_1_00
          mu_1_01_all[aantal] <- mu_1_01
          mu_1_10_all[aantal] <- mu_1_10
          mu_1_11_all[aantal] <- mu_1_11
          sigma_11_00_all[aantal] <- sigma_11_00
          sigma_11_01_all[aantal] <- sigma_11_01
          sigma_11_10_all[aantal] <- sigma_11_10
          sigma_11_11_all[aantal] <- sigma_11_11
          gamma_0_all[aantal] <- gamma_0
          gamma_1_all[aantal] <- gamma_1
          mean_Y_S0_all[aantal] <- mean_Y_S0
          mean_Y_S1_all[aantal] <- mean_Y_S1
          var_Y_S0_all[aantal] <- var_Y_S0
          var_Y_S1_all[aantal] <- var_Y_S1
          dev_S0_all[aantal] <- dev_S0
          dev_S1_all[aantal] <- dev_S1
          AIC_S0_all[aantal] <- AIC_S0
          AIC_S1_all[aantal] <- AIC_S1
          AIC_tot_all[aantal] <- AIC_tot
          code_nlm_0_all[aantal] <- code_nlm_0
          code_nlm_1_all[aantal] <- code_nlm_1

        } # min eigen and CN requirements
      } # mixture distribution available
    } # pi's requirements
  } # end while loop

  # 3. OUTPUT RESULTS ----------------------------------------------------
  fit <- list(R2_H=R2_H_all, R2_C=R2_C_all,
              pi_00=as.numeric(pi_00_all), pi_01=as.numeric(pi_01_all),
              pi_10=as.numeric(pi_10_all), pi_11=as.numeric(pi_11_all),
              G_rho_01_00=as.numeric(G_rho_01_00_all),
              G_rho_01_01=as.numeric(G_rho_01_01_all),
              G_rho_01_10=as.numeric(G_rho_01_10_all),
              G_rho_01_11=as.numeric(G_rho_01_11_all),
              pi_Delta_T_min1=as.numeric(pi_Delta_T_min1_all),
              pi_Delta_T_0=as.numeric(pi_Delta_T_0_all),
              pi_Delta_T_1=as.numeric(pi_Delta_T_1_all),
              pi_0_00=as.numeric(pi_0_00_e_all), pi_0_01=as.numeric(pi_0_01_e_all),
              pi_0_10=as.numeric(pi_0_10_e_all), pi_0_11=as.numeric(pi_0_11_e_all),
              mu_0_00=as.numeric(mu_0_00_all), mu_0_01=as.numeric(mu_0_01_all),
              mu_0_10=as.numeric(mu_0_10_all), mu_0_11=as.numeric(mu_0_11_all),
              sigma2_00_00=as.numeric(sigma_00_00_all), sigma2_00_01=as.numeric(sigma_00_01_all),
              sigma2_00_10=as.numeric(sigma_00_10_all), sigma2_00_11=as.numeric(sigma_00_11_all),
              pi_1_00=as.numeric(pi_1_00_e_all), pi_1_01=as.numeric(pi_1_01_e_all),
              pi_1_10=as.numeric(pi_1_10_e_all), pi_1_11=as.numeric(pi_1_11_e_all),
              mu_1_00=as.numeric(mu_1_00_all), mu_1_01=as.numeric(mu_1_01_all),
              mu_1_10=as.numeric(mu_1_10_all), mu_1_11=as.numeric(mu_1_11_all),
              sigma2_11_00=as.numeric(sigma_11_00_all), sigma2_11_01=as.numeric(sigma_11_01_all),
              sigma2_11_10=as.numeric(sigma_11_10_all), sigma2_11_11=as.numeric(sigma_11_11_all),
              gamma_0=as.numeric(gamma_0_all),         # covariate slope, control arm
              gamma_1=as.numeric(gamma_1_all),         # covariate slope, experimental arm
              cov_center=cov_center,                   # centering constant
              x_eval=x_eval,                           # covariate value used for the ICA
              mean_Y_S0=as.numeric(mean_Y_S0_all), mean_Y_S1=as.numeric(mean_Y_S1_all),
              var_Y_S0=as.numeric(var_Y_S0_all), var_Y_S1=as.numeric(var_Y_S1_all),
              dev_S0=as.numeric(dev_S0_all), dev_S1=as.numeric(dev_S1_all),
              AIC_S0=as.numeric(AIC_S0_all), AIC_S1=as.numeric(AIC_S1_all),
              AIC_tot=as.numeric(AIC_tot_all),       # joint AIC of the two arms
              npar=need,                             # estimated parameters PER arm (k)
              code_nlm_0=as.numeric(code_nlm_0_all), code_nlm_1=as.numeric(code_nlm_1_all),
              mean.S0=mean(S_0, na.rm=TRUE), var.S0=var(S_0, na.rm=TRUE),
              mean.S1=mean(S_1, na.rm=TRUE), var.S1=var(S_1, na.rm=TRUE),
              Call=match.call())

  class(fit) <- "ICA.BinCont"
  fit
} # end of function




## =====================================================================
## ICA.BinCont with TWO covariates in the component means mu_ij
## ---------------------------------------------------------------------
##
## Model (homoscedastic case, Diff.Sigma = FALSE):
##   Control arm:       S0_k ~ sum_ij pi_ij * N( beta0_ij + gamma0_1 * X1_k
##                                                        + gamma0_2 * X2_k , sigma00 )
##   Experimental arm:  S1_k ~ sum_ij pi_ij * N( beta1_ij + gamma1_1 * X1_k
##                                                        + gamma1_2 * X2_k , sigma11 )
##   where X1_k, X2_k are two covariates, each centered at its sample mean.
##
## This is a finite mixture of LINEAR REGRESSIONS (now with two covariate
## slopes per arm) and KNOWN mixing weights pi_ij (fixed in step 3 of the
## sensitivity analysis). The intercepts beta_ij replace the old means mu_ij;
## gamma_1, gamma_2 are the two covariate slopes (shared across the four
## mixture components within an arm, as in the one-covariate version).
## =====================================================================

ICA.BinCont.Cov2 <- function(Dataset, Surr, True, Treat, Cov1, Cov2,
                             G_pi_10 = c(0,1),
                             G_rho_01_00=c(-1,1),
                             G_rho_01_01=c(-1,1),
                             G_rho_01_10=c(-1,1),
                             G_rho_01_11=c(-1,1),
                             Theta.S_0, # 4 intercepts + 2 slopes + sigma(s)
                             Theta.S_1,
                             M=1000,        # number of sensitivity-analysis draws
                             Seed=123,      # integer for standalone use; NULL inside an outer parallel loop
                             N.MC=10000,    # Monte Carlo size for the integrals I_ij
                             Center.Cov=TRUE,
                             x_eval=c(0,0), # length-2: (x1, x2) on centered scale
                             Monotonicity=FALSE,
                             Independence=FALSE,
                             HAA=FALSE,
                             Cond_ind=FALSE,
                             Plots=TRUE, Save.Plots="No",
                             Show.Details=FALSE){
  # 1. PREPARATION --------------------------------------------------------------

  Diff.Sigma = FALSE # All 4 SD values in each group are the same
  
  ## EXTRACT VARIABLE
  Surr  <- Dataset[,paste(substitute(Surr))]
  True  <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Cov1  <- Dataset[,paste(substitute(Cov1))]
  Cov2  <- Dataset[,paste(substitute(Cov2))]
  
  ## CHECK DATA FORMAT 
  data_no_miss.original <- data.frame(na.exclude(cbind(Surr, True, Treat, Cov1, Cov2)))
  if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (max(na.exclude(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
  if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  
  ## CHECK x_eval LENGTH 
  if (length(x_eval)!=2) stop("\nWith two covariates, x_eval must be a length-2 vector c(x1, x2).")
  
  ## CHECK LENGTH OF STARTING VALUES 
  need <- if (Diff.Sigma) 10 else 7
  if (length(Theta.S_0)!=need | length(Theta.S_1)!=need){
    stop(paste0("\nWith two covariates, Theta.S_0 and Theta.S_1 must have ",
                need, " values (4 intercepts + 2 slopes",
                if (Diff.Sigma) " + 4 sigmas)." else " + 1 sigma)."))
  }
  
  ## SEED SETUP 
  # Scientific stream: set once here, then left to flow through the loop.
  if (!is.null(Seed)) set.seed(Seed)
  # Make sure a random seed already exist, if not create one
  if (!exists(".Random.seed", envir = .GlobalEnv)) stats::runif(1)
  # Store the current state of the random seed to use later
  .sci_state <- get(".Random.seed", envir = .GlobalEnv)
  
  # Fixed N(0,1) integration patterns (reproduce the original seeds 1..4),
  # drawn once and held constant for every iteration. The scientific stream
  # saved above is restored immediately after, so these draws do not perturb it.
  set.seed(1L); z_10 <- rnorm(N.MC)
  set.seed(2L); z_01 <- rnorm(N.MC)
  set.seed(3L); z_00 <- rnorm(N.MC)
  set.seed(4L); z_11 <- rnorm(N.MC)
  
  assign(".Random.seed", .sci_state, envir = .GlobalEnv)
  
  ## COVARIATE CENTERING CONSTANTS
  cov1_center <- if (Center.Cov) mean(data_no_miss.original$Cov1) else 0
  cov2_center <- if (Center.Cov) mean(data_no_miss.original$Cov2) else 0
  
  ## BUILD FUNCTIONS TO CALCULATE MIXTURE DEVIANCE
  # The mixing weights pi_ij are passed in as arguments (p00..p11).
  # covar1, covar2 are the two centered covariate vectors.
  if (Diff.Sigma==TRUE){
    # Unequal sigma's: theta = c(b00,b01,b10,b11, gamma1, gamma2, s00,s01,s10,s11)
    mixt.deviance <- function(theta, data, covar1, covar2, p00, p01, p10, p11){
      gc <- theta[5]*covar1 + theta[6]*covar2   # gamma1*X1 + gamma2*X2 (shared by all 4 means)
      # Specify the density function of normal mixture models
      pdf <- p00*dnorm(data, theta[1]+gc, theta[7]) +
        p01*dnorm(data, theta[2]+gc, theta[8]) +
        p10*dnorm(data, theta[3]+gc, theta[9]) +
        p11*dnorm(data, theta[4]+gc, theta[10])
      # Deviance = -2 log-likelihood
      -2*sum(log(pdf))
    }
  } else {
    ## Equal sigma's: theta = c(b00,b01,b10,b11, gamma1, gamma2, sigma)
    mixt.deviance <- function(theta, data, covar1, covar2, p00, p01, p10, p11){
      sigma_all <- theta[7]
      gc <- theta[5]*covar1 + theta[6]*covar2   # gamma1*X1 + gamma2*X2 (shared by all 4 means)
      pdf <- p00*dnorm(data, theta[1]+gc, sigma_all) +
        p01*dnorm(data, theta[2]+gc, sigma_all) +
        p10*dnorm(data, theta[3]+gc, sigma_all) +
        p11*dnorm(data, theta[4]+gc, sigma_all)
      -2*sum(log(pdf))
    }
  }
  
  ## PRE-ALLOCATE RESULT ACCUMULATORS -------------------------------------------
  R2_H_all <- R2_C_all <- numeric(M)
  pi_00_all <- pi_01_all <- pi_10_all <- pi_11_all <- numeric(M)
  G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- numeric(M)
  pi_Delta_T_min1_all <- pi_Delta_T_0_all <- pi_Delta_T_1_all <- numeric(M)
  pi_0_00_e_all <- pi_0_01_e_all <- pi_0_10_e_all <- pi_0_11_e_all <- numeric(M)
  mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- numeric(M)
  sigma_00_00_all <- sigma_00_01_all <- sigma_00_10_all <- sigma_00_11_all <- numeric(M)
  pi_1_00_e_all <- pi_1_01_e_all <- pi_1_10_e_all <- pi_1_11_e_all <- numeric(M)
  mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- numeric(M)
  sigma_11_00_all <- sigma_11_01_all <- sigma_11_10_all <- sigma_11_11_all <- numeric(M)
  gamma_0_1_all <- gamma_0_2_all <- numeric(M)                          # covariate slopes, control
  gamma_1_1_all <- gamma_1_2_all <- numeric(M)                          # covariate slopes, experimental
  mean_Y_S0_all <- mean_Y_S1_all <- numeric(M)
  var_Y_S0_all <- var_Y_S1_all <- numeric(M)
  dev_S0_all <- dev_S1_all <- numeric(M)
  AIC_S0_all <- AIC_S1_all <- AIC_tot_all <- numeric(M)   # AIC per arm and combined
  code_nlm_0_all <- code_nlm_1_all <- numeric(M)
  aantal <- 0
  
  ## COMPUTE PER-ARM VECTORS AND MARGINAL PROBABILITIES ------------------------
  dnm <- data_no_miss.original
  data_conttreat <- dnm[dnm$Treat == -1, , drop=FALSE]   # control arm
  data_exptreat  <- dnm[dnm$Treat ==  1, , drop=FALSE]   # experimental arm
  S_0 <- data_conttreat$Surr
  S_1 <- data_exptreat$Surr
  X1_0 <- data_conttreat$Cov1 - cov1_center # centered covariate 1, control
  X1_1 <- data_exptreat$Cov1  - cov1_center # centered covariate 1, experimental
  X2_0 <- data_conttreat$Cov2 - cov2_center # centered covariate 2, control
  X2_1 <- data_exptreat$Cov2  - cov2_center # centered covariate 2, experimental
  
  # Estimate the marginal probabilities of T
  pi_punt1 <- mean(data_exptreat$True)  # E(T|Z=1)
  pi_punt0 <- 1 - pi_punt1
  pi_1punt <- mean(data_conttreat$True) # E(T|Z=0)
  pi_0punt <- 1 - pi_1punt
  if (Monotonicity==TRUE){
    if (pi_punt0 > pi_0punt | pi_1punt > pi_punt1) {stop("\nThe monotonicity assumption is not compatible with the data")}
    G_pi_10 = c(0,0)
  }
  
  #  ---------------------------------------------------------------------------
  # 2. SENSITIVITY ANALYSIS ----------------------------------------------------
  # ----------------------------------------------------------------------------
  
  
  while (aantal < M){
    
    ## 2.1. Sample the unidentifiable parameters -------------------------------
    ## (no per-iteration set.seed; the stream set up above advances naturally)
    if(HAA==FALSE){
      G_rho_01_00_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
      G_rho_01_01_hier <- runif(1, min=G_rho_01_01[1], max=G_rho_01_01[2])
      G_rho_01_10_hier <- runif(1, min=G_rho_01_10[1], max=G_rho_01_10[2])
      G_rho_01_11_hier <- runif(1, min=G_rho_01_11[1], max=G_rho_01_11[2])
    }
    if(HAA==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
    }
    if(Cond_ind==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- 0
    }
    if(Show.Details==TRUE){
      flush.console()
      cat("\n\nG_rho:", cbind(G_rho_01_00_hier, G_rho_01_01_hier, G_rho_01_10_hier, G_rho_01_11_hier), "\n")
    }
    
    # Sample pi_10 and derive the remaining pi_ij compatible with the margins
    if(Independence==FALSE){
      pi_10_hier <- runif(1, min=G_pi_10[1], max=G_pi_10[2])
      pi_00_hier <- pi_punt0 - pi_10_hier
      pi_01_hier <- pi_0punt - pi_00_hier
      pi_11_hier <- pi_1punt - pi_10_hier
    }
    if(Independence==TRUE){
      pi_10_hier <- pi_1punt*pi_punt0
      pi_00_hier <- pi_0punt*pi_punt0
      pi_01_hier <- pi_0punt*pi_punt1
      pi_11_hier <- pi_1punt*pi_punt1
    }
    
    if ((pi_00_hier >= 0 & pi_01_hier >= 0 & pi_10_hier >= 0 & pi_11_hier >= 0) &
        (pi_00_hier <= 1 & pi_01_hier <= 1 & pi_10_hier <= 1 & pi_11_hier <= 1)){
      
      if(Show.Details==TRUE){
        flush.console()
        cat("\n\nPi's used here:", cbind(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), "\n")
      }
      
      ## 2.2. Calculate pi_delta_T -----------------------------------------------
      # individual causal effect of the treatment on the true endpoint: delta_T = T1 - T0
      pi_Delta_T_min1 <- pi_10_hier # Delta_T = -1
      pi_Delta_T_0    <- pi_00_hier + pi_11_hier # Delta_T = 0
      pi_Delta_T_1    <- pi_01_hier # Delta_T = 1
      
      ## 2.3. Estimate Gaussian mixture models of Y_S| Y_T -----------------------
      
      # Fit model: maximize likelihood = minimizing deviance
      # Iterative process: start with an initial guess of parameter values,
      # then iteratively updat them to minimize deviance, using Newton-Raphson type algorithm
      
      options(warn = -999)
      mix1 <- nlm(f = mixt.deviance, p = Theta.S_0, data = S_0, covar1 = X1_0, covar2 = X2_0,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000 #the maximum number of iterations to be performed before the program is terminated
      )
      mix2 <- nlm(f = mixt.deviance, p = Theta.S_1, data = S_1, covar1 = X1_1, covar2 = X2_1,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000)
      options(warn = 1)
      
      if (exists("mix1")==TRUE & exists("mix2")==TRUE & is.null(mix1)==FALSE & is.null(mix2)==FALSE &
          mix1$code<=3 & mix2$code<=3){
        
        ## 2.4. Extract fitted parameters (intercepts, slopes, residual variances) & deviance ----
        # Control arm (S0)
        mu_0_00 <- mix1$estimate[1]; mu_0_01 <- mix1$estimate[2]
        mu_0_10 <- mix1$estimate[3]; mu_0_11 <- mix1$estimate[4]
        gamma_0_1 <- mix1$estimate[5]; gamma_0_2 <- mix1$estimate[6]
        if (Diff.Sigma==FALSE){ sigma_00_00 <- sigma_00_01 <- sigma_00_10 <- sigma_00_11 <- mix1$estimate[7]**2 }
        if (Diff.Sigma==TRUE){
          sigma_00_00 <- mix1$estimate[7]**2; sigma_00_01 <- mix1$estimate[8]**2
          sigma_00_10 <- mix1$estimate[9]**2; sigma_00_11 <- mix1$estimate[10]**2 }
        pi_0_00_e <- pi_00_hier; pi_0_01_e <- pi_01_hier; pi_0_10_e <- pi_10_hier; pi_0_11_e <- pi_11_hier
        code_nlm_0 <- mix1$code
        
        # Experimental arm (S1)
        mu_1_00 <- mix2$estimate[1]; mu_1_01 <- mix2$estimate[2]
        mu_1_10 <- mix2$estimate[3]; mu_1_11 <- mix2$estimate[4]
        gamma_1_1 <- mix2$estimate[5]; gamma_1_2 <- mix2$estimate[6]
        if (Diff.Sigma==FALSE){ sigma_11_00 <- sigma_11_01 <- sigma_11_10 <- sigma_11_11 <- mix2$estimate[7]**2 }
        if (Diff.Sigma==TRUE){
          sigma_11_00 <- mix2$estimate[7]**2; sigma_11_01 <- mix2$estimate[8]**2
          sigma_11_10 <- mix2$estimate[9]**2; sigma_11_11 <- mix2$estimate[10]**2 }
        pi_1_00_e <- pi_00_hier; pi_1_01_e <- pi_01_hier; pi_1_10_e <- pi_10_hier; pi_1_11_e <- pi_11_hier
        code_nlm_1 <- mix2$code
        
        # Conditional deviance of each fitted mixture-of-regressions (for GOF check)
        dev_S0 <- mix1$minimum
        dev_S1 <- mix2$minimum
        
        # AIC = deviance + 2*k. k = number of parameters ESTIMATED per arm
        # (4 intercepts + slope + sigma(s) = length of Theta.S, i.e. `need`).
        # The mixing weights pi_ij are fixed inputs of the sensitivity draw,
        # not estimated by MLE, so they are NOT counted in k. Each arm is a
        # separate fit; AIC_tot is the joint criterion for the two independent
        # arms (deviances add, parameters add).
        k_arm   <- need
        AIC_S0  <- dev_S0 + 2*k_arm
        AIC_S1  <- dev_S1 + 2*k_arm
        AIC_tot <- (dev_S0 + dev_S1) + 2*(2*k_arm)
        
        if (Plots == TRUE){
          ## Marginal fitted density, averaged over the observed covariate values,
          ## overlaid on the marginal histogram of S (O(n^2), so only when plotting).
          sd_0 <- sqrt(sigma_00_00); sd_1 <- sqrt(sigma_11_00)
          comp_dens <- function(grid, b, g1, g2, covar1, covar2, sdv){
            sapply(grid, function(s) mean(dnorm(s, b + g1*covar1 + g2*covar2, sdv)))
          }
          x_S_0 <- sort(S_0)
          S0_M1 <- pi_00_hier * comp_dens(x_S_0, mu_0_00, gamma_0_1, gamma_0_2, X1_0, X2_0, sd_0)
          S0_M2 <- pi_01_hier * comp_dens(x_S_0, mu_0_01, gamma_0_1, gamma_0_2, X1_0, X2_0, sd_0)
          S0_M3 <- pi_10_hier * comp_dens(x_S_0, mu_0_10, gamma_0_1, gamma_0_2, X1_0, X2_0, sd_0)
          S0_M4 <- pi_11_hier * comp_dens(x_S_0, mu_0_11, gamma_0_1, gamma_0_2, X1_0, X2_0, sd_0)
          Sum_M_S_0 <- S0_M1 + S0_M2 + S0_M3 + S0_M4
          x_S_1 <- sort(S_1)
          S1_M1 <- pi_00_hier * comp_dens(x_S_1, mu_1_00, gamma_1_1, gamma_1_2, X1_1, X2_1, sd_1)
          S1_M2 <- pi_01_hier * comp_dens(x_S_1, mu_1_01, gamma_1_1, gamma_1_2, X1_1, X2_1, sd_1)
          S1_M3 <- pi_10_hier * comp_dens(x_S_1, mu_1_10, gamma_1_1, gamma_1_2, X1_1, X2_1, sd_1)
          S1_M4 <- pi_11_hier * comp_dens(x_S_1, mu_1_11, gamma_1_1, gamma_1_2, X1_1, X2_1, sd_1)
          Sum_M_S_1 <- S1_M1 + S1_M2 + S1_M3 + S1_M4
          
          plot_one <- function(){
            par(mfrow=c(1,2))
            max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
            hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]),
                 main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_0, Sum_M_S_0, col="red", lwd=4, lty=1)
            lines(x_S_0, S0_M1, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M2, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M3, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M4, col="blue", lwd=2, lty=3)
            max_val <- max(max(density(S_1)$y), max(Sum_M_S_1))
            hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]),
                 main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_1, Sum_M_S_1, col="red", lwd=4, lty=1)
            lines(x_S_1, S1_M1, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M2, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M3, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M4, col="blue", lwd=2, lty=3)
            par(mfrow=c(1,1))
          }
          if (Save.Plots!="No"){
            options(warn=-1); dir.create(Save.Plots); options(warn=0)
            pdf(paste(Save.Plots, "/Run", (aantal+1), ".pdf", sep=""),width=6,height=4,paper='special')
            plot_one(); dev.off()
          }
          plot_one()
        }
        
        ## 2.5 Covariance-matrix checks (PD and condition number) --------------
        # If checks do not pass, the results are discared (so we always have M valid results)
        mat_a <- matrix(c(sigma_00_00,
                          (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                          (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                          sigma_11_00), nrow = 2, byrow = TRUE)
        eigen_mat_a <- min(eigen((mat_a))$values)
        mat_b <- matrix(c(sigma_00_01,
                          (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                          (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                          sigma_11_01), nrow = 2, byrow = TRUE)
        eigen_mat_b <- min(eigen((mat_b))$values)
        mat_c <- matrix(c(sigma_00_10,
                          (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                          (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                          sigma_11_10), nrow = 2, byrow = TRUE)
        eigen_mat_c <- min(eigen((mat_c))$values)
        mat_d <- matrix(c(sigma_00_11,
                          (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                          (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                          sigma_11_11), nrow = 2, byrow = TRUE)
        eigen_mat_d <- min(eigen((mat_d))$values)
        
        CN_mat_a <- max(svd(mat_a)$d)/min(svd(mat_a)$d)
        CN_mat_b <- max(svd(mat_b)$d)/min(svd(mat_b)$d)
        CN_mat_c <- max(svd(mat_c)$d)/min(svd(mat_c)$d)
        CN_mat_d <- max(svd(mat_d)$d)/min(svd(mat_d)$d)
        
        if (eigen_mat_a > 0 & eigen_mat_b > 0 & eigen_mat_c > 0 & eigen_mat_d > 0 &
            CN_mat_a < 50 & CN_mat_b < 50 & CN_mat_c < 50 & CN_mat_d < 50){
          
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nMin Eigenvalues:", cbind(eigen_mat_a, eigen_mat_b, eigen_mat_c, eigen_mat_d), "\n")
            cat("\n\nCondition number:", cbind(CN_mat_a, CN_mat_b, CN_mat_c, CN_mat_d), "\n")
            cat("\n\nChecks passed...")
          }
          
          ## 2.6 Marginal moments of S0, S1 (for GOF check) --------------------
          ## Two covariates: means, variances and covariance per arm.
          m0_1 <- mean(X1_0); v0_1 <- var(X1_0)   # control-arm covariate-1 mean and variance
          m0_2 <- mean(X2_0); v0_2 <- var(X2_0)   # control-arm covariate-2 mean and variance
          c0_12 <- cov(X1_0, X2_0)                # control-arm covariance(X1, X2)
          m1_1 <- mean(X1_1); v1_1 <- var(X1_1)   # experimental-arm covariate-1 mean and variance
          m1_2 <- mean(X2_1); v1_2 <- var(X2_1)   # experimental-arm covariate-2 mean and variance
          c1_12 <- cov(X1_1, X2_1)                # experimental-arm covariance(X1, X2)
          
          mean_int_S0 <- (pi_0_00_e*mu_0_00) + (pi_0_01_e*mu_0_01) + (pi_0_10_e*mu_0_10) + (pi_0_11_e*mu_0_11)
          mean_int_S1 <- (pi_1_00_e*mu_1_00) + (pi_1_01_e*mu_1_01) + (pi_1_10_e*mu_1_10) + (pi_1_11_e*mu_1_11)
          
          mean_Y_S0 <- mean_int_S0 + gamma_0_1 * m0_1 + gamma_0_2 * m0_2   # E(S | Z=0)
          mean_Y_S1 <- mean_int_S1 + gamma_1_1 * m1_1 + gamma_1_2 * m1_2   # E(S | Z=1)
          
          var_Y_S0 <- ((pi_0_00_e*mu_0_00^2)+(pi_0_01_e*mu_0_01^2)+(pi_0_10_e*mu_0_10^2)+(pi_0_11_e*mu_0_11^2)) - mean_int_S0^2 +
            ((pi_0_00_e*sigma_00_00)+(pi_0_01_e*sigma_00_01)+(pi_0_10_e*sigma_00_10)+(pi_0_11_e*sigma_00_11)) +
            (gamma_0_1^2 * v0_1) + (gamma_0_2^2 * v0_2) + (2 * gamma_0_1 * gamma_0_2 * c0_12)
          
          var_Y_S1 <- ((pi_1_00_e*mu_1_00^2)+(pi_1_01_e*mu_1_01^2)+(pi_1_10_e*mu_1_10^2)+(pi_1_11_e*mu_1_11^2)) - mean_int_S1^2 +
            ((pi_1_00_e*sigma_11_00)+(pi_1_01_e*sigma_11_01)+(pi_1_10_e*sigma_11_10)+(pi_1_11_e*sigma_11_11)) +
            (gamma_1_1^2 * v1_1) + (gamma_1_2^2 * v1_2) + (2 * gamma_1_1 * gamma_1_2 * c1_12)
          
          ## 2.7 Calculate ICA via Monte Carlo ---------------------------------
          ## Treatment effect on the surrogate within each T-stratum, at covariate x_eval=(x1,x2).
          ## mu*_ij(x) = (beta1_ij - beta0_ij) + (gamma1_1 - gamma0_1)*x1 + (gamma1_2 - gamma0_2)*x2
          x1s <- x_eval[1]; x2s <- x_eval[2]
          dgamma1 <- gamma_1_1 - gamma_0_1
          dgamma2 <- gamma_1_2 - gamma_0_2
          mu_s00 <- (mu_1_00 - mu_0_00) + dgamma1*x1s + dgamma2*x2s
          mu_s01 <- (mu_1_01 - mu_0_01) + dgamma1*x1s + dgamma2*x2s
          mu_s10 <- (mu_1_10 - mu_0_10) + dgamma1*x1s + dgamma2*x2s
          mu_s11 <- (mu_1_11 - mu_0_11) + dgamma1*x1s + dgamma2*x2s
          
          # Weights of pi_00 and pi_11
          omega_1 <- pi_00_hier / (pi_00_hier + pi_11_hier); omega_2 <- 1 - omega_1
          
          sigma_s00 <- sigma_00_00 + sigma_11_00 - (2 * (sqrt(sigma_00_00 * sigma_11_00) * G_rho_01_00_hier))
          sigma_s01 <- sigma_00_01 + sigma_11_01 - (2 * (sqrt(sigma_00_01 * sigma_11_01) * G_rho_01_01_hier))
          sigma_s10 <- sigma_00_10 + sigma_11_10 - (2 * (sqrt(sigma_00_10 * sigma_11_10) * G_rho_01_10_hier))
          sigma_s11 <- sigma_00_11 + sigma_11_11 - (2 * (sqrt(sigma_00_11 * sigma_11_11) * G_rho_01_11_hier))
          
          f_Delta_S <- function(val){
            (pi_00_hier * dnorm(val, mu_s00, sd = sqrt(sigma_s00))) +
              (pi_01_hier * dnorm(val, mu_s01, sd = sqrt(sigma_s01))) +
              (pi_10_hier * dnorm(val, mu_s10, sd = sqrt(sigma_s10))) +
              (pi_11_hier * dnorm(val, mu_s11, sd = sqrt(sigma_s11)))
          }
          
          ## Sample S1, S2, S3, S4 from normal distribution
          ## S = mu + sd * z  is identical to rnorm(n, mu, sd) under the same seed.
          ## z values were sampled from N(0,1) in the SEED SETUP block above.
          S1mc <- mu_s10 + sqrt(sigma_s10) * z_10
          I_10 <- mean(log(dnorm(S1mc, mu_s10, sd = sqrt(sigma_s10)) / f_Delta_S(S1mc)))
          
          S2mc <- mu_s01 + sqrt(sigma_s01) * z_01
          I_01 <- mean(log(dnorm(S2mc, mu_s01, sd = sqrt(sigma_s01)) / f_Delta_S(S2mc)))
          
          S3mc <- mu_s00 + sqrt(sigma_s00) * z_00
          I_00 <- mean(log(((omega_1 * dnorm(S3mc, mu_s00, sd = sqrt(sigma_s00))) +
                              (omega_2 * dnorm(S3mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S3mc)))
          
          S4mc <- mu_s11 + sqrt(sigma_s11) * z_11
          I_11 <- mean(log(((omega_1 * dnorm(S4mc, mu_s00, sd = sqrt(sigma_s00))) +
                              (omega_2 * dnorm(S4mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S4mc)))
          
          # Calculate ICA
          I_Delta_T_Delta_S <- (pi_00_hier * I_00) + (pi_01_hier * I_01) + (pi_10_hier * I_10) + (pi_11_hier * I_11)
          H_Delta_T <- -((pi_Delta_T_min1 * log(pi_Delta_T_min1)) +
                           (pi_Delta_T_0 * log(pi_Delta_T_0)) +
                           (pi_Delta_T_1 * log(pi_Delta_T_1)))
          if(Monotonicity==TRUE){
            H_Delta_T <- -((pi_Delta_T_0 * log(pi_Delta_T_0)) +
                             (pi_Delta_T_1 * log(pi_Delta_T_1)))
          }
          R2_H <- I_Delta_T_Delta_S / H_Delta_T
          R2_C <- 1 - exp(-2*I_Delta_T_Delta_S)
          
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nR2H:", R2_H, "  R2C:", R2_C, "\n")
          }
          
          aantal <- aantal + 1
          if(Show.Details==TRUE){
            cat("\n", (aantal/M)*100, "% done. \n", sep="")
          }
          
          # save results (indexed assignment into pre-allocated vectors)
          R2_H_all[aantal] <- R2_H
          R2_C_all[aantal] <- R2_C
          G_rho_01_00_all[aantal] <- G_rho_01_00_hier
          G_rho_01_01_all[aantal] <- G_rho_01_01_hier
          G_rho_01_10_all[aantal] <- G_rho_01_10_hier
          G_rho_01_11_all[aantal] <- G_rho_01_11_hier
          pi_00_all[aantal] <- pi_00_hier
          pi_01_all[aantal] <- pi_01_hier
          pi_10_all[aantal] <- pi_10_hier
          pi_11_all[aantal] <- pi_11_hier
          pi_Delta_T_min1_all[aantal] <- pi_Delta_T_min1
          pi_Delta_T_0_all[aantal] <- pi_Delta_T_0
          pi_Delta_T_1_all[aantal] <- pi_Delta_T_1
          pi_0_00_e_all[aantal] <- pi_0_00_e
          pi_0_01_e_all[aantal] <- pi_0_01_e
          pi_0_10_e_all[aantal] <- pi_0_10_e
          pi_0_11_e_all[aantal] <- pi_0_11_e
          mu_0_00_all[aantal] <- mu_0_00
          mu_0_01_all[aantal] <- mu_0_01
          mu_0_10_all[aantal] <- mu_0_10
          mu_0_11_all[aantal] <- mu_0_11
          sigma_00_00_all[aantal] <- sigma_00_00
          sigma_00_01_all[aantal] <- sigma_00_01
          sigma_00_10_all[aantal] <- sigma_00_10
          sigma_00_11_all[aantal] <- sigma_00_11
          pi_1_00_e_all[aantal] <- pi_1_00_e
          pi_1_01_e_all[aantal] <- pi_1_01_e
          pi_1_10_e_all[aantal] <- pi_1_10_e
          pi_1_11_e_all[aantal] <- pi_1_11_e
          mu_1_00_all[aantal] <- mu_1_00
          mu_1_01_all[aantal] <- mu_1_01
          mu_1_10_all[aantal] <- mu_1_10
          mu_1_11_all[aantal] <- mu_1_11
          sigma_11_00_all[aantal] <- sigma_11_00
          sigma_11_01_all[aantal] <- sigma_11_01
          sigma_11_10_all[aantal] <- sigma_11_10
          sigma_11_11_all[aantal] <- sigma_11_11
          gamma_0_1_all[aantal] <- gamma_0_1
          gamma_0_2_all[aantal] <- gamma_0_2
          gamma_1_1_all[aantal] <- gamma_1_1
          gamma_1_2_all[aantal] <- gamma_1_2
          mean_Y_S0_all[aantal] <- mean_Y_S0
          mean_Y_S1_all[aantal] <- mean_Y_S1
          var_Y_S0_all[aantal] <- var_Y_S0
          var_Y_S1_all[aantal] <- var_Y_S1
          dev_S0_all[aantal] <- dev_S0
          dev_S1_all[aantal] <- dev_S1
          AIC_S0_all[aantal] <- AIC_S0
          AIC_S1_all[aantal] <- AIC_S1
          AIC_tot_all[aantal] <- AIC_tot
          code_nlm_0_all[aantal] <- code_nlm_0
          code_nlm_1_all[aantal] <- code_nlm_1
          
        } # min eigen and CN requirements
      } # mixture distribution available
    } # pi's requirements
  } # end while loop
  
  #  ---------------------------------------------------------------------------
  # 3. OUTPUT RESULTS ----------------------------------------------------
  # ----------------------------------------------------------------------------
  fit <- list(R2_H=R2_H_all, R2_C=R2_C_all,
              pi_00=as.numeric(pi_00_all), pi_01=as.numeric(pi_01_all),
              pi_10=as.numeric(pi_10_all), pi_11=as.numeric(pi_11_all),
              G_rho_01_00=as.numeric(G_rho_01_00_all),
              G_rho_01_01=as.numeric(G_rho_01_01_all),
              G_rho_01_10=as.numeric(G_rho_01_10_all),
              G_rho_01_11=as.numeric(G_rho_01_11_all),
              pi_Delta_T_min1=as.numeric(pi_Delta_T_min1_all),
              pi_Delta_T_0=as.numeric(pi_Delta_T_0_all),
              pi_Delta_T_1=as.numeric(pi_Delta_T_1_all),
              pi_0_00=as.numeric(pi_0_00_e_all), pi_0_01=as.numeric(pi_0_01_e_all),
              pi_0_10=as.numeric(pi_0_10_e_all), pi_0_11=as.numeric(pi_0_11_e_all),
              mu_0_00=as.numeric(mu_0_00_all), mu_0_01=as.numeric(mu_0_01_all),
              mu_0_10=as.numeric(mu_0_10_all), mu_0_11=as.numeric(mu_0_11_all),
              sigma2_00_00=as.numeric(sigma_00_00_all), sigma2_00_01=as.numeric(sigma_00_01_all),
              sigma2_00_10=as.numeric(sigma_00_10_all), sigma2_00_11=as.numeric(sigma_00_11_all),
              pi_1_00=as.numeric(pi_1_00_e_all), pi_1_01=as.numeric(pi_1_01_e_all),
              pi_1_10=as.numeric(pi_1_10_e_all), pi_1_11=as.numeric(pi_1_11_e_all),
              mu_1_00=as.numeric(mu_1_00_all), mu_1_01=as.numeric(mu_1_01_all),
              mu_1_10=as.numeric(mu_1_10_all), mu_1_11=as.numeric(mu_1_11_all),
              sigma2_11_00=as.numeric(sigma_11_00_all), sigma2_11_01=as.numeric(sigma_11_01_all),
              sigma2_11_10=as.numeric(sigma_11_10_all), sigma2_11_11=as.numeric(sigma_11_11_all),
              gamma_0_1=as.numeric(gamma_0_1_all),     # covariate-1 slope, control arm
              gamma_0_2=as.numeric(gamma_0_2_all),     # covariate-2 slope, control arm
              gamma_1_1=as.numeric(gamma_1_1_all),     # covariate-1 slope, experimental arm
              gamma_1_2=as.numeric(gamma_1_2_all),     # covariate-2 slope, experimental arm
              cov1_center=cov1_center,                 # centering constant, covariate 1
              cov2_center=cov2_center,                 # centering constant, covariate 2
              x_eval=x_eval,                           # covariate values used for the ICA
              mean_Y_S0=as.numeric(mean_Y_S0_all), mean_Y_S1=as.numeric(mean_Y_S1_all),
              var_Y_S0=as.numeric(var_Y_S0_all), var_Y_S1=as.numeric(var_Y_S1_all),
              dev_S0=as.numeric(dev_S0_all), dev_S1=as.numeric(dev_S1_all),
              AIC_S0=as.numeric(AIC_S0_all), AIC_S1=as.numeric(AIC_S1_all),
              AIC_tot=as.numeric(AIC_tot_all),       # joint AIC of the two arms
              code_nlm_0=as.numeric(code_nlm_0_all), code_nlm_1=as.numeric(code_nlm_1_all),
              mean.S0=mean(S_0, na.rm=TRUE), var.S0=var(S_0, na.rm=TRUE),
              mean.S1=mean(S_1, na.rm=TRUE), var.S1=var(S_1, na.rm=TRUE),
              Call=match.call())
  
  class(fit) <- "ICA.BinCont"
  fit
} # end of function


## =====================================================================
## ICA.BinCont with ONE CATEGORICAL covariate in the component means mu_ij
## ---------------------------------------------------------------------
##
## Model (homoscedastic case, Diff.Sigma = FALSE):
##   Control arm:       S0_k ~ sum_ij pi_ij * N( beta0_ij + sum_l gamma0_l * D_lk , sigma00 )
##   Experimental arm:  S1_k ~ sum_ij pi_ij * N( beta1_ij + sum_l gamma1_l * D_lk , sigma11 )
##   where D_1k,...,D_(K-1)k are the (K-1) dummy indicators of a K-level
##   categorical covariate, with the FIRST level used as reference (all
##   dummies = 0). gamma_l is the shift of level (l+1) relative to reference,
##   SHARED across the four T-strata within an arm.
##
## This is the multi-covariate linear predictor with p = K-1 slopes per arm. The intercepts
## beta_ij are the component means at the reference level.

ICA.BinCont.CovCat <- function(Dataset, Surr, True, Treat, Cov,
                               Ref.Level=NULL,
                               Eval.Level=NULL,
                               G_pi_10 = c(0,1),
                               G_rho_01_00=c(-1,1),
                               G_rho_01_01=c(-1,1),
                               G_rho_01_10=c(-1,1),
                               G_rho_01_11=c(-1,1),
                               Theta.S_0, # 4 intercepts + (K-1) slopes + sigma(s)
                               Theta.S_1,
                               M=1000,        # number of sensitivity-analysis draws
                               Seed=123,      # integer for standalone use; NULL inside an outer parallel loop
                               N.MC=10000,    # Monte Carlo size for the integrals I_ij
                               Center.Cov=FALSE,
                               Monotonicity=FALSE,
                               Independence=FALSE,
                               HAA=FALSE,
                               Cond_ind=FALSE,
                               Plots=TRUE, Save.Plots="No",
                               Show.Details=FALSE){
  # 1. PREPARATION --------------------------------------------------------------

  Diff.Sigma = FALSE # All 4 SD values in each group are the same
  
  ## EXTRACT VARIABLE ----
  Surr  <- Dataset[,paste(substitute(Surr))]
  True  <- Dataset[,paste(substitute(True))]
  Treat <- Dataset[,paste(substitute(Treat))]
  Covv  <- Dataset[,paste(substitute(Cov))]
  Covf  <- factor(Covv)                 # coerce categorical covariate to factor
  
  ## CHECK DATA FORMAT -------------------
  df <- data.frame(Surr=Surr, True=True, Treat=Treat, Cov=Covf,
                   stringsAsFactors = TRUE)
  data_no_miss.original <- df[stats::complete.cases(df), , drop=FALSE]
  data_no_miss.original$Cov <- droplevels(data_no_miss.original$Cov)  # drop empty levels
  if (min(na.exclude(Treat))!=c(-1)) {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (max(na.exclude(Treat))!=c(1))  {stop("\nTreatment should be coded as -1=control and 1=experimental treatment.")}
  if (length(unique(na.exclude(True)))>2) {stop("\nThe true endpoint should be binary.")}
  if (min(na.exclude(True))!=c(0)) {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  if (max(na.exclude(True))!=c(1))  {stop("\nThe true endpoint should be coded as 0=no response and 1=response.")}
  if (nlevels(data_no_miss.original$Cov) < 2) {stop("\nThe categorical covariate must have at least 2 non-empty levels.")}
  
  ## SET REFERENCE LEVEL ------------------
  if (!is.null(Ref.Level)){
    if (!(Ref.Level %in% levels(data_no_miss.original$Cov)))
      stop("\nRef.Level is not a level of the categorical covariate.")
    data_no_miss.original$Cov <- relevel(data_no_miss.original$Cov, ref = as.character(Ref.Level))
  }
  lev <- levels(data_no_miss.original$Cov)     # lev[1] = reference
  nd  <- length(lev) - 1                        # number of dummies = K-1
  
  ## DEFAULT EVALUATION LEVEL -------------
  if (is.null(Eval.Level)) Eval.Level <- lev[1]        # reference by default
  if (!(as.character(Eval.Level) %in% lev))
    stop("\nEval.Level is not a level of the categorical covariate.")
  
  ## BUILD DUMMY DESIGN MATRIX ------------
  # model.matrix(~ Cov)[, -1] gives the (K-1) reference-coded dummies, in the
  # order lev[2], lev[3], ..., lev[K].
  Dfull <- stats::model.matrix(~ Cov, data = data_no_miss.original)[, -1, drop=FALSE]
  colnames(Dfull) <- lev[-1]
  dummy_center <- if (Center.Cov) colMeans(Dfull) else rep(0, nd)
  Dfull_c <- sweep(Dfull, 2, dummy_center, "-")        # centered (or raw) dummies
  
  # Dummy row for the evaluation level (same scale as the fitted design)
  d_eval <- as.numeric(lev[-1] == as.character(Eval.Level)) - dummy_center  # length nd
  
  ## CHECK LENGTH OF STARTING VALUES ------
  need <- if (Diff.Sigma) (4 + nd + 4) else (4 + nd + 1)
  if (length(Theta.S_0)!=need | length(Theta.S_1)!=need){
    stop(paste0("\nWith a ", nd+1, "-level categorical covariate, Theta.S_0 and ",
                "Theta.S_1 must have ", need, " values (4 intercepts + ", nd,
                " slopes", if (Diff.Sigma) " + 4 sigmas)." else " + 1 sigma)."))
  }
  
  ## SEED SETUP ---------------------------
  # Scientific stream: set once here, then left to flow through the loop.
  if (!is.null(Seed)) set.seed(Seed)
  # Make sure a random seed already exist, if not create one
  if (!exists(".Random.seed", envir = .GlobalEnv)) stats::runif(1)
  # Store the current state of the random seed to use later
  .sci_state <- get(".Random.seed", envir = .GlobalEnv)
  
  # Fixed N(0,1) integration patterns (reproduce the original seeds 1..4),
  # drawn once and held constant for every iteration. The scientific stream
  # saved above is restored immediately after, so these draws do not perturb it.
  set.seed(1L); z_10 <- rnorm(N.MC)
  set.seed(2L); z_01 <- rnorm(N.MC)
  set.seed(3L); z_00 <- rnorm(N.MC)
  set.seed(4L); z_11 <- rnorm(N.MC)
  
  assign(".Random.seed", .sci_state, envir = .GlobalEnv)
  
  ## BUILD FUNCTIONS TO CALCULATE MIXTURE DEVIANCE -----------------------------
  # The mixing weights pi_ij are passed in as arguments (p00..p11).
  # covarmat is the (n x nd) dummy design; gamma = theta[5:(4+nd)].
  if (Diff.Sigma==TRUE){
    # Unequal sigma's: theta = c(b00,b01,b10,b11, gamma_1..gamma_nd, s00,s01,s10,s11)
    mixt.deviance <- function(theta, data, covarmat, p00, p01, p10, p11){
      p_ <- ncol(covarmat)
      gamma <- theta[5:(4+p_)]
      s <- theta[(5+p_):(8+p_)]
      gc <- as.vector(covarmat %*% gamma)   # sum_l gamma_l * D_l (shared by all 4 means)
      pdf <- p00*dnorm(data, theta[1]+gc, s[1]) +
        p01*dnorm(data, theta[2]+gc, s[2]) +
        p10*dnorm(data, theta[3]+gc, s[3]) +
        p11*dnorm(data, theta[4]+gc, s[4])
      -2*sum(log(pdf))
    }
  } else {
    ## Equal sigma's: theta = c(b00,b01,b10,b11, gamma_1..gamma_nd, sigma)
    mixt.deviance <- function(theta, data, covarmat, p00, p01, p10, p11){
      p_ <- ncol(covarmat)
      gamma <- theta[5:(4+p_)]
      sigma_all <- theta[5+p_]
      gc <- as.vector(covarmat %*% gamma)   # sum_l gamma_l * D_l (shared by all 4 means)
      pdf <- p00*dnorm(data, theta[1]+gc, sigma_all) +
        p01*dnorm(data, theta[2]+gc, sigma_all) +
        p10*dnorm(data, theta[3]+gc, sigma_all) +
        p11*dnorm(data, theta[4]+gc, sigma_all)
      -2*sum(log(pdf))
    }
  }
  
  ## PRE-ALLOCATE RESULT ACCUMULATORS -------------------------------------------
  R2_H_all <- R2_C_all <- numeric(M)
  pi_00_all <- pi_01_all <- pi_10_all <- pi_11_all <- numeric(M)
  G_rho_01_00_all <- G_rho_01_01_all <- G_rho_01_10_all <- G_rho_01_11_all <- numeric(M)
  pi_Delta_T_min1_all <- pi_Delta_T_0_all <- pi_Delta_T_1_all <- numeric(M)
  pi_0_00_e_all <- pi_0_01_e_all <- pi_0_10_e_all <- pi_0_11_e_all <- numeric(M)
  mu_0_00_all <- mu_0_01_all <- mu_0_10_all <- mu_0_11_all <- numeric(M)
  sigma_00_00_all <- sigma_00_01_all <- sigma_00_10_all <- sigma_00_11_all <- numeric(M)
  pi_1_00_e_all <- pi_1_01_e_all <- pi_1_10_e_all <- pi_1_11_e_all <- numeric(M)
  mu_1_00_all <- mu_1_01_all <- mu_1_10_all <- mu_1_11_all <- numeric(M)
  sigma_11_00_all <- sigma_11_01_all <- sigma_11_10_all <- sigma_11_11_all <- numeric(M)
  gamma_0_all <- matrix(0, nrow=M, ncol=nd, dimnames=list(NULL, lev[-1]))   # slopes, control
  gamma_1_all <- matrix(0, nrow=M, ncol=nd, dimnames=list(NULL, lev[-1]))   # slopes, experimental
  mean_Y_S0_all <- mean_Y_S1_all <- numeric(M)
  var_Y_S0_all <- var_Y_S1_all <- numeric(M)
  dev_S0_all <- dev_S1_all <- numeric(M)
  AIC_S0_all <- AIC_S1_all <- AIC_tot_all <- numeric(M)   # AIC per arm and combined
  code_nlm_0_all <- code_nlm_1_all <- numeric(M)
  aantal <- 0
  
  ## COMPUTE PER-ARM VECTORS AND MARGINAL PROBABILITIES ------------------------
  dnm <- data_no_miss.original
  ctrl <- dnm$Treat == -1                # control arm rows
  expt <- dnm$Treat ==  1                # experimental arm rows
  S_0 <- dnm$Surr[ctrl]
  S_1 <- dnm$Surr[expt]
  D_0 <- Dfull_c[ctrl, , drop=FALSE]     # dummy design, control arm
  D_1 <- Dfull_c[expt, , drop=FALSE]     # dummy design, experimental arm
  
  # Estimate the marginal probabilities of T
  pi_punt1 <- mean(dnm$True[expt])  # E(T|Z=1)
  pi_punt0 <- 1 - pi_punt1
  pi_1punt <- mean(dnm$True[ctrl])  # E(T|Z=0)
  pi_0punt <- 1 - pi_1punt
  if (Monotonicity==TRUE){
    if (pi_punt0 > pi_0punt | pi_1punt > pi_punt1) {stop("\nThe monotonicity assumption is not compatible with the data")}
    G_pi_10 = c(0,0)
  }
  
  # Arm-specific dummy means (proportions) and covariance matrices (for GOF)
  m0 <- colMeans(D_0); CovD0 <- stats::cov(D_0)   # control
  m1 <- colMeans(D_1); CovD1 <- stats::cov(D_1)   # experimental
  
  #  ---------------------------------------------------------------------------
  # 2. SENSITIVITY ANALYSIS ----------------------------------------------------
  # ----------------------------------------------------------------------------
  
  
  while (aantal < M){
    
    ## 2.1. Sample the unidentifiable parameters -------------------------------
    ## (no per-iteration set.seed; the stream set up above advances naturally)
    if(HAA==FALSE){
      G_rho_01_00_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
      G_rho_01_01_hier <- runif(1, min=G_rho_01_01[1], max=G_rho_01_01[2])
      G_rho_01_10_hier <- runif(1, min=G_rho_01_10[1], max=G_rho_01_10[2])
      G_rho_01_11_hier <- runif(1, min=G_rho_01_11[1], max=G_rho_01_11[2])
    }
    if(HAA==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- runif(1, min=G_rho_01_00[1], max=G_rho_01_00[2])
    }
    if(Cond_ind==TRUE){
      G_rho_01_00_hier <- G_rho_01_01_hier <- G_rho_01_10_hier <- G_rho_01_11_hier <- 0
    }
    if(Show.Details==TRUE){
      flush.console()
      cat("\n\nG_rho:", cbind(G_rho_01_00_hier, G_rho_01_01_hier, G_rho_01_10_hier, G_rho_01_11_hier), "\n")
    }
    
    # Sample pi_10 and derive the remaining pi_ij compatible with the margins
    if(Independence==FALSE){
      pi_10_hier <- runif(1, min=G_pi_10[1], max=G_pi_10[2])
      pi_00_hier <- pi_punt0 - pi_10_hier
      pi_01_hier <- pi_0punt - pi_00_hier
      pi_11_hier <- pi_1punt - pi_10_hier
    }
    if(Independence==TRUE){
      pi_10_hier <- pi_1punt*pi_punt0
      pi_00_hier <- pi_0punt*pi_punt0
      pi_01_hier <- pi_0punt*pi_punt1
      pi_11_hier <- pi_1punt*pi_punt1
    }
    
    if ((pi_00_hier >= 0 & pi_01_hier >= 0 & pi_10_hier >= 0 & pi_11_hier >= 0) &
        (pi_00_hier <= 1 & pi_01_hier <= 1 & pi_10_hier <= 1 & pi_11_hier <= 1)){
      
      if(Show.Details==TRUE){
        flush.console()
        cat("\n\nPi's used here:", cbind(pi_00_hier, pi_01_hier, pi_10_hier, pi_11_hier), "\n")
      }
      
      ## 2.2. Calculate pi_delta_T -----------------------------------------------
      # individual causal effect of the treatment on the true endpoint: delta_T = T1 - T0
      pi_Delta_T_min1 <- pi_10_hier # Delta_T = -1
      pi_Delta_T_0    <- pi_00_hier + pi_11_hier # Delta_T = 0
      pi_Delta_T_1    <- pi_01_hier # Delta_T = 1
      
      ## 2.3. Estimate Gaussian mixture models of Y_S| Y_T -----------------------
      
      # Fit model: maximize likelihood = minimizing deviance
      # Iterative process: start with an initial guess of parameter values,
      # then iteratively updat them to minimize deviance, using Newton-Raphson type algorithm
      
      options(warn = -999)
      mix1 <- nlm(f = mixt.deviance, p = Theta.S_0, data = S_0, covarmat = D_0,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000 #the maximum number of iterations to be performed before the program is terminated
      )
      mix2 <- nlm(f = mixt.deviance, p = Theta.S_1, data = S_1, covarmat = D_1,
                  p00 = pi_00_hier, p01 = pi_01_hier, p10 = pi_10_hier, p11 = pi_11_hier,
                  iterlim = 100000)
      options(warn = 1)
      
      if (exists("mix1")==TRUE & exists("mix2")==TRUE & is.null(mix1)==FALSE & is.null(mix2)==FALSE &
          mix1$code<=3 & mix2$code<=3){
        
        ## 2.4. Extract fitted parameters (intercepts, slopes, residual variances) & deviance ----
        # Control arm (S0)
        mu_0_00 <- mix1$estimate[1]; mu_0_01 <- mix1$estimate[2]
        mu_0_10 <- mix1$estimate[3]; mu_0_11 <- mix1$estimate[4]
        gamma_0 <- mix1$estimate[5:(4+nd)]           # length-nd slope vector
        if (Diff.Sigma==FALSE){ sigma_00_00 <- sigma_00_01 <- sigma_00_10 <- sigma_00_11 <- mix1$estimate[5+nd]**2 }
        if (Diff.Sigma==TRUE){
          sigma_00_00 <- mix1$estimate[5+nd]**2; sigma_00_01 <- mix1$estimate[6+nd]**2
          sigma_00_10 <- mix1$estimate[7+nd]**2; sigma_00_11 <- mix1$estimate[8+nd]**2 }
        pi_0_00_e <- pi_00_hier; pi_0_01_e <- pi_01_hier; pi_0_10_e <- pi_10_hier; pi_0_11_e <- pi_11_hier
        code_nlm_0 <- mix1$code
        
        # Experimental arm (S1)
        mu_1_00 <- mix2$estimate[1]; mu_1_01 <- mix2$estimate[2]
        mu_1_10 <- mix2$estimate[3]; mu_1_11 <- mix2$estimate[4]
        gamma_1 <- mix2$estimate[5:(4+nd)]           # length-nd slope vector
        if (Diff.Sigma==FALSE){ sigma_11_00 <- sigma_11_01 <- sigma_11_10 <- sigma_11_11 <- mix2$estimate[5+nd]**2 }
        if (Diff.Sigma==TRUE){
          sigma_11_00 <- mix2$estimate[5+nd]**2; sigma_11_01 <- mix2$estimate[6+nd]**2
          sigma_11_10 <- mix2$estimate[7+nd]**2; sigma_11_11 <- mix2$estimate[8+nd]**2 }
        pi_1_00_e <- pi_00_hier; pi_1_01_e <- pi_01_hier; pi_1_10_e <- pi_10_hier; pi_1_11_e <- pi_11_hier
        code_nlm_1 <- mix2$code
        
        # Conditional deviance of each fitted mixture-of-regressions (for GOF check)
        dev_S0 <- mix1$minimum
        dev_S1 <- mix2$minimum
        
        # AIC = deviance + 2*k. k = number of parameters ESTIMATED per arm
        # (4 intercepts + slope + sigma(s) = length of Theta.S, i.e. `need`).
        # The mixing weights pi_ij are fixed inputs of the sensitivity draw,
        # not estimated by MLE, so they are NOT counted in k. Each arm is a
        # separate fit; AIC_tot is the joint criterion for the two independent
        # arms (deviances add, parameters add).
        k_arm   <- need
        AIC_S0  <- dev_S0 + 2*k_arm
        AIC_S1  <- dev_S1 + 2*k_arm
        AIC_tot <- (dev_S0 + dev_S1) + 2*(2*k_arm)
        
        if (Plots == TRUE){
          ## Marginal fitted density, averaged over the observed dummy rows,
          ## overlaid on the marginal histogram of S (O(n^2), so only when plotting).
          sd_0 <- sqrt(sigma_00_00); sd_1 <- sqrt(sigma_11_00)
          comp_dens <- function(grid, b, gamma, covarmat, sdv){
            gc <- as.vector(covarmat %*% gamma)
            sapply(grid, function(s) mean(dnorm(s, b + gc, sdv)))
          }
          x_S_0 <- sort(S_0)
          S0_M1 <- pi_00_hier * comp_dens(x_S_0, mu_0_00, gamma_0, D_0, sd_0)
          S0_M2 <- pi_01_hier * comp_dens(x_S_0, mu_0_01, gamma_0, D_0, sd_0)
          S0_M3 <- pi_10_hier * comp_dens(x_S_0, mu_0_10, gamma_0, D_0, sd_0)
          S0_M4 <- pi_11_hier * comp_dens(x_S_0, mu_0_11, gamma_0, D_0, sd_0)
          Sum_M_S_0 <- S0_M1 + S0_M2 + S0_M3 + S0_M4
          x_S_1 <- sort(S_1)
          S1_M1 <- pi_00_hier * comp_dens(x_S_1, mu_1_00, gamma_1, D_1, sd_1)
          S1_M2 <- pi_01_hier * comp_dens(x_S_1, mu_1_01, gamma_1, D_1, sd_1)
          S1_M3 <- pi_10_hier * comp_dens(x_S_1, mu_1_10, gamma_1, D_1, sd_1)
          S1_M4 <- pi_11_hier * comp_dens(x_S_1, mu_1_11, gamma_1, D_1, sd_1)
          Sum_M_S_1 <- S1_M1 + S1_M2 + S1_M3 + S1_M4
          
          plot_one <- function(){
            par(mfrow=c(1,2))
            max_val <- max(max(density(S_0)$y), max(Sum_M_S_0))
            hist(S_0, freq = F, ylim=c(0, max_val), xlab=expression(S[0]),
                 main=bquote(paste(S[0], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_0, Sum_M_S_0, col="red", lwd=4, lty=1)
            lines(x_S_0, S0_M1, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M2, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M3, col="blue", lwd=2, lty=3)
            lines(x_S_0, S0_M4, col="blue", lwd=2, lty=3)
            max_val <- max(max(density(S_1)$y), max(Sum_M_S_1))
            hist(S_1, freq = F, ylim=c(0, max_val), xlab=expression(S[1]),
                 main=bquote(paste(S[1], " (run "~.(aantal+1), ")")), col="grey")
            lines(x_S_1, Sum_M_S_1, col="red", lwd=4, lty=1)
            lines(x_S_1, S1_M1, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M2, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M3, col="blue", lwd=2, lty=3)
            lines(x_S_1, S1_M4, col="blue", lwd=2, lty=3)
            par(mfrow=c(1,1))
          }
          if (Save.Plots!="No"){
            options(warn=-1); dir.create(Save.Plots); options(warn=0)
            pdf(paste(Save.Plots, "/Run", (aantal+1), ".pdf", sep=""),width=6,height=4,paper='special')
            plot_one(); dev.off()
          }
          plot_one()
        }
        
        ## 2.5 Covariance-matrix checks (PD and condition number) --------------
        # If checks do not pass, the results are discared (so we always have M valid results)
        mat_a <- matrix(c(sigma_00_00,
                          (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                          (G_rho_01_00_hier*(sqrt(sigma_00_00*sigma_11_00))),
                          sigma_11_00), nrow = 2, byrow = TRUE)
        eigen_mat_a <- min(eigen((mat_a))$values)
        mat_b <- matrix(c(sigma_00_01,
                          (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                          (G_rho_01_01_hier*(sqrt(sigma_00_01*sigma_11_01))),
                          sigma_11_01), nrow = 2, byrow = TRUE)
        eigen_mat_b <- min(eigen((mat_b))$values)
        mat_c <- matrix(c(sigma_00_10,
                          (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                          (G_rho_01_10_hier*(sqrt(sigma_00_10*sigma_11_10))),
                          sigma_11_10), nrow = 2, byrow = TRUE)
        eigen_mat_c <- min(eigen((mat_c))$values)
        mat_d <- matrix(c(sigma_00_11,
                          (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                          (G_rho_01_11_hier*(sqrt(sigma_00_11*sigma_11_11))),
                          sigma_11_11), nrow = 2, byrow = TRUE)
        eigen_mat_d <- min(eigen((mat_d))$values)
        
        CN_mat_a <- max(svd(mat_a)$d)/min(svd(mat_a)$d)
        CN_mat_b <- max(svd(mat_b)$d)/min(svd(mat_b)$d)
        CN_mat_c <- max(svd(mat_c)$d)/min(svd(mat_c)$d)
        CN_mat_d <- max(svd(mat_d)$d)/min(svd(mat_d)$d)
        
        if (eigen_mat_a > 0 & eigen_mat_b > 0 & eigen_mat_c > 0 & eigen_mat_d > 0 &
            CN_mat_a < 50 & CN_mat_b < 50 & CN_mat_c < 50 & CN_mat_d < 50){
          
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nMin Eigenvalues:", cbind(eigen_mat_a, eigen_mat_b, eigen_mat_c, eigen_mat_d), "\n")
            cat("\n\nCondition number:", cbind(CN_mat_a, CN_mat_b, CN_mat_c, CN_mat_d), "\n")
            cat("\n\nChecks passed...")
          }
          
          ## 2.6 Marginal moments of S0, S1 (for GOF check) --------------------
          ## Categorical covariate contribution via the dummy design:
          ##   mean:      sum_l gamma_l * qbar_l         (qbar = arm dummy proportions)
          ##   variance:  t(gamma) %*% Cov(D) %*% gamma  (between-category variance)
          mean_int_S0 <- (pi_0_00_e*mu_0_00) + (pi_0_01_e*mu_0_01) + (pi_0_10_e*mu_0_10) + (pi_0_11_e*mu_0_11)
          mean_int_S1 <- (pi_1_00_e*mu_1_00) + (pi_1_01_e*mu_1_01) + (pi_1_10_e*mu_1_10) + (pi_1_11_e*mu_1_11)
          
          mean_Y_S0 <- mean_int_S0 + sum(gamma_0 * m0)   # E(S | Z=0)
          mean_Y_S1 <- mean_int_S1 + sum(gamma_1 * m1)   # E(S | Z=1)
          
          covvar_S0 <- as.numeric(t(gamma_0) %*% CovD0 %*% gamma_0)   # Var(sum_l gamma0_l D_l)
          covvar_S1 <- as.numeric(t(gamma_1) %*% CovD1 %*% gamma_1)   # Var(sum_l gamma1_l D_l)
          
          var_Y_S0 <- ((pi_0_00_e*mu_0_00^2)+(pi_0_01_e*mu_0_01^2)+(pi_0_10_e*mu_0_10^2)+(pi_0_11_e*mu_0_11^2)) - mean_int_S0^2 +
            ((pi_0_00_e*sigma_00_00)+(pi_0_01_e*sigma_00_01)+(pi_0_10_e*sigma_00_10)+(pi_0_11_e*sigma_00_11)) +
            covvar_S0
          
          var_Y_S1 <- ((pi_1_00_e*mu_1_00^2)+(pi_1_01_e*mu_1_01^2)+(pi_1_10_e*mu_1_10^2)+(pi_1_11_e*mu_1_11^2)) - mean_int_S1^2 +
            ((pi_1_00_e*sigma_11_00)+(pi_1_01_e*sigma_11_01)+(pi_1_10_e*sigma_11_10)+(pi_1_11_e*sigma_11_11)) +
            covvar_S1
          
          ## 2.7 Calculate ICA via Monte Carlo ---------------------------------
          ## Treatment effect on the surrogate within each T-stratum, at level Eval.Level.
          ## mu*_ij = (beta1_ij - beta0_ij) + sum_l (gamma1_l - gamma0_l) * d_eval_l
          ## NOTE: sum_l (gamma1_l - gamma0_l) d_eval_l is the SAME constant for all
          ## four strata (shared slopes), a common location shift, so the chosen
          ## Eval.Level does not change R2_H / R2_C.
          delta_cov <- sum((gamma_1 - gamma_0) * d_eval)
          mu_s00 <- (mu_1_00 - mu_0_00) + delta_cov
          mu_s01 <- (mu_1_01 - mu_0_01) + delta_cov
          mu_s10 <- (mu_1_10 - mu_0_10) + delta_cov
          mu_s11 <- (mu_1_11 - mu_0_11) + delta_cov
          
          # Weights of pi_00 and pi_11
          omega_1 <- pi_00_hier / (pi_00_hier + pi_11_hier); omega_2 <- 1 - omega_1
          
          sigma_s00 <- sigma_00_00 + sigma_11_00 - (2 * (sqrt(sigma_00_00 * sigma_11_00) * G_rho_01_00_hier))
          sigma_s01 <- sigma_00_01 + sigma_11_01 - (2 * (sqrt(sigma_00_01 * sigma_11_01) * G_rho_01_01_hier))
          sigma_s10 <- sigma_00_10 + sigma_11_10 - (2 * (sqrt(sigma_00_10 * sigma_11_10) * G_rho_01_10_hier))
          sigma_s11 <- sigma_00_11 + sigma_11_11 - (2 * (sqrt(sigma_00_11 * sigma_11_11) * G_rho_01_11_hier))
          
          f_Delta_S <- function(val){
            (pi_00_hier * dnorm(val, mu_s00, sd = sqrt(sigma_s00))) +
              (pi_01_hier * dnorm(val, mu_s01, sd = sqrt(sigma_s01))) +
              (pi_10_hier * dnorm(val, mu_s10, sd = sqrt(sigma_s10))) +
              (pi_11_hier * dnorm(val, mu_s11, sd = sqrt(sigma_s11)))
          }
          
          ## Sample S1, S2, S3, S4 from normal distribution
          ## S = mu + sd * z  is identical to rnorm(n, mu, sd) under the same seed.
          ## z values were sampled from N(0,1) in the SEED SETUP block above.
          S1mc <- mu_s10 + sqrt(sigma_s10) * z_10
          I_10 <- mean(log(dnorm(S1mc, mu_s10, sd = sqrt(sigma_s10)) / f_Delta_S(S1mc)))
          
          S2mc <- mu_s01 + sqrt(sigma_s01) * z_01
          I_01 <- mean(log(dnorm(S2mc, mu_s01, sd = sqrt(sigma_s01)) / f_Delta_S(S2mc)))
          
          S3mc <- mu_s00 + sqrt(sigma_s00) * z_00
          I_00 <- mean(log(((omega_1 * dnorm(S3mc, mu_s00, sd = sqrt(sigma_s00))) +
                              (omega_2 * dnorm(S3mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S3mc)))
          
          S4mc <- mu_s11 + sqrt(sigma_s11) * z_11
          I_11 <- mean(log(((omega_1 * dnorm(S4mc, mu_s00, sd = sqrt(sigma_s00))) +
                              (omega_2 * dnorm(S4mc, mu_s11, sd = sqrt(sigma_s11)))) / f_Delta_S(S4mc)))
          
          # Calculate ICA
          I_Delta_T_Delta_S <- (pi_00_hier * I_00) + (pi_01_hier * I_01) + (pi_10_hier * I_10) + (pi_11_hier * I_11)
          H_Delta_T <- -((pi_Delta_T_min1 * log(pi_Delta_T_min1)) +
                           (pi_Delta_T_0 * log(pi_Delta_T_0)) +
                           (pi_Delta_T_1 * log(pi_Delta_T_1)))
          if(Monotonicity==TRUE){
            H_Delta_T <- -((pi_Delta_T_0 * log(pi_Delta_T_0)) +
                             (pi_Delta_T_1 * log(pi_Delta_T_1)))
          }
          R2_H <- I_Delta_T_Delta_S / H_Delta_T
          R2_C <- 1 - exp(-2*I_Delta_T_Delta_S)
          
          if(Show.Details==TRUE){
            flush.console()
            cat("\n\nR2H:", R2_H, "  R2C:", R2_C, "\n")
          }
          
          aantal <- aantal + 1
          if(Show.Details==TRUE){
            cat("\n", (aantal/M)*100, "% done. \n", sep="")
          }
          
          # save results (indexed assignment into pre-allocated vectors)
          R2_H_all[aantal] <- R2_H
          R2_C_all[aantal] <- R2_C
          G_rho_01_00_all[aantal] <- G_rho_01_00_hier
          G_rho_01_01_all[aantal] <- G_rho_01_01_hier
          G_rho_01_10_all[aantal] <- G_rho_01_10_hier
          G_rho_01_11_all[aantal] <- G_rho_01_11_hier
          pi_00_all[aantal] <- pi_00_hier
          pi_01_all[aantal] <- pi_01_hier
          pi_10_all[aantal] <- pi_10_hier
          pi_11_all[aantal] <- pi_11_hier
          pi_Delta_T_min1_all[aantal] <- pi_Delta_T_min1
          pi_Delta_T_0_all[aantal] <- pi_Delta_T_0
          pi_Delta_T_1_all[aantal] <- pi_Delta_T_1
          pi_0_00_e_all[aantal] <- pi_0_00_e
          pi_0_01_e_all[aantal] <- pi_0_01_e
          pi_0_10_e_all[aantal] <- pi_0_10_e
          pi_0_11_e_all[aantal] <- pi_0_11_e
          mu_0_00_all[aantal] <- mu_0_00
          mu_0_01_all[aantal] <- mu_0_01
          mu_0_10_all[aantal] <- mu_0_10
          mu_0_11_all[aantal] <- mu_0_11
          sigma_00_00_all[aantal] <- sigma_00_00
          sigma_00_01_all[aantal] <- sigma_00_01
          sigma_00_10_all[aantal] <- sigma_00_10
          sigma_00_11_all[aantal] <- sigma_00_11
          pi_1_00_e_all[aantal] <- pi_1_00_e
          pi_1_01_e_all[aantal] <- pi_1_01_e
          pi_1_10_e_all[aantal] <- pi_1_10_e
          pi_1_11_e_all[aantal] <- pi_1_11_e
          mu_1_00_all[aantal] <- mu_1_00
          mu_1_01_all[aantal] <- mu_1_01
          mu_1_10_all[aantal] <- mu_1_10
          mu_1_11_all[aantal] <- mu_1_11
          sigma_11_00_all[aantal] <- sigma_11_00
          sigma_11_01_all[aantal] <- sigma_11_01
          sigma_11_10_all[aantal] <- sigma_11_10
          sigma_11_11_all[aantal] <- sigma_11_11
          gamma_0_all[aantal, ] <- gamma_0
          gamma_1_all[aantal, ] <- gamma_1
          mean_Y_S0_all[aantal] <- mean_Y_S0
          mean_Y_S1_all[aantal] <- mean_Y_S1
          var_Y_S0_all[aantal] <- var_Y_S0
          var_Y_S1_all[aantal] <- var_Y_S1
          dev_S0_all[aantal] <- dev_S0
          dev_S1_all[aantal] <- dev_S1
          AIC_S0_all[aantal] <- AIC_S0
          AIC_S1_all[aantal] <- AIC_S1
          AIC_tot_all[aantal] <- AIC_tot
          code_nlm_0_all[aantal] <- code_nlm_0
          code_nlm_1_all[aantal] <- code_nlm_1
          
        } # min eigen and CN requirements
      } # mixture distribution available
    } # pi's requirements
  } # end while loop
  
  #  ---------------------------------------------------------------------------
  # 3. OUTPUT RESULTS ----------------------------------------------------
  # ----------------------------------------------------------------------------
  fit <- list(R2_H=R2_H_all, R2_C=R2_C_all,
              pi_00=as.numeric(pi_00_all), pi_01=as.numeric(pi_01_all),
              pi_10=as.numeric(pi_10_all), pi_11=as.numeric(pi_11_all),
              G_rho_01_00=as.numeric(G_rho_01_00_all),
              G_rho_01_01=as.numeric(G_rho_01_01_all),
              G_rho_01_10=as.numeric(G_rho_01_10_all),
              G_rho_01_11=as.numeric(G_rho_01_11_all),
              pi_Delta_T_min1=as.numeric(pi_Delta_T_min1_all),
              pi_Delta_T_0=as.numeric(pi_Delta_T_0_all),
              pi_Delta_T_1=as.numeric(pi_Delta_T_1_all),
              pi_0_00=as.numeric(pi_0_00_e_all), pi_0_01=as.numeric(pi_0_01_e_all),
              pi_0_10=as.numeric(pi_0_10_e_all), pi_0_11=as.numeric(pi_0_11_e_all),
              mu_0_00=as.numeric(mu_0_00_all), mu_0_01=as.numeric(mu_0_01_all),
              mu_0_10=as.numeric(mu_0_10_all), mu_0_11=as.numeric(mu_0_11_all),
              sigma2_00_00=as.numeric(sigma_00_00_all), sigma2_00_01=as.numeric(sigma_00_01_all),
              sigma2_00_10=as.numeric(sigma_00_10_all), sigma2_00_11=as.numeric(sigma_00_11_all),
              pi_1_00=as.numeric(pi_1_00_e_all), pi_1_01=as.numeric(pi_1_01_e_all),
              pi_1_10=as.numeric(pi_1_10_e_all), pi_1_11=as.numeric(pi_1_11_e_all),
              mu_1_00=as.numeric(mu_1_00_all), mu_1_01=as.numeric(mu_1_01_all),
              mu_1_10=as.numeric(mu_1_10_all), mu_1_11=as.numeric(mu_1_11_all),
              sigma2_11_00=as.numeric(sigma_11_00_all), sigma2_11_01=as.numeric(sigma_11_01_all),
              sigma2_11_10=as.numeric(sigma_11_10_all), sigma2_11_11=as.numeric(sigma_11_11_all),
              gamma_0=gamma_0_all,                     # (M x nd) dummy slopes, control arm
              gamma_1=gamma_1_all,                     # (M x nd) dummy slopes, experimental arm
              levels=lev,                              # factor levels, lev[1] = reference
              ref_level=lev[1],                        # reference level
              eval_level=as.character(Eval.Level),     # level used for the ICA evaluation
              dummy_center=dummy_center,               # centering constants for the dummies
              mean_Y_S0=as.numeric(mean_Y_S0_all), mean_Y_S1=as.numeric(mean_Y_S1_all),
              var_Y_S0=as.numeric(var_Y_S0_all), var_Y_S1=as.numeric(var_Y_S1_all),
              dev_S0=as.numeric(dev_S0_all), dev_S1=as.numeric(dev_S1_all),
              AIC_S0=as.numeric(AIC_S0_all), AIC_S1=as.numeric(AIC_S1_all),
              AIC_tot=as.numeric(AIC_tot_all),       # joint AIC of the two arms
              npar=need,                             # estimated parameters PER arm (k)
              code_nlm_0=as.numeric(code_nlm_0_all), code_nlm_1=as.numeric(code_nlm_1_all),
              mean.S0=mean(S_0, na.rm=TRUE), var.S0=var(S_0, na.rm=TRUE),
              mean.S1=mean(S_1, na.rm=TRUE), var.S1=var(S_1, na.rm=TRUE),
              Call=match.call())
  
  class(fit) <- "ICA.BinCont"
  fit
} # end of function
