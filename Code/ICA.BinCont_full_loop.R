# P_INFL_OO6                    ################################################
# Script: ICA_Bincont_full_loop ################################################
# Author: Trang Tu              ################################################


## =====================================================================
## NESTED BOOTSTRAP + SENSITIVITY ANALYSIS for ICA.BinCont.Cov
## ---------------------------------------------------------------------
## Implements the workflow of Alonso Abad et al. (Stat Med 2024):
##   OUTER loop : Nboot stratified bootstrap resamples of the data (imprecision)
##   INNER call : the full M-draw sensitivity sweep on each resample (ignorance)
##   POOL       : stack all Nboot * M ICA values into one frequency distribution
##
## The OUTER loop is parallelised with foreach() %dofuture%. The INNER loop is
## ICA.BinCont.Cov() from ICA.BinCont.Cov_inner.R.
##
## SEED DESIGN
##   * One global seed governs everything: set.seed(global_seed) before the
##     foreach, plus .options.future = list(seed = TRUE), gives each bootstrap
##     its own independent, reproducible L'Ecuyer-CMRG stream.
##   * The inner function is called with Seed = NULL, so it does NOT reseed and
##     instead draws the resample AND the theta_U draws from that worker stream.
##   * The inner function's four fixed N(0,1) integration patterns are generated
##     with a save/restore of the stream, so they do not disturb it.
##
## COMBINING THE INNER OUTPUT
##   Each inner call returns a list whose elements have different lengths
##   (length-M vectors, a few scalars, and Call). fit_to_df() reshapes one such
##   list into a tidy M-row data frame (scalars recycled, Call dropped, a `boot`
##   id and a `draw` id added). The outer loop returns one data frame per
##   bootstrap; do.call(rbind, .) stacks them into a single Nboot*M-row table.
## =====================================================================

library(pacman)
pacman::p_load(
  foreach,     # foreach loop
  rio,         # Export file
  here,
  future,
  doFuture,
  Surrogate,
  progressr,   # progress reporting that works with parallel futures
  progress     # back-end used by progressr::handler_progress for the format string
)

progressr::handlers(global = TRUE)

## ---------------------------------------------------------------------
## Function to reshape one inner result (a classed list) into a tidy M-row data frame.
## ---------------------------------------------------------------------
## - Length-M vectors -> columns (one row per accepted sensitivity draw).
## - Scalars          -> recycled to M rows as bootstrap-level columns.
## - Call             -> dropped.
## - boot, draw       -> identifiers added in front.
fit_to_df <- function(fit, boot_id){
  fit$Call <- NULL
  ## the six bootstrap-level scalars in the inner output
  scalar_names <- c("cov_center", "x_eval", "mean.S0", "var.S0", "mean.S1", "var.S1")
  scal <- fit[intersect(scalar_names, names(fit))]
  vecs <- fit[setdiff(names(fit), names(scal))]
  df <- as.data.frame(vecs, stringsAsFactors = FALSE)   # M rows
  for (nm in names(scal)) df[[nm]] <- scal[[nm]]        # length-1 recycles to M
  cbind(boot = boot_id, draw = seq_len(nrow(df)), df)
}


## -----------------------------------------------------------------------------
## Main Function ---------------------------------------------------------------
## ---------------------------------------------------------------------
## Dataset columns are referenced by NAME (character strings) to keep the
## non-standard evaluation out of the parallel workers. Internally the data are
## relabelled to Surr/True/Treat/Cov so the inner function's substitute()-based
## extraction works on the resampled data frame.
run_nested_ICA <- function(Dataset,
                           surr_name, true_name, treat_name, cov_name, # Need to put "" around the name
                           Theta.S_0, Theta.S_1,
                           Nboot = 300,
                           M = 1000,
                           N.MC = 10000,
                           global_seed = 12345,
                           workers = max(1L, parallel::detectCores() - 2L),
                           inner_path = here("Code","ICA.BinCont.Cov_inner.R"), # Path to the function of the inner loop
                           G_pi_10 = c(0, 1),
                           G_rho_01_00 = c(-1, 1), G_rho_01_01 = c(-1, 1),
                           G_rho_01_10 = c(-1, 1), G_rho_01_11 = c(-1, 1),
                           Center.Cov = TRUE, x_eval = 0,
                           Monotonicity = FALSE, Independence = FALSE,
                           HAA = FALSE, Cond_ind = FALSE,
                           progress = TRUE,           # print "Progress: x%: x/Nboot bootstrap done"
                           out_dir = ".",             # directory for the .csv and log files (created if missing)
                           file_prefix = "ICA_draws", # file-name stem; timestamp dd-mm-yy-hh-mm is appended
                           save_csv = TRUE){          # write draws + metadata to .csv (log is always written)
  # 1. Preparation ------------------------------------------------------------
  
  ## make ICA.BinCont.Cov() available (sourced into the global env so the
  ## parallel workers can pick it up as an exported global)
  source(inner_path)

  if (is.null(cov_name)){
    base <- data.frame(Surr  = Dataset[[surr_name]],
                       True  = Dataset[[true_name]],
                       Treat = Dataset[[treat_name]])
  } else {
    ## standardised, complete-case data with fixed column names
    base <- data.frame(Surr  = Dataset[[surr_name]],
                       True  = Dataset[[true_name]],
                       Treat = Dataset[[treat_name]],
                       Cov   = Dataset[[cov_name]])
  }

  base <- base[stats::complete.cases(base), , drop = FALSE]

  ## parallel backend (restored on exit)
  oplan <- plan(multisession, workers = workers)
  on.exit(plan(oplan), add = TRUE)

  ## one global seed -> reproducible, independent per-bootstrap streams
  set.seed(global_seed)

  ## Progress handler producing e.g. "Progress: 30%: 30/100 bootstrap done"
  prog_handler <- progressr::handler_progress(
    format = "Progress: :percent: :current/:total bootstrap done",
    clear  = FALSE, width = 80)

  ## The outer loop, wrapped so progress can be signalled once per bootstrap.
  ## p() advances the bar; it is a no-op when progress = FALSE. Each task returns
  ## list(boot, df, log, ok); warnings and errors are captured PER bootstrap in
  ## the worker (parallel-safe) and returned as strings, so one failure never
  ## stops the run. The master writes all messages to one log file afterwards.
  do_run <- function(){
    p <- if (progress) progressr::progressor(steps = Nboot) else function(...) invisible(NULL)
    foreach(
      b = seq_len(Nboot),
      .options.future = list(seed = TRUE) # independent L'Ecuyer stream per task
    ) %dofuture% {

      ## ---- OUTER: stratified bootstrap (arm sizes fixed), from the worker stream ----
      A  <- base[base$Treat == -1, , drop = FALSE]
      Bx <- base[base$Treat ==  1, , drop = FALSE]
      boot_dat <- rbind(A [sample(nrow(A),  replace = TRUE), , drop = FALSE],
                        Bx[sample(nrow(Bx), replace = TRUE), , drop = FALSE])

      ## message accumulator for this bootstrap (env avoids <<- scoping issues)
      acc <- new.env(parent = emptyenv()); acc$msg <- character(0)
      addmsg <- function(x) acc$msg <- c(acc$msg, x)

      ## ---- INNER: full M-draw sensitivity sweep on this resample ----
      ## Seed = NULL -> use this worker's stream (do not reseed).
      ## withCallingHandlers captures warnings (then muffles them and keeps going);
      ## tryCatch captures a fatal error and returns NULL for this bootstrap only.
      if (is.null(cov_name)){
        fit <- withCallingHandlers(
          tryCatch(
          ICA.BinCont.v2(
              Dataset = boot_dat,
              Surr = Surr,
              True = True,             
              Treat = Treat, 
              G_pi_10 = G_pi_10,
              G_rho_01_00 = G_rho_01_00, G_rho_01_01 = G_rho_01_01,
              G_rho_01_10 = G_rho_01_10, G_rho_01_11 = G_rho_01_11,
              Theta.S_0 = Theta.S_0, Theta.S_1 = Theta.S_1,
              M = M, Seed = NULL,
              N.MC = N.MC,
              Monotonicity = Monotonicity, Independence = Independence,
              HAA = HAA, Cond_ind = Cond_ind,
              Plots = FALSE, Save.Plots = "No", Show.Details = FALSE
            ),
            error = function(e){ addmsg(paste0("ERROR in fitting model without covariate: ", conditionMessage(e))); NULL }),
          warning = function(w){ addmsg(paste0("WARNING in fitting model without covariate: ", conditionMessage(w))); invokeRestart("muffleWarning") })
        
      } else {
        fit <- withCallingHandlers(
          tryCatch(
            ICA.BinCont.Cov(
              Dataset = boot_dat,
              Surr = Surr, True = True, Treat = Treat, Cov = Cov,   # bare names = the columns above
              G_pi_10 = G_pi_10,
              G_rho_01_00 = G_rho_01_00, G_rho_01_01 = G_rho_01_01,
              G_rho_01_10 = G_rho_01_10, G_rho_01_11 = G_rho_01_11,
              Theta.S_0 = Theta.S_0, Theta.S_1 = Theta.S_1,
              M = M, Seed = NULL, N.MC = N.MC,
              Center.Cov = Center.Cov, x_eval = x_eval,
              Monotonicity = Monotonicity, Independence = Independence,
              HAA = HAA, Cond_ind = Cond_ind,
              Plots = FALSE, Save.Plots = "No", Show.Details = FALSE),
            error = function(e){ addmsg(paste0("ERROR in fitting model with covariate: ", conditionMessage(e))); NULL }),
          warning = function(w){ addmsg(paste0("WARNING in fitting model with covariate: ", conditionMessage(w))); invokeRestart("muffleWarning") })
      }
      

      ## signal one bootstrap done (relayed to the master as each task resolves)
      p()

      ## Output to dataframe (df = NULL if this bootstrap failed)
      df <- if (is.null(fit)) NULL else fit_to_df(fit, boot_id = b)
      list(boot = b, df = df, log = acc$msg, ok = !is.null(fit))
    }
  }

  t0 <- Sys.time()
  if (progress) {
    progressr::with_progress(results <- do_run(), handlers = prog_handler)
  } else {
    results <- do_run()
  }
  runtime_min <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 3)

  # 3. Build output paths (prefix + dd-mm-yy-hh-mm) ---------------------------
  ts <- format(Sys.time(), "%y-%m-%d-%H-%M")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_path   <- file.path(out_dir, paste0(file_prefix, "_log_",  ts, ".txt"))
  draws_path <- file.path(out_dir, paste0(file_prefix, "_",      ts, ".csv"))
  meta_path  <- file.path(out_dir, paste0(file_prefix, "_meta_", ts, ".csv"))

  # 4. Collect and WRITE the log (warnings + errors) --------------------------
  log_lines <- unlist(lapply(results, function(r)
    if (length(r$log)) sprintf("[boot %d] %s", r$boot, r$log) else character(0)))
  n_failed <- sum(vapply(results, function(r) !isTRUE(r$ok), logical(1)))
  n_warn   <- length(grep("WARNING:", log_lines, fixed = TRUE))
  writeLines(c(
    sprintf("Run log           : %s", ts),
    sprintf("Nboot=%d  M=%d  N.MC=%d  seed=%d  workers=%d", Nboot, M, N.MC, global_seed, workers),
    sprintf("Runtime           : %s min", runtime_min),
    sprintf("Failed bootstraps : %d / %d", n_failed, Nboot),
    sprintf("Warning messages  : %d", n_warn),
    "----------------------------------------------------------------",
    if (length(log_lines)) log_lines else "No warnings or errors."),
    con = log_path)

  # 5. Stack the successful bootstraps ----------------------------------------
  res_list <- lapply(results, `[[`, "df")
  res_list <- res_list[!vapply(res_list, is.null, logical(1))]
  if (!length(res_list)){
    warning("All bootstraps failed. See log: ", log_path)
    return(invisible(NULL))
  }
  draws <- do.call(rbind, res_list)          # up to Nboot * M rows
  rownames(draws) <- NULL
  attr(draws, "Nboot")       <- Nboot
  attr(draws, "M")           <- M
  attr(draws, "global_seed") <- global_seed
  attr(draws, "log_file")    <- log_path

  # 6. Metadata ---------------------------------------------------------------
  ## helper: collapse a numeric vector (e.g. a range or the Theta vectors) to a string
  fmt <- function(v) paste(v, collapse = ", ")
  cov_name <- ifelse(is.null(cov_name), NA_character_, cov_name)
  
  meta <- data.frame(
    key = c(## --- run bookkeeping ---
            "timestamp","runtime_min","n_boot_returned","n_rows",
            "n_failed_boot","n_warnings",
            ## --- run size / reproducibility ---
            "Nboot","M","N.MC","global_seed","workers",
            ## --- data columns ---
            "surr","true","treat","cov",
            ## --- model configuration ---
            "G_pi_10","G_rho_01_00","G_rho_01_01","G_rho_01_10","G_rho_01_11",
            "Center.Cov","x_eval","Monotonicity","Independence","HAA","Cond_ind",
            "Theta.S_0","Theta.S_1","inner_path"),
    value = c(ts, runtime_min, length(unique(draws$boot)), nrow(draws),
              n_failed, n_warn,
              Nboot, M, N.MC, global_seed, workers,
              surr_name, true_name, treat_name, cov_name,
              fmt(G_pi_10), fmt(G_rho_01_00), fmt(G_rho_01_01), fmt(G_rho_01_10), fmt(G_rho_01_11),
              Center.Cov, x_eval, Monotonicity, Independence, HAA, Cond_ind,
              fmt(Theta.S_0), fmt(Theta.S_1), inner_path),
    stringsAsFactors = FALSE)

  # 7. Write the CSVs ---------------------------------------------------------
  if (isTRUE(save_csv)){
    utils::write.csv(as.data.frame(draws), draws_path, row.names = FALSE)
    utils::write.csv(meta, meta_path, row.names = FALSE)
    attr(draws, "draws_file") <- draws_path
    attr(draws, "meta_file")  <- meta_path
    message("Wrote draws : ", draws_path)
    message("Wrote meta  : ", meta_path)
  }
  message("Wrote log   : ", log_path)

  draws
}


