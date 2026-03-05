# =============================================================================
# main.R — Master Run Script for EACC-Meta Simulation Study
# =============================================================================
# Externally Anchored Covariate Completion for Meta-analysis with
# Incomplete Aggregate Data: A Simulation Study on Transportability
#
# Authors: Lingyao Sun, Mingye Zhao, Dachuang Zhou, Yi Rong,
#          Xiatong Ke, Lei Tian
# Manuscript: JCEPI-D-26-00191 (Journal of Clinical Epidemiology)
#
# Study designed and reported following the ADEMP framework
# (Morris et al., 2019; Siepe et al., 2024)
#
# Usage:
#   1. Edit R/00_config.R to set WORK_DIR and parameters
#   2. Place RWD files in WORK_DIR/real_world_data/
#   3. source("main.R")
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("
============================================================
  EACC-Meta Simulation Study
  ADEMP Framework: Aims | Data-generating | Estimands |
                   Methods | Performance
============================================================
\n")

# --- 0. Install and load packages --------------------------------------------
suppressPackageStartupMessages({
  pkgs <- c(
    "dplyr", "tibble", "readr", "purrr", "stringr", "tidyr",
    "MASS", "survey", "readxl", "data.table", "meta",
    "future", "future.apply", "progressr"
  )
  need <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(need)) install.packages(need, dependencies = TRUE)
  lapply(pkgs, library, character.only = TRUE)
})
options(survey.lonely.psu = "adjust")

# --- 1. Source all modules ----------------------------------------------------
source("R/00_config.R")
source("R/01_utilities.R")
source("R/02_data_loading.R")
source("R/03_data_generation.R")
source("R/04_raking_matching.R")
source("R/05_blup_imputation.R")
source("R/06_wls_extrapolation.R")
source("R/07_meta_pooling.R")
source("R/08_fh_bounds.R")
source("R/09_fullboot_pipeline.R")
source("R/10_run_simulation.R")
source("R/11_postprocessing.R")

# --- 2. Set up output directory -----------------------------------------------
OUT_ROOT_PARENT <- file.path(WORK_DIR, "outputs_meta_reps")
safe_mkdir(OUT_ROOT_PARENT)
RUN_ID  <- format(as.POSIXct(Sys.time(), tz = "Asia/Shanghai"),
                   "%Y%m%d_%H%M%S")
RUN_TAG <- sprintf("RUN_%s__REPS%d__FULLBOOT%d", RUN_ID, N_REP, B_BOOT)
OUT_ROOT_RUN <- file.path(OUT_ROOT_PARENT, RUN_TAG)
safe_mkdir(OUT_ROOT_RUN)
message("[Run folder] ", OUT_ROOT_RUN)

# --- 3. Configure parallel computing -----------------------------------------
future::plan(future::multisession, workers = N_WORKERS)
progressr::handlers(global = TRUE)
progressr::handlers("txtprogressbar")
options(future.globals.maxSize = 8 * 1024^3)
message("[Parallel] workers = ", N_WORKERS)

# --- 4. Load real-world data --------------------------------------------------
rwd_data <- load_rwd_pools(RWD_DIR, MATCH_POOLS, GEN_POOL)
rwd_std     <- rwd_data$rwd_std
MATCH_POOLS <- rwd_data$match_pools
pool_df_gen <- rwd_data$pool_df_gen
message("[Data] Loaded pools: ", paste(names(rwd_std), collapse = ", "))

# --- 5. Run simulation --------------------------------------------------------
N_SCENES        <- length(pL_SET) * length(pV_SET)
N_CASES_PER_REP <- N_SCENES * length(K_SUB)
N_TICKS         <- as.integer(N_REP * N_CASES_PER_REP)

message(sprintf(
  "\n[Simulation] %d reps x %d scenes x %d K-subsets = %d total cases",
  N_REP, N_SCENES, length(K_SUB), N_REP * N_CASES_PER_REP
))
message(sprintf("[Simulation] %d bootstrap draws per case\n", B_BOOT))

REP_ALL <- progressr::with_progress({
  p <- progressr::progressor(steps = N_TICKS)

  future.apply::future_lapply(seq_len(N_REP), function(r) {
    setwd(WORK_DIR)
    safe_mkdir(rep_dir(r))

    # Save seed manifest
    data.table::fwrite(
      data.table::data.table(
        rep = r,
        SEED_GEN = SEED_GEN, SEED_EXT = SEED_EXT, SEED_BOOT = SEED_BOOT
      ),
      file.path(rep_dir(r), "REP_MANIFEST.csv")
    )

    rep_dt <- run_one_rep(
      rep_id = r,
      rwd_std_local = rwd_std,
      pool_df_gen_local = pool_df_gen,
      p = p
    )

    data.table::fwrite(
      rep_dt,
      file.path(OUT_ROOT_RUN, sprintf("REP_LEVEL_rep%03d.csv", r))
    )
    data.table::fwrite(rep_dt, file.path(rep_dir(r), "REP_LEVEL.csv"))
    rep_dt
  }, future.seed = TRUE)
})

# --- 6. Aggregate results -----------------------------------------------------
REP_LEVEL_ALL <- data.table::rbindlist(REP_ALL, use.names = TRUE, fill = TRUE)
data.table::fwrite(REP_LEVEL_ALL,
                    file.path(OUT_ROOT_RUN, "REP_LEVEL_ALL.csv"))

CASE_KEYS <- c("scene_id", "pL", "pV", "K", "pool", "varset",
               "degree", "meta_strategy", "meta_model")

SUMMARY_200 <- compute_summary(REP_LEVEL_ALL, CASE_KEYS)
data.table::setcolorder(
  SUMMARY_200,
  c(CASE_KEYS, "n_rep", "mean_trut", "mean_est", "mean_bias",
    "mse", "rmse", "cover", "mcse_cover", "mean_wid")
)
data.table::fwrite(SUMMARY_200,
                    file.path(OUT_ROOT_RUN, "SUMMARY_200reps.csv"))

message("\n======================================")
message("  SIMULATION COMPLETE")
message("  Rep-level:  ", file.path(OUT_ROOT_RUN, "REP_LEVEL_ALL.csv"))
message("  Summary:    ", file.path(OUT_ROOT_RUN, "SUMMARY_200reps.csv"))
message("======================================\n")

# --- 7. Rubin post-processing -------------------------------------------------
run_rubin_postprocessing(OUT_ROOT_RUN)

message("\n============ ALL DONE ============\n")
