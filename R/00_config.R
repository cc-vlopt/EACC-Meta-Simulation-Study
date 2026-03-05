# =============================================================================
# 00_config.R — User-Configurable Parameters
# =============================================================================
# EACC-Meta Simulation Study
# ADEMP Component: Aims (A) + Data-generating mechanisms (D) — parameter setup
#
# All user-modifiable settings are consolidated here.
# Edit this file to change simulation scenarios; no other files need editing.
# =============================================================================

# --- 1. Random Seed Management ------------------------------------------------
# Three independent RNG streams to avoid cross-contamination
SEED_BASE <- 20260112L
SEED_GEN  <- SEED_BASE + 100000000L   # Data generation
SEED_EXT  <- SEED_BASE + 200000000L   # External matching
SEED_BOOT <- SEED_BASE + 300000000L   # FULLBOOT draws

# --- 2. Simulation Design (ADEMP: Aims) --------------------------------------
N_REP <- 200L                          # Number of simulation repetitions

# Base-case scenario grid: pL (large-trial proportion) x pV (Age65+ proportion)
pL_SET <- c(0.2, 0.5, 0.8)
pV_SET <- c(0.10, 0.30, 0.50)

# Trial universe and subsets
K_MAX <- 15L                           # Maximum number of trials per scenario
K_SUB <- c(5L, 10L, 15L)              # Subsets for meta-analysis

# Transport intensity (standardized shift between trial and target populations)
DEGREE_SET <- c(0.5, 0.8, 1.2)

# --- 3. FULLBOOT Settings (ADEMP: Methods) -----------------------------------
B_BOOT          <- 200L               # Bootstrap draws (use 20 for debugging)
N_LARGE_CUTOFF  <- 2000L              # Threshold for large-trial classification
PERTURB_TE      <- TRUE               # Perturb treatment effects in bootstrap
RANDOM_BLUP     <- TRUE               # Random BLUP imputation (vs. deterministic)
CLIP_RANGE      <- c(0, 1)            # Clipping bounds for imputed proportions
ETA_SIGMA       <- 0.9                # Shrinkage toward base Sigma
JITTER_EIG_MIN  <- 1e-8               # Minimum eigenvalue for PSD enforcement
MAX_TRY_MATCH   <- 2L                 # Maximum raking attempts before fallback

# --- 4. External Data Pools ---------------------------------------------------
MATCH_POOLS <- c("CVD_and_DM", "CVD_all")
GEN_POOL    <- "CVD_and_DM"           # Pool used for trial population generation

# --- 5. Covariate Definitions -------------------------------------------------
COVARS  <- c("age65", "htn", "male", "smoke")
EM_VARS <- c("age65", "htn")          # Effect modifiers

# Variable set configurations (full vs. omit each EM)
VARSETS <- list(
  Full       = c("age65", "htn", "male", "smoke"),
  Omit_age65 = c("htn", "male", "smoke"),
  Omit_htn   = c("age65", "male", "smoke"),
  Omit_male  = c("age65", "htn", "smoke"),
  Omit_smoke = c("age65", "htn", "male")
)

# --- 6. Data-Generating Process Parameters (ADEMP: DGP) ----------------------
tau_AB <- -0.20                        # Average treatment effect (probit scale)
beta   <- c(age65 = 0.35, htn = 0.25, male = 0.26, smoke = 0.30)  # Prognostic
delta  <- c(age65 = -0.15, htn = -0.10)                            # Effect modification
mu_a   <- 0                            # Mean of trial baseline intercept

# --- 7. Meta-Analysis Models --------------------------------------------------
META_MODELS              <- c("FE", "RE_REML")
RE_HAKN_TRY              <- TRUE
RE_METHOD_TAU_PRIMARY    <- "REML"
RE_METHOD_TAU_FALLBACK   <- "DL"

# --- 8. Output Controls -------------------------------------------------------
SAVE_IPD             <- FALSE          # Save individual participant data (large!)
SAVE_SCENE_INPUTS    <- TRUE           # Save universe-level inputs
SAVE_K_INPUTS        <- TRUE           # Save K-level plan/cov/agd9/truth
SAVE_BASE_MATCH      <- TRUE           # Save external matching diagnostics
SAVE_WLS_COEF_DRAWS  <- TRUE           # Save WLS coefficients per draw
SAVE_WLS_PRED_DRAWS  <- TRUE           # Save WLS predictions per draw
SAVE_META_DRAWS      <- TRUE           # Save meta-analysis draw-level results
SAVE_TRIAL_INPUTS    <- TRUE           # Save trial-level inputs per draw
SAVE_IMPUTE_DIAG     <- TRUE           # Save BLUP imputation diagnostics
SAVE_SIGMA_DIAG      <- TRUE           # Save covariance diagnostics per draw

# --- 9. Directory Configuration -----------------------------------------------
# Set WORK_DIR to the project root directory.
# RWD files should be in WORK_DIR/real_world_data/
WORK_DIR <- getwd()
RWD_DIR  <- file.path(WORK_DIR, "real_world_data")

# --- 10. Parallel Computing ---------------------------------------------------
N_WORKERS <- max(1L, parallel::detectCores(logical = TRUE) - 2L)
