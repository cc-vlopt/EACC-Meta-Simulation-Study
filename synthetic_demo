# ============================================================
# Minimal worked example for EACC-Meta vs Meta (single scenario)
# Outputs: all key result tables + required figures in one folder
# Scenario (main analysis): pL=0.5, pV(age65)=0.3, K=10
# Degrees: 0.5 / 0.8 / 1.2 SD
# FULLBOOT draws: 200
# ============================================================

rm(list=ls()); options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  pkgs <- c("dplyr","tibble","readr","purrr","stringr","tidyr",
            "survey","data.table","meta","ggplot2","patchwork","openxlsx")
  need <- pkgs[!pkgs %in% rownames(installed.packages())]
  if(length(need)) install.packages(need, dependencies = TRUE)
  lapply(pkgs, library, character.only = TRUE)
})
options(survey.lonely.psu = "adjust")

`%||%` <- function(a, b) if(!is.null(a)) a else b
safe_mkdir <- function(p) if(!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) User parameters (minimal example)
# -----------------------------
SEED_BASE <- 20260112
SEED_GEN  <- SEED_BASE + 100000000L   # data generation
SEED_EXT  <- SEED_BASE + 200000000L   # external matching base
SEED_BOOT <- SEED_BASE + 300000000L   # FULLBOOT draws

# Main analysis scenario
K        <- 10L
pL       <- 0.5
pV_age65 <- 0.30
DEGREE_SET <- c(0.5, 0.8, 1.2)

# FULLBOOT
B_BOOT <- 200L
N_LARGE_CUTOFF <- 2000L
PERTURB_TE <- TRUE
RANDOM_BLUP <- TRUE
CLIP_RANGE <- c(0, 1)
ETA_SIGMA <- 0.9
JITTER_EIG_MIN <- 1e-8

# External data pool to use in this minimal example
POOL_NAME <- "CVD_and_DM"  # change to "CVD_all" if you prefer

# Variable set (Full only for the minimal example)
VARSET_NAME <- "Full"
VARSET <- c("age65","htn","male","smoke")

# Output root
WORK_DIR <- getwd()
OUT_DIR <- file.path(WORK_DIR, "outputs_minimal_example_EACC_meta")
safe_mkdir(OUT_DIR)
message("[Output folder] ", OUT_DIR)

# If you have RWD pools, put them here:
#   ./真实世界数据/final_analysis_dataset_CVD_and_DM.csv
#   ./真实世界数据/final_analysis_dataset_CVD_all.csv
RWD_DIR <- file.path(WORK_DIR, "真实世界数据")

# -----------------------------
# 2) Utilities and standardisation
# -----------------------------
round_half_up <- function(x) floor(x + 0.5)

rtruncnorm_1 <- function(mu, sigma, a, b){
  pa <- pnorm((a - mu)/sigma)
  pb <- pnorm((b - mu)/sigma)
  u  <- runif(1, pa, pb)
  mu + sigma * qnorm(u)
}

col_has <- function(df, patt){
  nm <- names(df); hit <- nm[stringr::str_detect(tolower(nm), patt)]
  if(length(hit)) hit[1] else NA_character_
}

std_bin <- function(x){
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) return(as.integer(tolower(x) %in% c("1","y","yes","true","男","male","m","current","ever")))
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    if (all(na.omit(x) %in% c(0,1))) return(as.integer(x))
    if (is.finite(min(x,na.rm=TRUE)) && is.finite(max(x,na.rm=TRUE)) &&
        min(x,na.rm=TRUE)>=0 && max(x,na.rm=TRUE)<=1) return(as.integer(x>=0.5))
    if (all(na.omit(x) %in% c(1,2))) return(as.integer(x==1))
    return(as.integer(x==1))
  }
  stop("Cannot standardize binary variable")
}

prep_covars_4 <- function(df){
  nm_age65 <- col_has(df, "age\\s*65|age65|65\\+|65plus")
  nm_age   <- col_has(df, "^age$|age_year|agey")
  age65 <- if(!is.na(nm_age65)) std_bin(df[[nm_age65]]) else if(!is.na(nm_age) && is.numeric(df[[nm_age]])) as.integer(df[[nm_age]]>=65) else rbinom(nrow(df),1,0.5)

  nm_htn <- col_has(df,"htn|hypertens")
  htn <- if(!is.na(nm_htn)) std_bin(df[[nm_htn]]) else rbinom(nrow(df),1,0.5)

  nm_sex <- col_has(df,"^male$|sex|gender")
  male <- if(!is.na(nm_sex)){
    v <- df[[nm_sex]]
    if(is.character(v)||is.factor(v)) as.integer(tolower(v) %in% c("male","m","1","男")) else std_bin(v)
  } else rbinom(nrow(df),1,0.5)

  nm_smok <- col_has(df,"smok")
  smoke <- if(!is.na(nm_smok)) std_bin(df[[nm_smok]]) else rbinom(nrow(df),1,0.5)

  tibble::tibble(
    age65 = as.integer(age65),
    htn   = as.integer(htn),
    male  = as.integer(male),
    smoke = as.integer(smoke)
  )
}

read_any <- function(base){
  cand <- file.path(RWD_DIR, paste0(base, c(".csv",".xlsx",".xls")))
  if(file.exists(cand[1])) return(readr::read_csv(cand[1], show_col_types = FALSE))
  if(file.exists(cand[2])) return(readxl::read_xlsx(cand[2]))
  if(file.exists(cand[3])) return(readxl::read_xls(cand[3]))
  stop("Cannot find: ", base, " (.csv/.xlsx/.xls) in ", RWD_DIR)
}

# -----------------------------
# 3) Load RWD pools (or generate a synthetic pool)
# -----------------------------
rwd_map <- list()
if (dir.exists(RWD_DIR)) {
  if (file.exists(file.path(RWD_DIR, "final_analysis_dataset_CVD_all.csv")) ||
      file.exists(file.path(RWD_DIR, "final_analysis_dataset_CVD_all.xlsx")) ||
      file.exists(file.path(RWD_DIR, "final_analysis_dataset_CVD_all.xls"))) {
    rwd_map[["CVD_all"]] <- read_any("final_analysis_dataset_CVD_all")
  }
  if (file.exists(file.path(RWD_DIR, "final_analysis_dataset_CVD_and_DM.csv")) ||
      file.exists(file.path(RWD_DIR, "final_analysis_dataset_CVD_and_DM.xlsx")) ||
      file.exists(file.path(RWD_DIR, "final_analysis_dataset_CVD_and_DM.xls"))) {
    rwd_map[["CVD_and_DM"]] <- read_any("final_analysis_dataset_CVD_and_DM")
  }
}

if (!length(rwd_map)) {
  message("[Info] No RWD found under ./真实世界数据. Using a synthetic pool instead.")
  set.seed(12345)
  n_pool <- 200000
  # correlated latent normals -> binary
  Z <- matrix(rnorm(n_pool*4), ncol=4)
  # induce mild correlations
  Z[,2] <- 0.35*Z[,1] + sqrt(1-0.35^2)*Z[,2]
  Z[,3] <- 0.20*Z[,1] + 0.25*Z[,2] + sqrt(1-0.20^2-0.25^2)*Z[,3]
  Z[,4] <- 0.10*Z[,1] + 0.15*Z[,2] + 0.20*Z[,3] + sqrt(1-0.10^2-0.15^2-0.20^2)*Z[,4]
  rwd_map[["CVD_all"]] <- tibble::tibble(
    age65 = as.integer(Z[,1] > qnorm(0.45)),
    htn   = as.integer(Z[,2] > qnorm(0.50)),
    male  = as.integer(Z[,3] > qnorm(0.52)),
    smoke = as.integer(Z[,4] > qnorm(0.48))
  )
  rwd_map[["CVD_and_DM"]] <- rwd_map[["CVD_all"]]
}

rwd_std <- lapply(rwd_map, prep_covars_4)

if(!POOL_NAME %in% names(rwd_std)) stop("POOL_NAME not available: ", POOL_NAME)
pool_df_gen <- rwd_std[[POOL_NAME]]  # also use for generation in this minimal example

COVARS  <- c("age65","htn","male","smoke")
EM_VARS <- c("age65","htn")

# -----------------------------
# 4) DGP (AB, probit)
# -----------------------------
tau_AB <- -0.20
beta   <- c(age65=0.35, htn=0.25, male=0.26, smoke=0.30)
delta  <- c(age65=-0.15, htn=-0.10)

inv_probit <- function(z) pnorm(z)

draw_alpha_k <- function(is_large){
  mu_a <- 0
  sigma <- if (isTRUE(is_large)) 0.10 else 0.05
  a <- mu_a - 2.5*sigma
  b <- mu_a + 2.5*sigma
  rtruncnorm_1(mu = mu_a, sigma = sigma, a = a, b = b)
}

te_probit_from_counts <- function(r1, n1, r2, n2, cc = 0.5, eps = 1e-6){
  # returns te = z1 - z2 where z=Phi^{-1}(p), and delta-method se
  if (any(!is.finite(c(r1,n1,r2,n2))) || any(c(n1,n2) <= 0)) {
    return(list(te=NA_real_, se=NA_real_))
  }
  p1 <- (r1 + cc) / (n1 + 2*cc); p2 <- (r2 + cc) / (n2 + 2*cc)
  p1 <- pmin(pmax(p1, eps), 1 - eps); p2 <- pmin(pmax(p2, eps), 1 - eps)
  z1 <- qnorm(p1); z2 <- qnorm(p2); te <- z1 - z2
  phi1 <- dnorm(z1); phi2 <- dnorm(z2)
  var1 <- p1 * (1 - p1) / (n1 * phi1^2)
  var2 <- p2 * (1 - p2) / (n2 * phi2^2)
  se <- sqrt(var1 + var2)
  list(te = as.numeric(te), se = as.numeric(se))
}

# -----------------------------
# 5) Raking sampling with tolerance ladder
# -----------------------------
# Ladder as described in eMethods:
# (0.05,0.10) -> (0.07,0.12) -> (0.10,0.15) for age65 vs others
TOL_LADDER <- list(
  c(age65=0.05, htn=0.10, male=0.10, smoke=0.10),
  c(age65=0.07, htn=0.12, male=0.12, smoke=0.12),
  c(age65=0.10, htn=0.15, male=0.15, smoke=0.15)
)

sample_with_raking_4 <- function(pool_df, N, targets,
                                 max_try = 10,
                                 tol = c(age65=0.05, htn=0.10, male=0.10, smoke=0.10),
                                 ipf_maxit = 300, ipf_tol = 1e-7) {
  vars <- COVARS
  stopifnot(all(vars %in% names(pool_df)))

  df <- pool_df[, vars, drop = FALSE]
  for (v in vars) df[[v]] <- as.integer(df[[v]] %in% c(1L, TRUE))
  df <- stats::na.omit(df)
  if (nrow(df) < 100) stop("pool_df too small after omitting missing")

  base <- df; base$.id <- seq_len(nrow(base))
  for (v in vars) base[[v]] <- factor(base[[v]], levels = c(0,1))
  dsgn <- survey::svydesign(ids = ~.id, data = base, weights = ~1)
  Npop <- nrow(base)

  pop_margins <- lapply(vars, function(v){
    data.frame(tmp = factor(0:1), Freq = c(1 - targets[[v]], targets[[v]]) * Npop) %>%
      dplyr::rename(!!v := tmp)
  })

  adj <- try(survey::rake(
    dsgn,
    sample.margins = lapply(vars, function(v) as.formula(paste0("~", v))),
    population.margins = pop_margins,
    control = list(maxit = ipf_maxit, epsilon = ipf_tol)
  ), silent = TRUE)

  use_uniform <- inherits(adj, "try-error")
  if (use_uniform) {
    w <- rep(1, nrow(df))
  } else {
    w <- as.numeric(weights(adj))
    w[!is.finite(w) | w <= 0] <- min(w[w > 0], na.rm = TRUE)
    # trim 1% and 99%
    qq <- stats::quantile(w, c(0.01, 0.99), na.rm = TRUE)
    w <- pmin(pmax(w, qq[1]), qq[2])
  }

  best <- NULL; best_dev <- Inf
  for (t in seq_len(max_try)) {
    id <- sample.int(nrow(df), size = N, replace = TRUE, prob = w)
    samp <- df[id, , drop = FALSE]
    obs  <- colMeans(samp[, vars, drop = FALSE])
    devs <- abs(obs - as.numeric(targets[vars]))
    max_dev <- max(devs)
    if (max_dev < best_dev) {
      best_dev <- max_dev
      best <- list(samp = samp, obs = obs, devs = devs)
    }
    if (all(devs <= as.numeric(tol[vars]))) return(best$samp)
  }
  best$samp
}

sample_with_raking_ladder <- function(pool_df, N, targets, max_try=12){
  for(step in seq_along(TOL_LADDER)){
    tol <- TOL_LADDER[[step]]
    samp <- sample_with_raking_4(pool_df=pool_df, N=N, targets=targets, max_try=max_try, tol=tol)
    obs <- colMeans(samp[, COVARS, drop=FALSE])
    devs <- abs(obs - as.numeric(targets[COVARS]))
    if(all(devs <= as.numeric(tol[COVARS]))){
      return(list(sample=samp, ok=TRUE, tol=tol, obs=obs, devs=devs, step=step))
    }
  }
  # return best-effort at the last step
  tol <- TOL_LADDER[[length(TOL_LADDER)]]
  samp <- sample_with_raking_4(pool_df=pool_df, N=N, targets=targets, max_try=max_try, tol=tol)
  obs <- colMeans(samp[, COVARS, drop=FALSE])
  devs <- abs(obs - as.numeric(targets[COVARS]))
  list(sample=samp, ok=FALSE, tol=tol, obs=obs, devs=devs, step=length(TOL_LADDER))
}

# -----------------------------
# 6) Trial simulation (fast)
# -----------------------------
gen_one_trial_AB_fast <- function(trial_id, is_large, N_total, alpha_k, targets, pool_df, sparse_small=TRUE){
  N_A <- floor(N_total/2)
  N_B <- N_total - N_A

  rr <- sample_with_raking_ladder(pool_df = pool_df, N = N_total, targets = targets, max_try = 20)
  samp <- rr$sample
  samp <- as.data.table(samp)

  # randomise A/B
  perm <- sample.int(N_total, N_total, replace = FALSE)
  samp <- samp[perm]
  samp[, Trt := ifelse(.I <= N_A, "A", "B")]
  samp[, t := ifelse(Trt=="B", 1L, 0L)]

  # generate y
  X <- as.matrix(samp[, ..COVARS])
  lp_main <- as.numeric(alpha_k + X %*% beta[COVARS])
  lp_trt  <- samp$t * (tau_AB + as.numeric(as.matrix(samp[, ..EM_VARS]) %*% delta[EM_VARS]))
  lp <- lp_main + lp_trt
  p  <- inv_probit(lp)
  y  <- rbinom(N_total, 1, p)
  samp[, y := y]
  samp[, trial := trial_id]

  # observed covariate means
  cov_obs <- samp[, lapply(.SD, mean), .SDcols = COVARS][1]
  cov_obs[, `:=`(trial = trial_id, N = N_total, is_large = is_large)]
  setcolorder(cov_obs, c("trial","N","is_large",COVARS))

  # arm counts
  arm_counts <- samp[, .(n = .N, events = sum(y), risk = mean(y)), by = .(trial, Trt)]

  # helper: 2-arm counts
  agg_two_arms_dt <- function(dt){
    arm <- dt[, .(r=sum(y), n=.N), by=.(Trt)]
    if(nrow(arm) < 2) return(NULL)
    list(
      rB = arm[Trt=="B"]$r, nB = arm[Trt=="B"]$n,
      rA = arm[Trt=="A"]$r, nA = arm[Trt=="A"]$n
    )
  }

  make_one_row <- function(type, dt_sub, x_age65, x_htn, x_male, x_smoke){
    cnt <- agg_two_arms_dt(dt_sub)
    if(is.null(cnt)) {
      te <- list(te=NA_real_, se=NA_real_)
      cnt <- list(rB=NA_real_, nB=NA_integer_, rA=NA_real_, nA=NA_integer_)
    } else {
      te <- te_probit_from_counts(cnt$rB, cnt$nB, cnt$rA, cnt$nA)
    }
    data.table(
      trial = trial_id, type = type,
      x_age65 = x_age65, x_htn = x_htn, x_male = x_male, x_smoke = x_smoke,
      TE = te$te, se = te$se,
      rB = cnt$rB, nB = cnt$nB, rA = cnt$rA, nA = cnt$nA
    )
  }

  agd_list <- list()
  agd_list[[1]] <- make_one_row(
    "overall", samp,
    mean(samp$age65), mean(samp$htn), mean(samp$male), mean(samp$smoke)
  )

  add_cov_rows <- function(v){
    for(val in c(1,0)){
      dt_sub <- samp[get(v)==val]
      x_age65 <- if(v=="age65") val else NA_real_
      x_htn   <- if(v=="htn")   val else NA_real_
      x_male  <- if(v=="male")  val else NA_real_
      x_smoke <- if(v=="smoke") val else NA_real_
      agd_list[[length(agd_list)+1]] <<- make_one_row(
        paste0(v,"=",val), dt_sub, x_age65, x_htn, x_male, x_smoke
      )
    }
  }
  add_cov_rows("age65"); add_cov_rows("htn"); add_cov_rows("male"); add_cov_rows("smoke")
  agd9 <- rbindlist(agd_list, use.names = TRUE, fill = TRUE)

  # emulate "small trials only report overall": blank out subgroup TE/se
  if(isTRUE(sparse_small) && !isTRUE(is_large)){
    agd9[type!="overall", `:=`(TE=NA_real_, se=NA_real_, rB=NA_real_, nB=NA_integer_, rA=NA_real_, nA=NA_integer_)]
  }

  list(cov_obs=cov_obs, arm_counts=arm_counts, agd9=agd9, raking_ok = rr$ok, raking_devs = rr$devs, raking_step = rr$step)
}

# -----------------------------
# 7) Target generation and truth
# -----------------------------
make_target_probs_equal <- function(base_p, delta_sd, clip=c(0.02,0.98)){
  if (is.data.frame(base_p)) {
    p0 <- as.numeric(unlist(base_p[1, COVARS], use.names = FALSE))
    names(p0) <- COVARS
  } else {
    p0 <- as.numeric(base_p[COVARS])
    names(p0) <- COVARS
  }
  p0 <- pmin(pmax(p0, 1e-6), 1-1e-6)
  z0 <- qnorm(p0)
  # common direction (1,1,1,1) normalized
  d  <- rep(1, length(COVARS)); d  <- d / sqrt(sum(d^2))
  z_star <- z0 + delta_sd * d
  p_star <- pnorm(z_star)
  p_star <- pmin(pmax(p_star, clip[1]), clip[2])
  stats::setNames(p_star, COVARS)
}

truth_TE_target <- function(p_star){
  tau_AB + sum(delta[EM_VARS] * p_star[EM_VARS])
}

# -----------------------------
# 8) External matching and PSD repair
# -----------------------------
choose_N_match <- function(N_total, mode="equal_trial_n", factor=5, min_n=200L, max_n=200000L){
  N_total <- as.integer(N_total)
  if(!is.finite(N_total) || N_total <= 0) return(min_n)
  if(mode == "equal_trial_n") return(N_total)
  as.integer(max(min_n, min(max_n, round(N_total * factor))))
}

make_psd <- function(S, eig_min = JITTER_EIG_MIN){
  S <- as.matrix(S)
  dn <- dimnames(S)
  S <- (S + t(S)) / 2
  S[!is.finite(S)] <- 0
  d <- diag(S)
  d[!is.finite(d)] <- 0
  d <- pmax(d, eig_min)
  diag(S) <- d
  eg <- eigen(S, symmetric = TRUE)
  vals <- eg$values
  vals[!is.finite(vals)] <- eig_min
  vals <- pmax(vals, eig_min)
  S2 <- eg$vectors %*% diag(vals, length(vals)) %*% t(eg$vectors)
  S2 <- (S2 + t(S2)) / 2
  dimnames(S2) <- dn
  S2
}

fh_corr_bounds <- function(p1, p2){
  v1 <- p1*(1-p1); v2 <- p2*(1-p2)
  if(v1 <= 0 || v2 <= 0) return(c(lower=NA_real_, upper=NA_real_))
  p11_lo <- max(0, p1 + p2 - 1)
  p11_hi <- min(p1, p2)
  cov_lo <- p11_lo - p1*p2
  cov_hi <- p11_hi - p1*p2
  s <- sqrt(v1*v2)
  c(lower = cov_lo/s, upper = cov_hi/s)
}

truncate_corr_to_fh <- function(mu_vec, Sigma){
  # Sigma is covariance of (0/1) variables; enforce FH bounds on implied correlations
  mu <- as.numeric(mu_vec); names(mu) <- COVARS
  S <- as.matrix(Sigma)
  vars <- diag(S)
  sdv <- sqrt(pmax(vars, 1e-12))
  R <- S / (sdv %o% sdv)
  diag(R) <- 1
  n_fix <- 0L
  for(i in 1:(length(COVARS)-1)){
    for(j in (i+1):length(COVARS)){
      b <- fh_corr_bounds(mu[i], mu[j])
      if(is.na(b["lower"]) || is.na(b["upper"])) next
      if(!is.finite(R[i,j])) next
      r0 <- R[i,j]
      r1 <- min(max(r0, b["lower"]), b["upper"])
      if(abs(r1 - r0) > 1e-12){
        R[i,j] <- R[j,i] <- r1
        n_fix <- n_fix + 1L
      }
    }
  }
  # reconstruct covariance
  S2 <- (sdv %o% sdv) * R
  rownames(S2) <- colnames(S2) <- COVARS
  list(Sigma = make_psd(S2), n_fix = n_fix)
}

match_external_base_4 <- function(pool_df, target_p, N_match, seed=NULL){
  if(!is.null(seed)) set.seed(seed)
  rr <- sample_with_raking_ladder(pool_df = pool_df, N = as.integer(N_match), targets = target_p, max_try = 25)
  samp <- rr$sample
  fallback <- !isTRUE(rr$ok)
  samp_use <- samp[, COVARS, drop=FALSE]
  mu0 <- colMeans(samp_use)
  S0  <- stats::cov(as.matrix(samp_use))
  S0  <- (S0 + t(S0))/2
  rownames(S0) <- colnames(S0) <- COVARS
  # FH truncation + PSD
  tr <- truncate_corr_to_fh(mu0, S0)
  list(sample=samp_use, mu=mu0, Sigma=tr$Sigma, fallback=fallback,
       match_ok=rr$ok, match_step=rr$step, match_devs=rr$devs,
       fh_fix_pairs=tr$n_fix)
}

# -----------------------------
# 9) BLUP completion (draw-wise)
# -----------------------------
# Uses MVN conditional mean/cov on covariate means. Random draws optional.
# Adds light observation variance on observed components: Var(p_hat) = p(1-p)/n,
# where n is taken from nA+nB in each row (if available). This follows the eMethods description.
blup_impute_agd9_draw <- function(agd_trial, Sigma, mu_vec, random=TRUE, clip=c(0,1)){
  xcols <- paste0("x_", COVARS)
  if(!all(xcols %in% names(agd_trial))) stop("agd_trial missing x_ cols")
  if(!"type" %in% names(agd_trial)) stop("agd_trial missing type col")

  Sigma <- as.matrix(Sigma)
  rownames(Sigma) <- colnames(Sigma) <- COVARS
  Sigma <- make_psd(Sigma)

  mu_vec <- as.numeric(mu_vec); names(mu_vec) <- COVARS
  if(any(!is.finite(mu_vec))) stop("mu_vec has NA/Inf")

  out <- agd_trial
  for(v in COVARS) out[[paste0("x_",v,"_raw")]] <- out[[paste0("x_",v)]]

  Xraw <- as.matrix(out[, xcols, drop=FALSE])
  colnames(Xraw) <- COVARS

  imp_n <- setNames(rep(0L, length(COVARS)), COVARS)
  clip_n <- setNames(rep(0L, length(COVARS)), COVARS)

  # row-wise sample size for observed variance
  n_row <- rep(NA_real_, nrow(out))
  if(all(c("nA","nB") %in% names(out))){
    n_row <- as.numeric(out$nA) + as.numeric(out$nB)
  }

  do_one_row <- function(x_raw, n_eff){
    obs <- which(is.finite(x_raw))
    mis <- which(!is.finite(x_raw))
    x_imp <- x_raw

    if(length(mis)==0) return(list(x_final=x_imp, n_imp=rep(0L,4), n_clip=rep(0L,4)))

    # observation variance on observed components
    # only applied when obs are proportions in (0,1); fixed 0/1 vars give var 0
    R_obs <- matrix(0, nrow=length(obs), ncol=length(obs))
    if(is.finite(n_eff) && n_eff > 0 && length(obs)>0){
      p_obs <- pmin(pmax(x_raw[obs], 0), 1)
      v_obs <- p_obs*(1-p_obs)/n_eff
      v_obs[!is.finite(v_obs)] <- 0
      R_obs <- diag(v_obs, nrow=length(obs))
    }

    if(length(obs)==0){
      mu_m <- mu_vec[mis]
      if(random){
        S_mm <- make_psd(Sigma[mis,mis,drop=FALSE])
        mu_m <- as.numeric(MASS::mvrnorm(1, mu=mu_m, Sigma=S_mm))
      }
      mu0 <- mu_m
      mu_m <- pmin(pmax(mu_m, clip[1]), clip[2])
      n_clip <- rep(0L,4); n_clip[mis] <- as.integer(mu_m != mu0)
      n_imp  <- rep(0L,4); n_imp[mis] <- 1L
      x_imp[mis] <- mu_m
      return(list(x_final=x_imp, n_imp=n_imp, n_clip=n_clip))
    }

    mu_o <- mu_vec[obs]; mu_m <- mu_vec[mis]
    x_o  <- x_raw[obs]

    S_oo <- Sigma[obs,obs,drop=FALSE] + R_obs
    S_mo <- Sigma[mis,obs,drop=FALSE]
    S_mm <- Sigma[mis,mis,drop=FALSE]
    inv_Soo <- tryCatch(solve(S_oo), error=function(e) MASS::ginv(S_oo))

    mu_cond <- as.numeric(mu_m + S_mo %*% inv_Soo %*% (x_o - mu_o))

    if(random){
      S_cond <- S_mm - S_mo %*% inv_Soo %*% t(S_mo)
      S_cond <- make_psd(S_cond)
      x_m <- as.numeric(MASS::mvrnorm(1, mu=mu_cond, Sigma=S_cond))
    } else {
      x_m <- mu_cond
    }

    x0 <- x_m
    x_m <- pmin(pmax(x_m, clip[1]), clip[2])

    n_clip <- rep(0L,4); n_clip[mis] <- as.integer(x_m != x0)
    n_imp  <- rep(0L,4); n_imp[mis] <- 1L

    x_imp[mis] <- x_m
    list(x_final=x_imp, n_imp=n_imp, n_clip=n_clip)
  }

  imp_list <- lapply(seq_len(nrow(Xraw)), function(i) do_one_row(Xraw[i,], n_row[i]))
  Xfinal <- do.call(rbind, lapply(imp_list, `[[`, "x_final"))
  n_imp  <- do.call(rbind, lapply(imp_list, `[[`, "n_imp"))
  n_clip <- do.call(rbind, lapply(imp_list, `[[`, "n_clip"))
  colnames(Xfinal) <- COVARS

  for(j in seq_along(COVARS)){
    out[[paste0("x_",COVARS[j],"_final")]] <- Xfinal[,j]
    imp_n[COVARS[j]] <- imp_n[COVARS[j]] + sum(n_imp[,j], na.rm=TRUE)
    clip_n[COVARS[j]]<- clip_n[COVARS[j]]+ sum(n_clip[,j],na.rm=TRUE)
  }
  for(v in COVARS) out[[paste0("x_",v)]] <- out[[paste0("x_",v,"_final")]]

  list(agd_imp=out, imp_n=imp_n, clip_n=clip_n)
}

# -----------------------------
# 10) WLS (trial-wise) and prediction
# -----------------------------
fit_wls_trial_draw <- function(agd_imp, varset, perturb_te=FALSE){
  dat <- agd_imp %>%
    dplyr::mutate(
      TE_use = dplyr::if_else(isTRUE(perturb_te) & is.finite(TE_draw), TE_draw, TE),
      w = 1/(se^2),
      w = dplyr::if_else(is.finite(w) & w>0, w, NA_real_)
    ) %>%
    dplyr::filter(is.finite(TE_use), is.finite(se), is.finite(w))

  if(nrow(dat) < 2) return(list(ok=FALSE))

  rhs <- if(length(varset)==0) "1" else paste(paste0("x_",varset,"_final"), collapse=" + ")
  fml <- as.formula(paste0("TE_use ~ ", rhs))
  fit <- tryCatch(lm(fml, data=dat, weights=w), error=function(e) NULL)
  if(is.null(fit)) return(list(ok=FALSE))

  V <- tryCatch(vcov(fit), error=function(e) NULL)
  if(is.null(V)) return(list(ok=FALSE))

  list(ok=TRUE, fit=fit, V=V)
}

predict_wls_at_target <- function(fit_obj, target_vec_named){
  if(!isTRUE(fit_obj$ok)) return(list(te=NA_real_, se=NA_real_))
  fit <- fit_obj$fit; V <- fit_obj$V
  terms <- names(coef(fit))
  xrow <- rep(0, length(terms)); names(xrow) <- terms
  if("(Intercept)" %in% terms) xrow["(Intercept)"] <- 1
  for(v in COVARS){
    nm <- paste0("x_",v,"_final")
    if(nm %in% terms) xrow[nm] <- as.numeric(target_vec_named[[v]])
  }
  te_hat <- as.numeric(sum(xrow * coef(fit)))
  se_hat <- sqrt(as.numeric(t(xrow) %*% V %*% xrow))
  list(te=te_hat, se=se_hat)
}

# -----------------------------
# 11) Meta pooling
# -----------------------------
fe_pool_closed <- function(te, se){
  ok <- is.finite(te) & is.finite(se) & se > 0
  te <- te[ok]; se <- se[ok]
  k_used <- length(te)
  if(k_used < 2){
    return(list(ok=FALSE, k_used=k_used,
                TE=NA_real_, lo=NA_real_, hi=NA_real_,
                tau2=NA_real_, I2=NA_real_, Q=NA_real_))
  }
  w <- 1/(se^2)
  TE_hat <- sum(w*te)/sum(w)
  se_hat <- sqrt(1/sum(w))
  z <- qnorm(0.975)
  lo <- TE_hat - z*se_hat
  hi <- TE_hat + z*se_hat
  Q <- sum(w * (te - TE_hat)^2)
  df <- k_used - 1
  I2 <- ifelse(is.finite(Q) && Q>0, max(0, (Q-df)/Q), 0)
  list(ok=TRUE, k_used=k_used, TE=TE_hat, lo=lo, hi=hi, tau2=0, I2=I2, Q=Q)
}

re_pool_metagen_safe <- function(te, se, method_tau_primary="REML", method_tau_fallback="DL", hakn_try=TRUE){
  ok <- is.finite(te) & is.finite(se) & se > 0
  te <- te[ok]; se <- se[ok]
  k_used <- length(te)
  if(k_used < 2){
    return(list(ok=FALSE, k_used=k_used,
                TE=NA_real_, lo=NA_real_, hi=NA_real_,
                tau2=NA_real_, I2=NA_real_, Q=NA_real_))
  }
  run_one <- function(method_tau, hakn){
    tryCatch(
      meta::metagen(TE=te, seTE=se, sm="MD", comb.fixed=FALSE, comb.random=TRUE,
                    method.tau=method_tau, hakn=hakn),
      error=function(e) NULL
    )
  }
  m <- NULL
  if(isTRUE(hakn_try)){
    m <- run_one(method_tau_primary, TRUE)
    if(is.null(m)) m <- run_one(method_tau_primary, FALSE)
  } else {
    m <- run_one(method_tau_primary, FALSE)
  }
  if(is.null(m)) {
    m <- run_one(method_tau_fallback, FALSE)
  }
  if(is.null(m)) return(list(ok=FALSE, k_used=k_used, TE=NA_real_, lo=NA_real_, hi=NA_real_, tau2=NA_real_, I2=NA_real_, Q=NA_real_))
  list(ok=TRUE, k_used=k_used,
       TE = as.numeric(m$TE.random),
       lo = as.numeric(m$lower.random),
       hi = as.numeric(m$upper.random),
       tau2 = as.numeric(m$tau^2),
       I2 = as.numeric(m$I2)/100,
       Q = as.numeric(m$Q))
}

# -----------------------------
# 12) Build one minimal scenario and run FULLBOOT
# -----------------------------
set.seed(SEED_GEN)

# Trial size plan
nL <- round_half_up(K * pL)
nS <- K - nL
N_large <- sample(8000:16000, nL, replace=TRUE)
N_small <- sample(200:300, nS, replace=TRUE)
N_total <- c(N_large, N_small)

plan <- tibble::tibble(
  trial = sprintf("T%02d", seq_len(K)),
  N_total = N_total,
  is_large = N_total >= N_LARGE_CUTOFF
)

# Trial covariate targets
targets <- plan %>%
  rowwise() %>%
  mutate(
    age65 = pV_age65,
    htn = runif(1, 0.4, 0.6),
    male = runif(1, 0.4, 0.6),
    smoke = runif(1, 0.4, 0.6)
  ) %>% ungroup()

# Trial baselines
alpha_tbl <- plan %>%
  mutate(alpha_k = purrr::map_dbl(is_large, draw_alpha_k))

# Generate trials
trial_out <- vector("list", K)
cov_obs_k <- list()
agd9_k <- list()
arm_counts_k <- list()

for(i in seq_len(K)){
  tr <- plan$trial[i]
  isL <- plan$is_large[i]
  Ni <- plan$N_total[i]
  al <- alpha_tbl$alpha_k[i]
  tg <- as.list(targets[i, COVARS])
  trial_out[[i]] <- gen_one_trial_AB_fast(
    trial_id=tr, is_large=isL, N_total=Ni, alpha_k=al, targets=tg, pool_df=pool_df_gen, sparse_small=TRUE
  )
  cov_obs_k[[i]] <- trial_out[[i]]$cov_obs
  agd9_k[[i]] <- trial_out[[i]]$agd9
  arm_counts_k[[i]] <- trial_out[[i]]$arm_counts
}
cov_obs_k <- dplyr::bind_rows(cov_obs_k)
agd9_k <- data.table::rbindlist(agd9_k, use.names=TRUE, fill=TRUE)
arm_counts_k <- dplyr::bind_rows(arm_counts_k)

# Build target populations and truth
ref_p <- cov_obs_k %>%
  mutate(N_w = N) %>%
  summarise(across(all_of(COVARS), ~ weighted.mean(.x, w=N_w, na.rm=TRUE)))

targets_tbl <- purrr::map_dfr(DEGREE_SET, function(deg){
  p_star <- make_target_probs_equal(ref_p, deg)
  tibble::tibble(degree=deg, !!!as.list(p_star), TE_true = truth_TE_target(p_star))
})

# Save inputs
readr::write_csv(plan, file.path(OUT_DIR, "trial_plan.csv"))
readr::write_csv(targets, file.path(OUT_DIR, "trial_targets.csv"))
readr::write_csv(alpha_tbl, file.path(OUT_DIR, "trial_alpha.csv"))
readr::write_csv(cov_obs_k, file.path(OUT_DIR, "cov_obs_generated.csv"))
readr::write_csv(arm_counts_k, file.path(OUT_DIR, "arm_counts_generated.csv"))
readr::write_csv(as.data.frame(agd9_k), file.path(OUT_DIR, "agd9_generated.csv"))
readr::write_csv(targets_tbl, file.path(OUT_DIR, "targets_and_truth.csv"))

# -----------------------------
# 13) FULLBOOT: precompute base external mu/Sigma for large trials
# -----------------------------
set.seed(SEED_EXT)
pool_df_ext <- rwd_std[[POOL_NAME]]
base_match <- list()
match_diag <- list()

for(i in seq_len(K)){
  tr <- plan$trial[i]
  if(!isTRUE(plan$is_large[i])) next
  N_match <- choose_N_match(plan$N_total[i], mode="equal_trial_n")
  target_p <- as.numeric(cov_obs_k %>% filter(trial==tr) %>% dplyr::select(dplyr::all_of(COVARS)) %>% slice(1))
  names(target_p) <- COVARS
  bm <- match_external_base_4(pool_df_ext, target_p=target_p, N_match=N_match, seed=SEED_EXT + i*1000L)
  base_match[[tr]] <- bm

  ess <- {
    # approximate ESS for equal-weight resample; for diagnostics use sample size
    nrow(bm$sample)
  }
  match_diag[[tr]] <- tibble::tibble(
    trial=tr,
    N_match=nrow(bm$sample),
    fallback=bm$fallback,
    match_ok=bm$match_ok,
    match_step=bm$match_step,
    max_abs_dev=max(abs(bm$match_devs)),
    fh_fix_pairs=bm$fh_fix_pairs,
    mu_age65=bm$mu["age65"], mu_htn=bm$mu["htn"], mu_male=bm$mu["male"], mu_smoke=bm$mu["smoke"]
  )
}
match_diag <- dplyr::bind_rows(match_diag)
readr::write_csv(match_diag, file.path(OUT_DIR, "match_diagnostics_base.csv"))

# -----------------------------
# 14) FULLBOOT draws
# -----------------------------
set.seed(SEED_BOOT)

# storage
meta_draws <- list()
trial_inputs_draws <- list()
wls_coef_draws <- list()
impute_diag <- list()

# helper: fetch agd9 for a trial
get_agd9_trial <- function(tr){
  as.data.frame(agd9_k[agd9_k$trial==tr, ])
}

# helper: overall TE/se row for a trial
get_overall_te_se <- function(tr){
  d <- agd9_k[trial==tr & type=="overall"]
  c(te=as.numeric(d$TE[1]), se=as.numeric(d$se[1]))
}

# per-draw loop
for(b in seq_len(B_BOOT)){
  # per-degree storage of trial-level inputs
  for(deg in DEGREE_SET){
    tstar <- targets_tbl %>% filter(degree==deg) %>% slice(1)
    target_vec <- as.list(tstar[, COVARS])
    TE_true <- as.numeric(tstar$TE_true)

    # trial-level te/se for two strategies
    te_overall <- rep(NA_real_, K); se_overall <- rep(NA_real_, K)
    te_eacc    <- rep(NA_real_, K); se_eacc    <- rep(NA_real_, K)

    coef_rows <- list()
    imp_rows <- list()

    for(i in seq_len(K)){
      tr <- plan$trial[i]
      ov <- get_overall_te_se(tr)
      te_overall[i] <- ov["te"]; se_overall[i] <- ov["se"]

      if(!isTRUE(plan$is_large[i])){
        # small trials: fallback to overall
        te_eacc[i] <- ov["te"]; se_eacc[i] <- ov["se"]
        next
      }

      # draw external mu/Sigma by bootstrap, then shrink and PSD repair
      bm <- base_match[[tr]]
      samp <- bm$sample
      idx <- sample.int(nrow(samp), size=nrow(samp), replace=TRUE)
      samp_b <- samp[idx, , drop=FALSE]
      mu_b <- colMeans(samp_b)
      S_b  <- stats::cov(as.matrix(samp_b))
      S_b <- (S_b + t(S_b))/2
      rownames(S_b) <- colnames(S_b) <- COVARS
      # shrink toward base
      S_b <- ETA_SIGMA * S_b + (1-ETA_SIGMA) * bm$Sigma
      S_b <- make_psd(S_b)
      # FH truncation + PSD
      trc <- truncate_corr_to_fh(mu_b, S_b)
      S_b <- trc$Sigma

      # take agd9, create TE_draw (optional perturb)
      agd_tr <- get_agd9_trial(tr)
      if(isTRUE(PERTURB_TE)){
        agd_tr$TE_draw <- ifelse(is.finite(agd_tr$TE) & is.finite(agd_tr$se) & agd_tr$se>0,
                                 rnorm(nrow(agd_tr), agd_tr$TE, agd_tr$se),
                                 NA_real_)
      } else {
        agd_tr$TE_draw <- NA_real_
      }

      # BLUP completion
      imp <- blup_impute_agd9_draw(agd_trial = agd_tr, Sigma = S_b, mu_vec = mu_b,
                                  random = RANDOM_BLUP, clip = CLIP_RANGE)
      agd_imp <- imp$agd_imp

      # WLS
      fit <- fit_wls_trial_draw(agd_imp = agd_imp, varset = VARSET, perturb_te = PERTURB_TE)
      if(isTRUE(fit$ok)){
        pred <- predict_wls_at_target(fit, target_vec_named = target_vec)
        te_eacc[i] <- pred$te; se_eacc[i] <- pred$se

        # save coef row
        cf <- coef(fit$fit)
        coef_rows[[length(coef_rows)+1]] <- tibble::tibble(
          draw=b, degree=deg, trial=tr,
          term = names(cf),
          estimate = as.numeric(cf)
        )
      } else {
        # fallback: overall
        te_eacc[i] <- ov["te"]; se_eacc[i] <- ov["se"]
      }

      imp_rows[[length(imp_rows)+1]] <- tibble::tibble(
        draw=b, degree=deg, trial=tr,
        imp_age65=imp$imp_n["age65"], imp_htn=imp$imp_n["htn"], imp_male=imp$imp_n["male"], imp_smoke=imp$imp_n["smoke"],
        clip_age65=imp$clip_n["age65"], clip_htn=imp$clip_n["htn"], clip_male=imp$clip_n["male"], clip_smoke=imp$clip_n["smoke"]
      )
    } # end trial loop

    # Meta pooling: 4 method labels
    # 1) FE overall-only
    fe0 <- fe_pool_closed(te_overall, se_overall)
    # 2) RE overall-only
    re0 <- re_pool_metagen_safe(te_overall, se_overall)
    # 3) FE EACC (Hybrid_WLS_large)
    fe1 <- fe_pool_closed(te_eacc, se_eacc)
    # 4) RE EACC (Hybrid_WLS_large)
    re1 <- re_pool_metagen_safe(te_eacc, se_eacc)

    one_row <- function(strategy, model, obj){
      tibble::tibble(
        draw=b, degree=deg,
        meta_strategy=strategy,
        meta_model=model,
        TE_true=TE_true,
        TE=obj$TE, lo=obj$lo, hi=obj$hi,
        tau2=obj$tau2, I2=obj$I2, Q=obj$Q,
        k_used=obj$k_used,
        ok=obj$ok
      )
    }

    meta_draws[[length(meta_draws)+1]] <- dplyr::bind_rows(
      one_row("Overall_only","FE",fe0),
      one_row("Overall_only","RE",re0),
      one_row("Hybrid_WLS_large","FE",fe1),
      one_row("Hybrid_WLS_large","RE",re1)
    )

    # save trial-level inputs used in meta for debugging
    trial_inputs_draws[[length(trial_inputs_draws)+1]] <- tibble::tibble(
      draw=b, degree=deg, trial=plan$trial,
      te_overall=te_overall, se_overall=se_overall,
      te_eacc=te_eacc, se_eacc=se_eacc,
      is_large=plan$is_large
    )

    if(length(coef_rows)) wls_coef_draws[[length(wls_coef_draws)+1]] <- dplyr::bind_rows(coef_rows)
    if(length(imp_rows)) impute_diag[[length(impute_diag)+1]] <- dplyr::bind_rows(imp_rows)
  } # end degree loop
} # end draw loop

meta_draws <- dplyr::bind_rows(meta_draws)
trial_inputs_draws <- dplyr::bind_rows(trial_inputs_draws)
wls_coef_draws <- dplyr::bind_rows(wls_coef_draws)
impute_diag <- dplyr::bind_rows(impute_diag)

readr::write_csv(meta_draws, file.path(OUT_DIR, "META_draws.csv"))
readr::write_csv(trial_inputs_draws, file.path(OUT_DIR, "META_trial_inputs_draws.csv"))
readr::write_csv(wls_coef_draws, file.path(OUT_DIR, "WLS_coef_draws.csv"))
readr::write_csv(impute_diag, file.path(OUT_DIR, "BLUP_impute_diag.csv"))

# -----------------------------
# 15) Summaries and figures
# -----------------------------
# Map to display labels
method_label <- function(strategy, model){
  if(strategy=="Hybrid_WLS_large" && model=="FE") "EACC-FE"
  else if(strategy=="Overall_only" && model=="FE") "FE"
  else if(strategy=="Hybrid_WLS_large" && model=="RE") "EACC-RE"
  else if(strategy=="Overall_only" && model=="RE") "RE"
  else paste(strategy, model)
}

meta_draws2 <- meta_draws %>%
  mutate(method = purrr::pmap_chr(list(meta_strategy, meta_model), method_label),
         cover = as.integer(isTRUE(ok) & is.finite(lo) & is.finite(hi) & (TE_true >= lo) & (TE_true <= hi)),
         err = TE - TE_true)

summary_tbl <- meta_draws2 %>%
  group_by(degree, method) %>%
  summarise(
    TE_true   = first(TE_true),
    mean_est  = mean(TE, na.rm = TRUE),

    # Empirical percentiles across FULLBOOT draws (used for EACC by default)
    lo_draw   = quantile(TE, 0.025, na.rm = TRUE),
    hi_draw   = quantile(TE, 0.975, na.rm = TRUE),

    # Conventional meta-analytic CI within each draw, then averaged across draws
    # (for Overall_only Meta, these are identical across draws because inputs do not change)
    lo_meta   = mean(lo, na.rm = TRUE),
    hi_meta   = mean(hi, na.rm = TRUE),

    mean_bias = mean(err, na.rm = TRUE),
    mse       = mean(err^2, na.rm = TRUE),
    rmse      = sqrt(mse),

    # Coverage computed from conventional within-draw meta CIs
    cover_meta    = mean(cover, na.rm = TRUE),
    mean_wid_meta = mean(hi - lo, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    # Key fix requested:
    # For traditional Meta (FE/RE): use the analytic meta CI, not the across-draw percentiles.
    # For EACC: keep empirical percentiles across FULLBOOT draws (consistent with the bootstrap-style workflow).
    lo = dplyr::if_else(method %in% c("FE", "RE"), lo_meta, lo_draw),
    hi = dplyr::if_else(method %in% c("FE", "RE"), hi_meta, hi_draw),

    # Keep the main coverage definition as cover_meta.
    cover = cover_meta,
    mean_wid = hi - lo
  ) %>%
  select(
    degree, method, TE_true,
    mean_est, lo, hi,
    mean_bias, mse, rmse, cover, mean_wid,
    lo_draw, hi_draw, lo_meta, hi_meta, cover_meta, mean_wid_meta
  ) %>%
  arrange(degree, factor(method, levels = c("EACC-FE", "FE", "EACC-RE", "RE")))

readr::write_csv(summary_tbl, file.path(OUT_DIR, "SUMMARY_minimal_example.csv"))

# matching diag plot
p_match <- match_diag %>%
  ggplot(aes(x=max_abs_dev, y=fh_fix_pairs)) +
  geom_point() +
  labs(x="Max abs marginal deviation in base match", y="Number of FH-truncated pairs",
       title="Base external matching diagnostics (large trials)")
ggsave(file.path(OUT_DIR, "fig_matching_diagnostics.png"), p_match, width=7.5, height=4.2, dpi=300)

# imputation diag plot
imp_sum <- impute_diag %>%
  group_by(trial) %>%
  summarise(
    clip_total = mean(clip_age65 + clip_htn + clip_male + clip_smoke),
    imp_total  = mean(imp_age65 + imp_htn + imp_male + imp_smoke),
    .groups="drop"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_imp <- imp_sum %>%
    ggplot(aes(x=imp_total, y=clip_total, label=trial)) +
    geom_point() +
    ggrepel::geom_text_repel(max.overlaps = 50) +
    labs(
      x="Mean number of imputed cells per draw",
      y="Mean number of clipped cells per draw",
      title="BLUP completion diagnostics (large trials)"
    ) +
    theme_minimal()
} else {
  p_imp <- imp_sum %>%
    ggplot(aes(x=imp_total, y=clip_total, label=trial)) +
    geom_point() +
    geom_text(vjust = -0.6, check_overlap = TRUE) +
    labs(
      x="Mean number of imputed cells per draw",
      y="Mean number of clipped cells per draw",
      title="BLUP completion diagnostics (large trials)"
    ) +
    theme_minimal()
}
ggsave(file.path(OUT_DIR, "fig_blup_diagnostics.png"), p_imp, width=7.5, height=4.2, dpi=300)

# -----------------------------
# 16) Figures: forest + 3-panel metrics
# -----------------------------
plot_metric <- function(metric, title, ref_line = NULL, digits = 4, xlim = NULL) {
  df <- summary_tbl %>%
    select(degree, method, value = dplyr::all_of(metric)) %>%
    mutate(method = factor(method, levels = c("EACC-FE","FE","EACC-RE","RE")))

  panels <- lapply(seq_along(DEGREE_SET), function(ii) {
    deg <- DEGREE_SET[ii]
    d <- df %>% filter(degree == deg) %>%
      mutate(label = sprintf(paste0("%.", digits, "f"), value))

    g <- ggplot(d, aes(y = method, x = value)) +
      geom_point() +
      geom_text(aes(label = label), hjust = -0.15, size = 3) +
      labs(
        title = paste0(title, " (degree = ", deg, " SD)"),
        x = title, y = NULL
      ) +
      theme_minimal() +
      coord_cartesian(clip = "off") +
      theme(
        plot.margin = margin(5.5, 35, 5.5, 5.5),
        plot.title = element_text(size = 10)
      )

    if (!is.null(ref_line)) {
      g <- g + geom_vline(xintercept = ref_line, linetype = "dashed", alpha = 0.6)
    }
    if (!is.null(xlim)) {
      g <- g + coord_cartesian(xlim = xlim, clip = "off")
    }
    if (ii > 1) {
      g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    g
  })

  patchwork::wrap_plots(panels, nrow = 1)
}

p_bias <- plot_metric("mean_bias", "Mean bias", ref_line = 0, digits = 4)
p_rmse <- plot_metric("rmse", "RMSE", ref_line = 0, digits = 4)
p_cov  <- plot_metric("cover", "Coverage", ref_line = 0.95, digits = 3, xlim = c(0,1))

ggsave(file.path(OUT_DIR, "fig_metric_bias.png"), p_bias, width=10.5, height=3.6, dpi=300)
ggsave(file.path(OUT_DIR, "fig_metric_rmse.png"), p_rmse, width=10.5, height=3.6, dpi=300)
ggsave(file.path(OUT_DIR, "fig_metric_cover.png"), p_cov,  width=10.5, height=3.6, dpi=300)

# forest plot for degree=0.8 (main)
deg_main <- 0.8
forest_tbl <- summary_tbl %>%
  filter(degree == deg_main) %>%
  mutate(method = factor(method, levels = c("EACC-FE","FE","EACC-RE","RE"))) %>%
  arrange(method)

p_forest <- forest_tbl %>%
  ggplot(aes(y = method, x = mean_est)) +
  geom_vline(xintercept = unique(forest_tbl$TE_true), linetype = "dashed") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18) +
  geom_point() +
  labs(
    x = "Pooled TE (probit scale)",
    y = NULL,
    title = paste0("Pooled results at degree = ", deg_main, " SD"),
    subtitle = paste0(
      "Points: mean pooled TE across draws. ",
      "Whiskers: for EACC, 2.5% and 97.5% percentiles across FULLBOOT draws; ",
      "for FE/RE, conventional 95% meta-analytic CI."
    )
  ) +
  theme_minimal()

ggsave(file.path(OUT_DIR, "fig_forest_degree08.png"), p_forest, width = 7.5, height = 3.8, dpi = 300)

# -----------------------------
# 17) One Excel workbook (convenience)
# -----------------------------
wb <- openxlsx::createWorkbook()
addSheet <- function(name, df){
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df)
}
addSheet("trial_plan", plan)
addSheet("trial_targets", targets)
addSheet("trial_alpha", alpha_tbl)
addSheet("cov_obs_generated", cov_obs_k)
addSheet("targets_truth", targets_tbl)
addSheet("agd9_generated", as.data.frame(agd9_k))
addSheet("match_diag_base", match_diag)
addSheet("META_draws", meta_draws)
addSheet("META_trial_inputs_draws", trial_inputs_draws)
addSheet("SUMMARY", summary_tbl)
openxlsx::saveWorkbook(wb, file.path(OUT_DIR, "minimal_example_outputs.xlsx"), overwrite=TRUE)

message("Done. Key outputs written to: ", OUT_DIR)
