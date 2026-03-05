# =============================================================================
# 03_data_generation.R — Data-Generating Process and Trial Simulation
# =============================================================================
# ADEMP Component: Data-generating mechanisms (D)
#
# Generates individual-level binary outcomes from a probit model with
# prognostic and effect-modifying components, then compresses to
# trial-level 9-row incomplete aggregate tables (agd9).
# =============================================================================

# --- Draw trial-specific baseline intercept -----------------------------------
draw_alpha_k <- function(is_large) {
  sigma <- if (isTRUE(is_large)) 0.10 else 0.05
  rtruncnorm_1(mu = mu_a, sigma = sigma,
               a = mu_a - 2.5 * sigma,
               b = mu_a + 2.5 * sigma)
}

# --- Generate one trial (A vs. B) and return aggregate data ------------------
gen_one_trial_AB <- function(trial_id, is_large, N_total, alpha_k,
                             targets, pool_df) {
  N_A <- floor(N_total / 2)
  N_B <- N_total - N_A

  # Sample baseline covariates via raking
  samp <- sample_with_raking_4(pool_df = pool_df, N = N_total, targets = targets)
  samp <- data.table::as.data.table(samp)

  # Randomize treatment assignment
  perm <- sample.int(N_total, N_total, replace = FALSE)
  samp <- samp[perm]
  samp[, Trt := ifelse(.I <= N_A, "A", "B")]
  samp[, t := ifelse(Trt == "B", 1L, 0L)]

  # Generate binary outcomes from probit model
  X <- as.matrix(samp[, ..COVARS])
  lp_main <- as.numeric(alpha_k + X %*% beta[COVARS])
  lp_trt  <- samp$t * (tau_AB + as.numeric(
    as.matrix(samp[, ..EM_VARS]) %*% delta[EM_VARS]
  ))
  lp <- lp_main + lp_trt
  p  <- inv_probit(lp)
  samp[, y := rbinom(.N, 1, p)]

  # Observed covariate means
  cov_obs <- samp[, lapply(.SD, mean), .SDcols = COVARS][1]
  cov_obs[, `:=`(trial = trial_id, N = N_total)]
  data.table::setcolorder(cov_obs, c("trial", "N", COVARS))

  # Arm-level counts
  samp[, trial := trial_id]
  arm_counts <- samp[, .(n = .N, events = sum(y), risk = mean(y)),
                     by = .(trial, Trt)]

  # Construct 9-row aggregate data table
  agd9 <- build_agd9(samp, trial_id)

  # Optional: return IPD
  ipd_out <- NULL
  if (isTRUE(SAVE_IPD)) {
    ipd_out <- samp[, c("trial", "Trt", "t", COVARS, "y"), with = FALSE]
  }

  list(cov_obs = cov_obs, arm_counts = arm_counts, agd9 = agd9, ipd = ipd_out)
}

# --- Build 9-row incomplete aggregate table for one trial ---------------------
build_agd9 <- function(samp, trial_id) {
  agg_two_arms <- function(dt) {
    arm <- dt[, .(r = sum(y), n = .N), by = .(Trt)]
    if (nrow(arm) < 2) return(NULL)
    list(rB = arm[Trt == "B"]$r, nB = arm[Trt == "B"]$n,
         rA = arm[Trt == "A"]$r, nA = arm[Trt == "A"]$n)
  }

  make_row <- function(type, dt_sub, x_age65, x_htn, x_male, x_smoke) {
    cnt <- agg_two_arms(dt_sub)
    if (is.null(cnt)) {
      te  <- list(te = NA_real_, se = NA_real_)
      cnt <- list(rB = NA_real_, nB = NA_integer_, rA = NA_real_, nA = NA_integer_)
    } else {
      te <- te_probit_from_counts(cnt$rB, cnt$nB, cnt$rA, cnt$nA)
    }
    data.table::data.table(
      trial = trial_id, type = type,
      x_age65 = x_age65, x_htn = x_htn, x_male = x_male, x_smoke = x_smoke,
      TE = te$te, se = te$se,
      rB = cnt$rB, nB = cnt$nB, rA = cnt$rA, nA = cnt$nA
    )
  }

  # Row 1: Overall
  rows <- list(
    make_row("overall", samp,
             mean(samp$age65), mean(samp$htn),
             mean(samp$male),  mean(samp$smoke))
  )

  # Rows 2–9: Single-covariate strata (1/0 for each covariate)
  for (v in COVARS) {
    for (val in c(1, 0)) {
      dt_sub  <- samp[get(v) == val]
      x_vals  <- stats::setNames(rep(NA_real_, 4), COVARS)
      x_vals[v] <- val
      rows[[length(rows) + 1]] <- make_row(
        paste0(v, "=", val), dt_sub,
        x_vals["age65"], x_vals["htn"], x_vals["male"], x_vals["smoke"]
      )
    }
  }

  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

# --- Target population generation and true estimand ---------------------------
make_target_probs_equal <- function(base_p, delta_sd, clip = c(0.02, 0.98)) {
  if (is.data.frame(base_p)) {
    p0 <- as.numeric(unlist(base_p[1, COVARS], use.names = FALSE))
  } else {
    p0 <- as.numeric(base_p[COVARS])
  }
  names(p0) <- COVARS
  p0 <- pmin(pmax(p0, 1e-6), 1 - 1e-6)
  z0 <- qnorm(p0)
  d  <- rep(1, length(COVARS))
  d  <- d / sqrt(sum(d^2))
  z_star <- z0 + delta_sd * d
  p_star <- pnorm(z_star)
  p_star <- pmin(pmax(p_star, clip[1]), clip[2])
  stats::setNames(p_star, COVARS)
}

truth_TE_target <- function(p_star) {
  tau_AB + sum(delta[EM_VARS] * p_star[EM_VARS])
}

compute_targets_truth_for_K <- function(cov_obs_k, plan_k, meta) {
  cov_w <- cov_obs_k %>%
    dplyr::left_join(
      plan_k %>% dplyr::select(trial, N_total_plan = N_total),
      by = "trial"
    ) %>%
    dplyr::mutate(N_w = dplyr::coalesce(.data$N_total_plan, as.numeric(.data$N)))

  ref_p <- cov_w %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(COVARS),
                     ~ weighted.mean(.x, w = .data$N_w, na.rm = TRUE))
    )

  out_targets <- list()
  out_truth   <- list()
  for (deg in DEGREE_SET) {
    p_star  <- make_target_probs_equal(ref_p[1, ], deg)
    te_true <- truth_TE_target(p_star)
    out_targets[[as.character(deg)]] <- tibble::tibble(
      rep = meta$rep, scene_id = meta$scene_id,
      pL = meta$pL, pV = meta$pV, K = meta$K,
      degree = deg, !!!as.list(p_star)
    )
    out_truth[[as.character(deg)]] <- tibble::tibble(
      rep = meta$rep, scene_id = meta$scene_id,
      pL = meta$pL, pV = meta$pV, K = meta$K,
      degree = deg, TE_true = te_true
    )
  }
  list(
    targets = dplyr::bind_rows(out_targets),
    truth   = dplyr::bind_rows(out_truth)
  )
}

# --- Helper: draw trial sizes and targets for a scenario ----------------------
calc_nL_nS <- function(K, pL) {
  nL <- as.integer(round_half_up(K * pL))
  nL <- max(0L, min(nL, as.integer(K)))
  nS <- as.integer(K - nL)
  list(nL = nL, nS = nS)
}

draw_trial_sizes <- function(nL, nS) {
  N_large <- if (nL > 0) sample(8000:16000, size = nL, replace = TRUE) else integer(0)
  N_small <- if (nS > 0) sample(200:300,    size = nS, replace = TRUE) else integer(0)
  c(N_large, N_small)
}

draw_trial_targets <- function(K, pV) {
  tibble::tibble(
    trial = paste0("T", seq_len(K)),
    age65 = rep(pV, K),
    htn   = runif(K, 0.4, 0.6),
    male  = runif(K, 0.4, 0.6),
    smoke = runif(K, 0.4, 0.6)
  )
}
