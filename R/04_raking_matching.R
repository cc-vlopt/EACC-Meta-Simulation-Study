# =============================================================================
# 04_raking_matching.R — Raking-Based Sampling and External Matching
# =============================================================================
# ADEMP Component: Methods (M) — Step 1 of EACC-Meta workflow
# =============================================================================

# --- Raking-weighted sampling from external pool ------------------------------
sample_with_raking_4 <- function(pool_df, N, targets,
                                 max_try = 10,
                                 tol = c(age65 = 0.05, htn = 0.10,
                                         male = 0.10, smoke = 0.10),
                                 ipf_maxit = 300, ipf_tol = 1e-7) {
  vars <- COVARS
  stopifnot(all(vars %in% names(pool_df)))

  df <- pool_df[, vars, drop = FALSE]
  for (v in vars) df[[v]] <- as.integer(df[[v]] %in% c(1L, TRUE))
  df <- stats::na.omit(df)
  if (nrow(df) < 100) stop("External pool has fewer than 100 usable observations.")

  base <- df
  base$.id <- seq_len(nrow(base))
  for (v in vars) base[[v]] <- factor(base[[v]], levels = c(0, 1))
  dsgn <- survey::svydesign(ids = ~.id, data = base, weights = ~1)
  Npop <- nrow(base)

  pop_margins <- lapply(vars, function(v) {
    data.frame(
      tmp  = factor(0:1),
      Freq = c(1 - targets[[v]], targets[[v]]) * Npop
    ) %>% dplyr::rename(!!v := tmp)
  })

  adj <- try(survey::rake(
    dsgn,
    sample.margins  = lapply(vars, function(v) as.formula(paste0("~", v))),
    population.margins = pop_margins,
    control = list(maxit = ipf_maxit, epsilon = ipf_tol)
  ), silent = TRUE)

  use_uniform <- inherits(adj, "try-error")
  if (use_uniform) {
    w <- rep(1, nrow(df))
  } else {
    w <- as.numeric(weights(adj))
    w[!is.finite(w) | w <= 0] <- min(w[w > 0], na.rm = TRUE)
    qq <- stats::quantile(w, c(0.01, 0.99), na.rm = TRUE)
    w  <- pmin(pmax(w, qq[1]), qq[2])
  }

  best <- NULL; best_dev <- Inf
  for (t in seq_len(max_try)) {
    id   <- sample.int(nrow(df), size = N, replace = TRUE, prob = w)
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

# --- Choose external sample size for matching ---------------------------------
choose_N_match <- function(N_total, mode = "equal_trial_n",
                           factor = 5, min_n = 200L, max_n = 200000L) {
  N_total <- as.integer(N_total)
  if (!is.finite(N_total) || N_total <= 0) return(min_n)
  if (mode == "equal_trial_n") return(N_total)
  as.integer(max(min_n, min(max_n, round(N_total * factor))))
}

# --- External base matching with FH bounds and tolerance relaxation -----------
match_external_base_4 <- function(pool_df, target_p, N_match, seed = NULL,
                                  tol0 = c(age65 = 0.05, htn = 0.10,
                                           male = 0.10, smoke = 0.10),
                                  relax_mult = c(1.0, 1.2, 1.5, 2.0),
                                  fh_eps = 1e-8) {
  if (!is.null(seed)) set.seed(seed)

  samp      <- NULL
  fh_ok     <- FALSE
  used_mult <- NA_real_
  fh_report <- NULL

  # Attempt matching with tolerance relaxation until FH bounds pass
  for (mult in relax_mult) {
    tol_use <- tol0 * mult
    attempt <- 0L
    samp    <- NULL

    while (attempt < MAX_TRY_MATCH) {
      attempt <- attempt + 1L
      samp <- tryCatch(
        sample_with_raking_4(pool_df = pool_df,
                             N = as.integer(N_match),
                             targets = target_p,
                             tol = tol_use),
        error = function(e) NULL
      )
      if (!is.null(samp)) break
    }
    if (is.null(samp)) next

    samp_use <- samp[, COVARS, drop = FALSE]
    mu0 <- colMeans(samp_use)
    S0  <- stats::cov(as.matrix(samp_use))
    S0  <- (S0 + t(S0)) / 2
    rownames(S0) <- colnames(S0) <- COVARS

    cor0 <- suppressWarnings(stats::cov2cor(S0))
    chk  <- fh_check_cor_matrix(cor0, target_p = target_p,
                                 vars = COVARS, eps = fh_eps)
    fh_ok     <- isTRUE(chk$ok)
    fh_report <- chk$report
    used_mult <- mult

    if (fh_ok) {
      S0 <- make_psd(S0)
      return(list(sample = samp_use, mu = mu0, Sigma = S0,
                  fallback = FALSE, fh_ok = TRUE,
                  fh_relax_mult = used_mult,
                  fh_action = "pass", fh_report = fh_report))
    }
  }

  # Fallback: random sample (no raking)
  fallback <- FALSE
  if (is.null(samp)) {
    fallback <- TRUE
    set.seed((seed %||% 1L) + 99991L)
    idx  <- sample.int(nrow(pool_df), size = as.integer(N_match),
                       replace = (N_match > nrow(pool_df)))
    samp <- pool_df[idx, , drop = FALSE]
  }

  samp_use <- samp[, COVARS, drop = FALSE]
  mu0 <- colMeans(samp_use)
  S0  <- stats::cov(as.matrix(samp_use))
  S0  <- (S0 + t(S0)) / 2
  rownames(S0) <- colnames(S0) <- COVARS
  cor0 <- suppressWarnings(stats::cov2cor(S0))

  chk2 <- fh_check_cor_matrix(cor0, target_p = target_p,
                                vars = COVARS, eps = fh_eps)
  fh_report <- chk2$report

  if (!isTRUE(chk2$ok)) {
    # Clip correlation to FH bounds and rebuild Sigma
    cor_clip <- fh_clip_cor_to_bounds(cor0, target_p = target_p,
                                      vars = COVARS, eps = fh_eps)
    S_fix <- sigma_from_cor_and_targetp(cor_clip, target_p = target_p,
                                         vars = COVARS)
    return(list(sample = samp_use, mu = mu0, Sigma = S_fix,
                fallback = fallback, fh_ok = FALSE,
                fh_relax_mult = used_mult,
                fh_action = "clip_fix", fh_report = fh_report))
  } else {
    S0 <- make_psd(S0)
    return(list(sample = samp_use, mu = mu0, Sigma = S0,
                fallback = fallback, fh_ok = TRUE,
                fh_relax_mult = used_mult,
                fh_action = "pass_after_fallback", fh_report = fh_report))
  }
}
