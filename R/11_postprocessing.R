# =============================================================================
# 11_postprocessing.R — Rubin's Rules Post-Processing
# =============================================================================
# ADEMP Component: Performance measures (P)
#
# Combines within-draw SE (from trial-level inputs) with between-draw
# variance using Rubin's multiple imputation rules.
# =============================================================================

# --- Apply Rubin's rules to one case directory --------------------------------
rubin_one_case_dir <- function(meta_path) {
  vdir     <- dirname(meta_path)
  tin_path <- file.path(vdir, "META_trial_inputs_draws.csv")
  if (!file.exists(tin_path)) return(NULL)

  meta_dt <- tryCatch(
    data.table::fread(meta_path, select = c(
      "rep", "scene_id", "pL", "pV", "K", "pool", "varset",
      "degree", "meta_strategy", "meta_model",
      "draw", "TE_true", "TE_est_draw", "tau2"
    )),
    error = function(e) NULL
  )
  tin_dt <- tryCatch(
    data.table::fread(tin_path, select = c(
      "draw", "degree", "meta_strategy", "trial", "se"
    )),
    error = function(e) NULL
  )
  if (is.null(meta_dt) || is.null(tin_dt)) return(NULL)

  tin_dt <- tin_dt[is.finite(se) & se > 0]
  if (nrow(tin_dt) == 0) return(NULL)

  # Within-draw SE for FE
  se_fe <- tin_dt[, .(
    SE_est_draw = sqrt(1 / sum(1 / (se^2)))
  ), by = .(draw, degree, meta_strategy)]
  se_fe[, meta_model := "FE"]

  # Within-draw SE for RE (using tau2)
  tau_dt <- meta_dt[
    meta_model != "FE" & is.finite(tau2),
    .(draw, degree, meta_strategy, tau2)
  ]
  if (nrow(tau_dt) > 0) {
    data.table::setkey(tau_dt, draw, degree, meta_strategy)
    tau_dt <- unique(tau_dt, by = data.table::key(tau_dt))
  }

  tin_tau <- merge(tin_dt, tau_dt,
                   by = c("draw", "degree", "meta_strategy"),
                   all.x = TRUE)
  tin_tau[!is.finite(tau2), tau2 := 0]

  se_re <- tin_tau[, .(
    SE_est_draw = sqrt(1 / sum(1 / (se^2 + tau2)))
  ), by = .(draw, degree, meta_strategy)]
  se_re[, meta_model := "RE_REML"]

  se_map <- rbind(se_fe, se_re, fill = TRUE)

  # Merge SE back
  meta_dt2 <- merge(meta_dt, se_map,
                    by = c("draw", "degree", "meta_strategy", "meta_model"),
                    all.x = TRUE)

  # Filter valid rows
  meta_dt2 <- meta_dt2[
    is.finite(TE_est_draw) &
    is.finite(SE_est_draw) & SE_est_draw > 0 &
    is.finite(TE_true)
  ]
  if (nrow(meta_dt2) == 0) return(NULL)

  # Rubin's combining rules
  rep_rubin <- meta_dt2[, {
    m    <- .N
    Qbar <- mean(TE_est_draw)
    Ubar <- mean(SE_est_draw^2)
    Bvar <- if (m > 1) var(TE_est_draw) else 0

    Tvar   <- Ubar + (1 + 1 / m) * Bvar
    se_tot <- sqrt(Tvar)

    nu <- if (Bvar > 0) {
      (m - 1) * (1 + Ubar / ((1 + 1 / m) * Bvar))^2
    } else {
      1e9
    }

    crit <- qt(0.975, df = nu)
    lo   <- Qbar - crit * se_tot
    hi   <- Qbar + crit * se_tot
    tt   <- TE_true[1]

    list(
      point_est   = Qbar,
      ci_lo       = lo,
      ci_hi       = hi,
      width       = hi - lo,
      cover       = as.integer(lo <= tt && tt <= hi),
      n_draw_used = m
    )
  }, by = .(rep, scene_id, pL, pV, K, pool, varset,
            degree, meta_strategy, meta_model, TE_true)]

  rep_rubin
}

# --- Run Rubin post-processing across all case directories --------------------
run_rubin_postprocessing <- function(root_dir) {
  message("\n[POST] Rubin/MI post-processing: scanning META_draws.csv ...")

  meta_files <- list.files(
    root_dir,
    pattern = "META_draws\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  message("[POST] Found ", length(meta_files), " META_draws files")
  if (length(meta_files) == 0) {
    stop("No META_draws.csv files found. Check SAVE_META_DRAWS setting.")
  }

  out_rep_path <- file.path(root_dir, "REP_LEVEL_ALL_RUBIN.csv")
  if (file.exists(out_rep_path)) file.remove(out_rep_path)

  n_ok <- 0L
  for (i in seq_along(meta_files)) {
    dt <- rubin_one_case_dir(meta_files[i])
    if (!is.null(dt) && nrow(dt) > 0) {
      data.table::fwrite(
        dt, out_rep_path,
        append = file.exists(out_rep_path),
        col.names = !file.exists(out_rep_path)
      )
      n_ok <- n_ok + 1L
    }
    if (i %% 2000 == 0) {
      message(sprintf("[POST] %d / %d processed, ok_cases = %d",
                      i, length(meta_files), n_ok))
    }
  }

  if (!file.exists(out_rep_path)) {
    stop("Rubin post-processing produced no output. ",
         "Check META_trial_inputs_draws.csv files.")
  }

  REP_LEVEL_RUBIN <- data.table::fread(out_rep_path)
  message("[POST] Rubin REP_LEVEL rows = ", nrow(REP_LEVEL_RUBIN))

  # Aggregate summary
  CASE_KEYS <- c("scene_id", "pL", "pV", "K", "pool", "varset",
                 "degree", "meta_strategy", "meta_model")

  SUMMARY_RUBIN <- compute_summary(REP_LEVEL_RUBIN, CASE_KEYS)
  data.table::setcolorder(
    SUMMARY_RUBIN,
    c(CASE_KEYS, "n_rep", "mean_trut", "mean_est", "mean_bias",
      "mse", "rmse", "cover", "mcse_cover", "mean_wid")
  )

  out_sum_path <- file.path(root_dir, "SUMMARY_200reps_RUBIN.csv")
  data.table::fwrite(SUMMARY_RUBIN, out_sum_path)

  message("\n[POST] Rubin post-processing complete:")
  message("  - ", out_rep_path)
  message("  - ", out_sum_path)

  invisible(list(rep_level = REP_LEVEL_RUBIN, summary = SUMMARY_RUBIN))
}
