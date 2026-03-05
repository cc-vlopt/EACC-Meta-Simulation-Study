# =============================================================================
# 10_run_simulation.R — Main Simulation Loop
# =============================================================================
# ADEMP Component: All — Orchestrates repetitions x scenarios x K-subsets
# =============================================================================

# --- Single repetition: generate data and run FULLBOOT for all scenarios ------
run_one_rep <- function(rep_id, rwd_std_local, pool_df_gen_local, p = NULL) {

  message(sprintf("\n==== [REP %d / %d] ====", rep_id, N_REP))

  # 9 base-case scenarios (pL x pV)
  base_grid <- tidyr::expand_grid(pL = pL_SET, pV = pV_SET) %>%
    dplyr::mutate(scene_id = dplyr::row_number())

  rep_leaves <- list()
  jj <- 0L

  for (i in seq_len(nrow(base_grid))) {
    g        <- base_grid[i, ]
    scene_id <- g$scene_id
    pL       <- g$pL
    pV       <- g$pV

    # Scene-level seed (generation stream)
    set.seed(as.integer(
      SEED_GEN + rep_id * 1000000 + scene_id * 10000 +
      as.integer(round(100 * pL)) * 10 + as.integer(round(100 * pV))
    ))

    ns  <- calc_nL_nS(K_MAX, pL)
    nL  <- ns$nL; nS <- ns$nS
    trial_ids    <- paste0("T", seq_len(K_MAX))
    is_large_vec <- c(rep(TRUE, nL), rep(FALSE, nS))
    N_vec        <- draw_trial_sizes(nL, nS)
    alpha_vec    <- purrr::map_dbl(is_large_vec, draw_alpha_k)

    plan_univ <- tibble::tibble(
      rep = rep_id, scene_id = scene_id, pL = pL, pV = pV,
      trial = trial_ids,
      trial_order = seq_len(K_MAX),
      N_total = as.integer(N_vec),
      alpha_k = as.numeric(alpha_vec),
      is_large_design = is_large_vec
    )

    targets_univ <- draw_trial_targets(K_MAX, pV) %>%
      dplyr::left_join(
        plan_univ %>% dplyr::select(trial, N_total, alpha_k, is_large_design),
        by = "trial"
      ) %>%
      dplyr::mutate(rep = rep_id, scene_id = scene_id,
                    pL = pL, pV = pV) %>%
      dplyr::relocate(rep, scene_id, pL, pV, trial)

    # Save scene-level inputs
    sdir <- scene_dir(rep_id, scene_id, pL, pV)
    safe_mkdir(sdir)
    if (isTRUE(SAVE_SCENE_INPUTS)) {
      data.table::fwrite(data.table::as.data.table(plan_univ),
                         file.path(sdir, "universe_plan.csv"))
      data.table::fwrite(data.table::as.data.table(targets_univ),
                         file.path(sdir, "universe_trial_targets.csv"))
    }

    # Simulate all K_MAX trials
    cov_list <- list(); arm_list <- list(); agd_list <- list()
    ipd_list <- list()

    for (tr in trial_ids) {
      rowt <- targets_univ %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)
      targ <- c(age65 = rowt$age65, htn = rowt$htn,
                male = rowt$male, smoke = rowt$smoke)
      rowp <- plan_univ %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)

      sim <- gen_one_trial_AB(
        trial_id = tr,
        is_large = isTRUE(rowp$is_large_design),
        N_total  = rowp$N_total,
        alpha_k  = rowp$alpha_k,
        targets  = targ,
        pool_df  = pool_df_gen_local
      )

      cov_list[[tr]] <- sim$cov_obs
      arm_list[[tr]] <- sim$arm_counts
      agd_list[[tr]] <- sim$agd9
      if (isTRUE(SAVE_IPD) && !is.null(sim$ipd)) ipd_list[[tr]] <- sim$ipd
    }

    cov_obs_univ <- data.table::rbindlist(cov_list, use.names = TRUE,
                                          fill = TRUE)
    cov_obs_univ <- as.data.frame(cov_obs_univ) %>%
      dplyr::mutate(rep = rep_id, scene_id = scene_id,
                    pL = pL, pV = pV) %>%
      dplyr::relocate(rep, scene_id, pL, pV, trial)

    arm_counts_univ <- data.table::rbindlist(arm_list, use.names = TRUE,
                                             fill = TRUE)
    arm_counts_univ <- as.data.frame(arm_counts_univ) %>%
      dplyr::mutate(rep = rep_id, scene_id = scene_id,
                    pL = pL, pV = pV) %>%
      dplyr::relocate(rep, scene_id, pL, pV, trial)

    agd9_univ <- data.table::rbindlist(agd_list, use.names = TRUE, fill = TRUE)
    agd9_univ <- as.data.frame(agd9_univ) %>%
      dplyr::mutate(rep = rep_id, scene_id = scene_id,
                    pL = pL, pV = pV) %>%
      dplyr::relocate(rep, scene_id, pL, pV, trial)

    if (isTRUE(SAVE_SCENE_INPUTS)) {
      data.table::fwrite(data.table::as.data.table(cov_obs_univ),
                         file.path(sdir, "universe_cov_obs.csv"))
      data.table::fwrite(data.table::as.data.table(arm_counts_univ),
                         file.path(sdir, "universe_arm_counts.csv"))
      data.table::fwrite(data.table::as.data.table(agd9_univ),
                         file.path(sdir, "universe_agd9.csv"))
    }

    if (isTRUE(SAVE_IPD) && length(ipd_list)) {
      ipd_univ <- data.table::rbindlist(ipd_list, use.names = TRUE,
                                        fill = TRUE)
      data.table::fwrite(
        ipd_univ,
        file.path(OUT_ROOT_RUN,
                  sprintf("IPD_rep%03d_scene%02d.csv", rep_id, scene_id))
      )
    }

    # K-subset loop
    for (K in K_SUB) {
      trials_k <- paste0("T", seq_len(K))
      plan_k <- plan_univ %>% dplyr::filter(trial %in% trials_k)
      cov_k  <- cov_obs_univ %>% dplyr::filter(trial %in% trials_k)
      agd_k  <- agd9_univ %>% dplyr::filter(trial %in% trials_k)

      tt <- compute_targets_truth_for_K(
        cov_obs_k = cov_k, plan_k = plan_k,
        meta = list(rep = rep_id, scene_id = scene_id,
                    pL = pL, pV = pV, K = K)
      )

      kdir0 <- k_dir(rep_id, scene_id, pL, pV, K)
      safe_mkdir(kdir0)
      if (isTRUE(SAVE_K_INPUTS)) {
        data.table::fwrite(data.table::as.data.table(plan_k),
                           file.path(kdir0, "plan_k.csv"))
        data.table::fwrite(data.table::as.data.table(cov_k),
                           file.path(kdir0, "cov_obs_k.csv"))
        data.table::fwrite(data.table::as.data.table(agd_k),
                           file.path(kdir0, "agd9_k.csv"))
        data.table::fwrite(data.table::as.data.table(tt$targets),
                           file.path(kdir0, "targets_by_degree.csv"))
        data.table::fwrite(data.table::as.data.table(tt$truth),
                           file.path(kdir0, "truth_by_degree.csv"))
      }

      jj <- jj + 1L
      message(sprintf("  scene=%d pL=%.2f pV=%.2f K=%d (%d/%d)",
                      scene_id, pL, pV, K, jj,
                      nrow(base_grid) * length(K_SUB)))

      if (!is.null(p)) {
        p(sprintf("rep %d  scene %d  K=%d  (%d/%d)",
                  rep_id, scene_id, K, jj,
                  nrow(base_grid) * length(K_SUB)))
      }

      rep_case <- run_fullboot_meta_one_case(
        rep_id = rep_id, scene_id = scene_id, pL = pL, pV = pV, K = K,
        plan_k = plan_k, cov_obs_k = cov_k, agd9_k = agd_k,
        targets_by_delta = tt$targets, truth_by_delta = tt$truth,
        rwd_std_local = rwd_std_local
      )
      rep_leaves[[length(rep_leaves) + 1]] <- rep_case
    }
  }

  data.table::rbindlist(rep_leaves, use.names = TRUE, fill = TRUE)
}

# --- Aggregate performance measures across repetitions ------------------------
# ADEMP Component: Performance measures (P)
compute_summary <- function(REP_LEVEL_ALL, case_keys) {
  REP_LEVEL_ALL[, {
    ok <- is.finite(point_est) & is.finite(TE_true)
    dd <- .SD[ok]
    n_eff <- nrow(dd)

    if (n_eff == 0) {
      list(n_rep = 0L,
           mean_trut = NA_real_, mean_est = NA_real_, mean_bias = NA_real_,
           mse = NA_real_, rmse = NA_real_,
           cover = NA_real_, mcse_cover = NA_real_, mean_wid = NA_real_)
    } else {
      tr  <- dd$TE_true
      est <- dd$point_est
      cv  <- dd$cover
      wd  <- dd$width

      mean_trut  <- mean(tr)
      mean_est   <- mean(est)
      mean_bias  <- mean(est - tr)
      mse_val    <- mean((est - tr)^2)
      rmse_val   <- sqrt(mse_val)
      cover_rate <- mean(cv, na.rm = TRUE)
      n_cover    <- sum(is.finite(cv))
      mcse <- if (n_cover > 0) {
        sqrt(cover_rate * (1 - cover_rate) / n_cover)
      } else {
        NA_real_
      }
      mean_wid <- mean(wd, na.rm = TRUE)

      list(n_rep = n_eff,
           mean_trut = mean_trut, mean_est = mean_est,
           mean_bias = mean_bias,
           mse = mse_val, rmse = rmse_val,
           cover = cover_rate, mcse_cover = mcse,
           mean_wid = mean_wid)
    }
  }, by = case_keys]
}
