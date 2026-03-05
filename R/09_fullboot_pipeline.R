# =============================================================================
# 09_fullboot_pipeline.R — Full-Process Bootstrap (FULLBOOT) Pipeline
# =============================================================================
# ADEMP Component: Methods (M) — Step 5 of EACC-Meta workflow
#
# Captures chained uncertainty from external matching through covariate
# completion, regression-based extrapolation, and meta-analytic pooling.
# =============================================================================

run_fullboot_meta_one_case <- function(rep_id, scene_id, pL, pV, K,
                                       plan_k, cov_obs_k, agd9_k,
                                       targets_by_delta, truth_by_delta,
                                       rwd_std_local) {

  plan_k <- plan_k %>%
    dplyr::mutate(is_large = (as.integer(N_total) >= as.integer(N_LARGE_CUTOFF)))

  trials <- plan_k$trial
  xcols  <- paste0("x_", COVARS)

  overall_tbl <- agd9_k %>%
    dplyr::filter(type == "overall") %>%
    dplyr::select(trial, TE_overall = TE, se_overall = se,
                  dplyr::all_of(xcols))

  kroot <- k_dir(rep_id, scene_id, pL, pV, K)
  safe_mkdir(kroot)

  degrees      <- sort(unique(targets_by_delta$degree))
  out_rep_list <- list()

  for (pool in MATCH_POOLS) {
    pool_df <- rwd_std_local[[pool]]
    pdir    <- file.path(kroot, sprintf("pool_%s", pool))
    safe_mkdir(pdir)

    # ------ Base matching for large trials ------------------------------------
    base_match     <- list()
    base_sigma_ref <- list()

    for (tr in trials) {
      rowp <- plan_k %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)
      if (!isTRUE(rowp$is_large)) next

      rowc     <- cov_obs_k %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)
      target_p <- stats::setNames(
        as.numeric(unlist(rowc[1, COVARS], use.names = FALSE)), COVARS
      )
      N_match_tr <- choose_N_match(rowp$N_total, mode = "equal_trial_n")
      seed_use   <- as.integer(
        SEED_EXT + rep_id * 1000000 + scene_id * 10000 +
        K * 100 + as.integer(gsub("\\D", "", tr)) + nchar(pool) * 1000
      )

      bm <- match_external_base_4(pool_df, target_p, N_match_tr, seed = seed_use)
      base_match[[tr]]     <- bm
      base_sigma_ref[[tr]] <- bm$Sigma
    }

    # Save base match diagnostics
    if (isTRUE(SAVE_BASE_MATCH)) {
      save_base_match_diagnostics(base_match, pdir)
    }

    # ------ Variable set loop -------------------------------------------------
    for (vn in names(VARSETS)) {
      varset <- VARSETS[[vn]]
      vdir   <- file.path(pdir, sprintf("varset_%s", vn))
      safe_mkdir(vdir)

      # Collectors for diagnostics
      coef_rows        <- list(); cc    <- 0L
      pred_rows        <- list(); pp    <- 0L
      trial_input_rows <- list(); ti    <- 0L
      imp_diag_rows    <- list(); i_imp <- 0L
      sigma_diag_rows  <- list(); ss    <- 0L
      draw_rows        <- list(); rr    <- 0L

      # ------ Bootstrap draw loop ---------------------------------------------
      for (b in seq_len(B_BOOT)) {
        set.seed(as.integer(
          SEED_BOOT + rep_id * 1000000 + scene_id * 10000 +
          K * 100 + b + nchar(pool) * 777 + nchar(vn) * 33
        ))

        # Bootstrap Sigma for each large trial
        Sigma_b_map <- list()
        for (tr in trials) {
          rowp <- plan_k %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)
          if (!isTRUE(rowp$is_large)) next
          bm <- base_match[[tr]]
          if (is.null(bm)) next

          S0 <- bm$sample
          n0 <- nrow(S0)
          if (!is.finite(n0) || n0 < 5) next

          idx_ext <- sample.int(n0, size = n0, replace = TRUE)
          Sb    <- S0[idx_ext, , drop = FALSE]
          mu_b  <- colMeans(Sb)
          Sig_b <- cov(as.matrix(Sb))
          Sig_b <- (Sig_b + t(Sig_b)) / 2
          rownames(Sig_b) <- colnames(Sig_b) <- COVARS

          # Shrinkage toward base Sigma
          Sig_ref <- base_sigma_ref[[tr]]
          if (!is.null(Sig_ref) && is.finite(ETA_SIGMA) && ETA_SIGMA < 1) {
            Sig_b <- ETA_SIGMA * Sig_b + (1 - ETA_SIGMA) * Sig_ref
            Sig_b <- (Sig_b + t(Sig_b)) / 2
          }
          Sig_b <- make_psd(Sig_b)

          # FH bounds check on bootstrap Sigma
          rowc <- cov_obs_k %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)
          target_p_tr <- stats::setNames(
            as.numeric(unlist(rowc[1, COVARS], use.names = FALSE)), COVARS
          )
          corS <- suppressWarnings(stats::cov2cor(Sig_b))
          chkS <- fh_check_cor_matrix(corS, target_p = target_p_tr,
                                       vars = COVARS, eps = 1e-8)
          if (!isTRUE(chkS$ok)) {
            cor_clip <- fh_clip_cor_to_bounds(corS, target_p = target_p_tr,
                                              vars = COVARS, eps = 1e-8)
            Sig_b <- sigma_from_cor_and_targetp(cor_clip,
                                                 target_p = target_p_tr,
                                                 vars = COVARS)
          }

          Sigma_b_map[[tr]] <- list(mu = mu_b, Sigma = Sig_b)

          if (isTRUE(SAVE_SIGMA_DIAG)) {
            S    <- Sig_b
            corD <- suppressWarnings(stats::cov2cor(S))
            ss   <- ss + 1L
            sigma_diag_rows[[ss]] <- data.table::data.table(
              draw = b, trial = tr,
              var_age65 = S["age65", "age65"], var_htn = S["htn", "htn"],
              var_male  = S["male", "male"],   var_smoke = S["smoke", "smoke"],
              cor_age65_htn   = corD["age65", "htn"],
              cor_age65_male  = corD["age65", "male"],
              cor_age65_smoke = corD["age65", "smoke"],
              cor_htn_male    = corD["htn", "male"],
              cor_htn_smoke   = corD["htn", "smoke"],
              cor_male_smoke  = corD["male", "smoke"]
            )
          }
        }

        # Fit WLS for each trial and draw perturbed TE
        trial_fit_map    <- list()
        overall_draw_map <- list()

        for (tr in trials) {
          ov <- overall_tbl %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)
          TE_overall <- ov$TE_overall
          se_overall <- ov$se_overall

          TE_overall_draw <- TE_overall
          if (isTRUE(PERTURB_TE) && is.finite(TE_overall) &&
              is.finite(se_overall) && se_overall > 0) {
            TE_overall_draw <- rnorm(1, mean = TE_overall, sd = se_overall)
          }
          overall_draw_map[[tr]] <- list(
            TE_overall_draw = TE_overall_draw, se_overall = se_overall
          )

          rowp <- plan_k %>% dplyr::filter(trial == tr) %>% dplyr::slice(1)
          if (!isTRUE(rowp$is_large)) {
            trial_fit_map[[tr]] <- list(ok = FALSE, small = TRUE)
            next
          }

          sb <- Sigma_b_map[[tr]]
          if (is.null(sb)) {
            trial_fit_map[[tr]] <- list(ok = FALSE, small = FALSE,
                                        reason = "no_sigma_draw")
            next
          }

          agd_tr <- agd9_k %>% dplyr::filter(trial == tr)

          # Row-level TE perturbation
          if (isTRUE(PERTURB_TE)) {
            agd_tr <- agd_tr %>%
              dplyr::mutate(
                TE_draw = dplyr::if_else(
                  is.finite(TE) & is.finite(se) & se > 0,
                  rnorm(dplyr::n(), mean = TE, sd = se),
                  NA_real_
                )
              )
          } else {
            agd_tr <- agd_tr %>% dplyr::mutate(TE_draw = NA_real_)
          }

          # Mu vector: prefer overall means, fallback to external
          mu_over <- stats::setNames(
            as.numeric(unlist(ov[1, xcols], use.names = FALSE)), COVARS
          )
          mu_use <- mu_over
          if (any(!is.finite(mu_use))) {
            mu_use <- sb$mu
            names(mu_use) <- COVARS
          }

          imp_res <- tryCatch(
            blup_impute_agd9_draw(agd_tr, Sigma = sb$Sigma, mu_vec = mu_use,
                                  random = isTRUE(RANDOM_BLUP),
                                  clip = CLIP_RANGE),
            error = function(e) e
          )
          if (inherits(imp_res, "error")) {
            trial_fit_map[[tr]] <- list(
              ok = FALSE, small = FALSE,
              reason = paste0("impute_fail:", imp_res$message)
            )
            next
          }

          fit_obj <- tryCatch(
            fit_wls_trial_draw(imp_res$agd_imp, varset = varset,
                               perturb_te = isTRUE(PERTURB_TE)),
            error = function(e) list(ok = FALSE)
          )
          if (!isTRUE(fit_obj$ok)) {
            trial_fit_map[[tr]] <- list(ok = FALSE, small = FALSE,
                                        reason = "wls_fail")
            next
          }

          trial_fit_map[[tr]] <- list(ok = TRUE, small = FALSE, fit = fit_obj)

          # Save WLS coefficients
          if (isTRUE(SAVE_WLS_COEF_DRAWS)) {
            cf    <- coef(fit_obj$fit)
            se_cf <- suppressWarnings(sqrt(diag(fit_obj$V)))
            cc    <- cc + 1L
            coef_rows[[cc]] <- data.table::data.table(
              draw = b, trial = tr, term = names(cf),
              coef = as.numeric(cf), se = as.numeric(se_cf)
            )
          }

          # Save imputation diagnostics
          if (isTRUE(SAVE_IMPUTE_DIAG)) {
            i_imp <- i_imp + 1L
            imp_diag_rows[[i_imp]] <- data.table::data.table(
              draw = b, trial = tr,
              imp_age65  = imp_res$imp_n["age65"],
              imp_htn    = imp_res$imp_n["htn"],
              imp_male   = imp_res$imp_n["male"],
              imp_smoke  = imp_res$imp_n["smoke"],
              clip_age65 = imp_res$clip_n["age65"],
              clip_htn   = imp_res$clip_n["htn"],
              clip_male  = imp_res$clip_n["male"],
              clip_smoke = imp_res$clip_n["smoke"]
            )
          }
        }

        # ------ Per-degree: meta-analytic pooling -----------------------------
        for (deg in degrees) {
          tgt_row <- targets_by_delta %>%
            dplyr::filter(degree == deg) %>% dplyr::slice(1)
          p_star <- stats::setNames(
            as.numeric(unlist(tgt_row[1, COVARS], use.names = FALSE)), COVARS
          )
          tr_row  <- truth_by_delta %>%
            dplyr::filter(degree == deg) %>% dplyr::slice(1)
          TE_true <- tr_row$TE_true

          for (strategy in c("Hybrid_WLS_large", "Overall_only")) {
            te_vec <- numeric(length(trials))
            se_vec <- numeric(length(trials))
            names(te_vec) <- names(se_vec) <- trials

            for (tr in trials) {
              od  <- overall_draw_map[[tr]]
              TEo <- od$TE_overall_draw
              seo <- od$se_overall

              if (strategy == "Overall_only") {
                te_vec[tr] <- TEo
                se_vec[tr] <- seo
              } else {
                rowp <- plan_k %>%
                  dplyr::filter(trial == tr) %>% dplyr::slice(1)
                if (!isTRUE(rowp$is_large)) {
                  te_vec[tr] <- TEo
                  se_vec[tr] <- seo
                } else {
                  tf <- trial_fit_map[[tr]]
                  if (isTRUE(tf$ok)) {
                    pr <- predict_wls_at_target(tf$fit, p_star)
                    if (isTRUE(SAVE_WLS_PRED_DRAWS)) {
                      pp <- pp + 1L
                      pred_rows[[pp]] <- data.table::data.table(
                        draw = b, degree = deg, trial = tr,
                        te_pred = pr$te, se_pred = pr$se, varset = vn
                      )
                    }
                    if (is.finite(pr$te) && is.finite(pr$se) && pr$se > 0) {
                      te_vec[tr] <- pr$te
                      se_vec[tr] <- pr$se
                    } else {
                      te_vec[tr] <- TEo
                      se_vec[tr] <- seo
                    }
                  } else {
                    te_vec[tr] <- TEo
                    se_vec[tr] <- seo
                  }
                }
              }
            }

            # Save trial-level inputs
            if (isTRUE(SAVE_TRIAL_INPUTS)) {
              ti <- ti + 1L
              trial_input_rows[[ti]] <- data.table::data.table(
                draw = b, degree = deg, meta_strategy = strategy,
                trial = trials,
                te = as.numeric(te_vec[trials]),
                se = as.numeric(se_vec[trials])
              )
            }

            # Run meta models
            for (mdl in META_MODELS) {
              if (mdl == "FE") {
                obj <- fe_pool_closed(te_vec, se_vec)
              } else {
                obj <- metagen_safe_RE(
                  te_vec, se_vec, studlab = trials,
                  hakn_try = RE_HAKN_TRY,
                  method_tau_primary = RE_METHOD_TAU_PRIMARY,
                  method_tau_fallback = RE_METHOD_TAU_FALLBACK
                )
              }

              rr <- rr + 1L
              draw_rows[[rr]] <- data.table::data.table(
                rep = rep_id,
                scene_id = scene_id, pL = pL, pV = pV, K = K,
                pool = pool, varset = vn,
                degree = deg,
                meta_strategy = strategy,
                meta_model = mdl,
                draw = b,
                TE_true = TE_true,
                TE_est_draw = obj$TE,
                ok = isTRUE(obj$ok),
                tau2 = obj$tau2 %||% NA_real_,
                I2   = obj$I2   %||% NA_real_,
                Q    = obj$Q    %||% NA_real_
              )
            }
          }
        }
      }
      # end draw loop

      # ------ Summarize draws and save diagnostics ----------------------------
      draw_dt <- data.table::rbindlist(draw_rows, use.names = TRUE, fill = TRUE)

      sum_dt <- draw_dt[, {
        dd <- TE_est_draw[is.finite(TE_est_draw)]
        tt <- TE_true[which(is.finite(TE_true))[1]] %||% NA_real_
        if (length(dd) == 0 || !is.finite(tt)) {
          list(point_est = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
               width = NA_real_, cover = NA_integer_)
        } else {
          pe <- mean(dd)
          lo <- as.numeric(quantile(dd, 0.025, na.rm = TRUE, names = FALSE))
          hi <- as.numeric(quantile(dd, 0.975, na.rm = TRUE, names = FALSE))
          list(point_est = pe, ci_lo = lo, ci_hi = hi,
               width = hi - lo,
               cover = as.integer(lo <= tt && tt <= hi))
        }
      }, by = .(rep, scene_id, pL, pV, K, pool, varset,
                degree, meta_strategy, meta_model)]

      # Write diagnostic files
      save_case_diagnostics(vdir, draw_dt, sum_dt,
                            coef_rows, pred_rows, trial_input_rows,
                            imp_diag_rows, sigma_diag_rows)

      # Attach TE_true to summary
      sum_dt[, TE_true := draw_dt[
        , TE_true[which(is.finite(TE_true))[1]] %||% NA_real_,
        by = .(rep, scene_id, pL, pV, K, pool, varset,
               degree, meta_strategy, meta_model)
      ]$V1]

      out_rep_list[[length(out_rep_list) + 1]] <- sum_dt
    }
  }

  data.table::rbindlist(out_rep_list, use.names = TRUE, fill = TRUE)
}

# --- Helper: save base match diagnostics --------------------------------------
save_base_match_diagnostics <- function(base_match, pdir) {
  bm_sum <- data.table::rbindlist(lapply(names(base_match), function(tr) {
    bm <- base_match[[tr]]
    if (is.null(bm)) return(NULL)
    data.table::data.table(
      trial = tr,
      fallback = as.integer(isTRUE(bm$fallback)),
      fh_ok = as.integer(isTRUE(bm$fh_ok)),
      fh_action = bm$fh_action %||% NA_character_,
      fh_relax_mult = bm$fh_relax_mult %||% NA_real_,
      N_match = nrow(bm$sample),
      mu_age65 = bm$mu["age65"], mu_htn = bm$mu["htn"],
      mu_male  = bm$mu["male"],  mu_smoke = bm$mu["smoke"]
    )
  }), fill = TRUE)
  if (!is.null(bm_sum) && nrow(bm_sum) > 0) {
    data.table::fwrite(bm_sum, file.path(pdir, "base_match_summary.csv"))
  }

  bm_sig <- data.table::rbindlist(lapply(names(base_match), function(tr) {
    bm <- base_match[[tr]]
    if (is.null(bm) || is.null(bm$Sigma)) return(NULL)
    S  <- bm$Sigma
    dt <- data.table::as.data.table(as.table(S))
    if (ncol(dt) < 3L) return(NULL)
    data.table::setnames(dt, old = names(dt)[1:3],
                         new = c("var1", "var2", "value"))
    if (ncol(dt) > 3L) dt <- dt[, 1:3]
    dt[, trial := tr]
    dt
  }), fill = TRUE)
  if (!is.null(bm_sig) && nrow(bm_sig) > 0) {
    data.table::fwrite(bm_sig, file.path(pdir, "base_match_Sigma_long.csv"))
  }
}

# --- Helper: save per-case diagnostic files -----------------------------------
save_case_diagnostics <- function(vdir, draw_dt, sum_dt,
                                  coef_rows, pred_rows, trial_input_rows,
                                  imp_diag_rows, sigma_diag_rows) {
  if (isTRUE(SAVE_META_DRAWS)) {
    data.table::fwrite(draw_dt, file.path(vdir, "META_draws.csv"))
    data.table::fwrite(sum_dt,  file.path(vdir, "META_summary.csv"))
  }
  if (isTRUE(SAVE_WLS_COEF_DRAWS) && length(coef_rows)) {
    data.table::fwrite(
      data.table::rbindlist(coef_rows, fill = TRUE),
      file.path(vdir, "WLS_coef_draws.csv")
    )
  }
  if (isTRUE(SAVE_WLS_PRED_DRAWS) && length(pred_rows)) {
    data.table::fwrite(
      data.table::rbindlist(pred_rows, fill = TRUE),
      file.path(vdir, "WLS_pred_draws.csv")
    )
  }
  if (isTRUE(SAVE_TRIAL_INPUTS) && length(trial_input_rows)) {
    data.table::fwrite(
      data.table::rbindlist(trial_input_rows, fill = TRUE),
      file.path(vdir, "META_trial_inputs_draws.csv")
    )
  }
  if (isTRUE(SAVE_IMPUTE_DIAG) && length(imp_diag_rows)) {
    data.table::fwrite(
      data.table::rbindlist(imp_diag_rows, fill = TRUE),
      file.path(vdir, "IMPUTE_diag_draws.csv")
    )
  }
  if (isTRUE(SAVE_SIGMA_DIAG) && length(sigma_diag_rows)) {
    data.table::fwrite(
      data.table::rbindlist(sigma_diag_rows, fill = TRUE),
      file.path(vdir, "Sigma_diag_cor_draws.csv")
    )
  }
}
