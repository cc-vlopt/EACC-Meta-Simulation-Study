# =============================================================================
# 05_blup_imputation.R — BLUP-Based Covariate Completion
# =============================================================================
# ADEMP Component: Methods (M) — Step 2 of EACC-Meta workflow
#
# Uses externally anchored covariance structure to impute missing
# cross-classified subgroup proportions via Best Linear Unbiased Prediction.
# =============================================================================

blup_impute_agd9_draw <- function(agd_trial, Sigma, mu_vec,
                                  random = TRUE, clip = c(0, 1)) {
  xcols <- paste0("x_", COVARS)
  stopifnot(all(xcols %in% names(agd_trial)), "type" %in% names(agd_trial))

  Sigma <- make_psd(as.matrix(Sigma))
  rownames(Sigma) <- colnames(Sigma) <- COVARS

  mu_vec <- as.numeric(mu_vec)
  names(mu_vec) <- COVARS
  if (any(!is.finite(mu_vec))) stop("mu_vec contains NA or Inf")

  out <- agd_trial
  for (v in COVARS) out[[paste0("x_", v, "_raw")]] <- out[[paste0("x_", v)]]

  Xraw <- as.matrix(out[, xcols, drop = FALSE])
  colnames(Xraw) <- COVARS

  imp_n  <- stats::setNames(rep(0L, length(COVARS)), COVARS)
  clip_n <- stats::setNames(rep(0L, length(COVARS)), COVARS)

  # Impute one row using conditional normal distribution
  impute_row <- function(x_raw) {
    obs <- which(is.finite(x_raw))
    mis <- which(!is.finite(x_raw))
    x_imp <- x_raw

    if (length(mis) == 0) {
      return(list(x_final = x_imp,
                  n_imp = rep(0L, 4), n_clip = rep(0L, 4)))
    }

    if (length(obs) == 0) {
      mu_m <- mu_vec[mis]
      if (random) {
        S_mm <- make_psd(Sigma[mis, mis, drop = FALSE])
        mu_m <- as.numeric(MASS::mvrnorm(1, mu = mu_m, Sigma = S_mm))
      }
      mu0 <- mu_m
      mu_m <- pmin(pmax(mu_m, clip[1]), clip[2])
      n_clip <- rep(0L, 4); n_clip[mis] <- as.integer(mu_m != mu0)
      n_imp  <- rep(0L, 4); n_imp[mis]  <- 1L
      x_imp[mis] <- mu_m
      return(list(x_final = x_imp, n_imp = n_imp, n_clip = n_clip))
    }

    # Conditional distribution: X_mis | X_obs
    mu_o <- mu_vec[obs]; mu_m <- mu_vec[mis]
    x_o  <- x_raw[obs]
    S_oo <- Sigma[obs, obs, drop = FALSE]
    S_mo <- Sigma[mis, obs, drop = FALSE]
    S_mm <- Sigma[mis, mis, drop = FALSE]
    inv_Soo <- tryCatch(solve(S_oo), error = function(e) MASS::ginv(S_oo))

    mu_cond <- as.numeric(mu_m + S_mo %*% inv_Soo %*% (x_o - mu_o))

    if (random) {
      S_cond <- make_psd(S_mm - S_mo %*% inv_Soo %*% t(S_mo))
      x_m <- as.numeric(MASS::mvrnorm(1, mu = mu_cond, Sigma = S_cond))
    } else {
      x_m <- mu_cond
    }

    x0  <- x_m
    x_m <- pmin(pmax(x_m, clip[1]), clip[2])

    n_clip <- rep(0L, 4); n_clip[mis] <- as.integer(x_m != x0)
    n_imp  <- rep(0L, 4); n_imp[mis]  <- 1L
    x_imp[mis] <- x_m
    list(x_final = x_imp, n_imp = n_imp, n_clip = n_clip)
  }

  imp_list <- apply(Xraw, 1, impute_row)
  Xfinal <- do.call(rbind, lapply(imp_list, `[[`, "x_final"))
  n_imp  <- do.call(rbind, lapply(imp_list, `[[`, "n_imp"))
  n_clip <- do.call(rbind, lapply(imp_list, `[[`, "n_clip"))
  colnames(Xfinal) <- COVARS

  for (j in seq_along(COVARS)) {
    out[[paste0("x_", COVARS[j], "_final")]] <- Xfinal[, j]
    imp_n[COVARS[j]]  <- imp_n[COVARS[j]]  + sum(n_imp[, j],  na.rm = TRUE)
    clip_n[COVARS[j]] <- clip_n[COVARS[j]] + sum(n_clip[, j], na.rm = TRUE)
  }
  for (v in COVARS) out[[paste0("x_", v)]] <- out[[paste0("x_", v, "_final")]]

  list(agd_imp = out, imp_n = imp_n, clip_n = clip_n)
}
