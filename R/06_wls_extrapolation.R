# =============================================================================
# 06_wls_extrapolation.R — Within-Trial WLS Regression and Extrapolation
# =============================================================================
# ADEMP Component: Methods (M) — Step 3 of EACC-Meta workflow
# =============================================================================

fit_wls_trial_draw <- function(agd_imp, varset, perturb_te = FALSE) {
  dat <- agd_imp %>%
    dplyr::mutate(
      TE_use = dplyr::if_else(isTRUE(perturb_te) & is.finite(TE_draw),
                               TE_draw, TE),
      w = 1 / (se^2),
      w = dplyr::if_else(is.finite(w) & w > 0, w, NA_real_)
    ) %>%
    dplyr::filter(is.finite(TE_use), is.finite(se), is.finite(w))

  if (nrow(dat) < 2) return(list(ok = FALSE))

  rhs <- if (length(varset) == 0) {
    "1"
  } else {
    paste(paste0("x_", varset, "_final"), collapse = " + ")
  }
  fml <- as.formula(paste0("TE_use ~ ", rhs))
  fit <- tryCatch(lm(fml, data = dat, weights = w), error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE))

  V <- tryCatch(vcov(fit), error = function(e) NULL)
  if (is.null(V)) return(list(ok = FALSE))

  list(ok = TRUE, fit = fit, V = V)
}

predict_wls_at_target <- function(fit_obj, target_vec_named) {
  if (!isTRUE(fit_obj$ok)) return(list(te = NA_real_, se = NA_real_))
  fit   <- fit_obj$fit
  V     <- fit_obj$V
  terms <- names(coef(fit))

  xrow <- rep(0, length(terms))
  names(xrow) <- terms
  if ("(Intercept)" %in% terms) xrow["(Intercept)"] <- 1
  for (v in COVARS) {
    nm <- paste0("x_", v, "_final")
    if (nm %in% terms) xrow[nm] <- as.numeric(target_vec_named[[v]])
  }

  te_hat <- as.numeric(sum(xrow * coef(fit)))
  se_hat <- sqrt(as.numeric(t(xrow) %*% V %*% xrow))
  list(te = te_hat, se = se_hat)
}
