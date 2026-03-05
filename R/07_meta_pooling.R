# =============================================================================
# 07_meta_pooling.R — Fixed-Effect and Random-Effects Meta-Analytic Pooling
# =============================================================================
# ADEMP Component: Methods (M) — Step 4 of EACC-Meta workflow
# =============================================================================

# --- Fixed-effect pooling (closed-form inverse-variance) ----------------------
fe_pool_closed <- function(te, se) {
  ok <- is.finite(te) & is.finite(se) & se > 0
  te <- te[ok]; se <- se[ok]
  k_used <- length(te)

  if (k_used < 2) {
    return(list(ok = FALSE, k_used = k_used,
                TE = NA_real_, lo = NA_real_, hi = NA_real_,
                tau2 = NA_real_, I2 = NA_real_, Q = NA_real_))
  }

  w      <- 1 / (se^2)
  TE_hat <- sum(w * te) / sum(w)
  se_hat <- sqrt(1 / sum(w))
  z      <- qnorm(0.975)
  lo     <- TE_hat - z * se_hat
  hi     <- TE_hat + z * se_hat

  # Cochran Q and I-squared for diagnostics
  Q  <- sum(w * (te - TE_hat)^2)
  df <- k_used - 1
  I2 <- if (is.finite(Q) && Q > 0) max(0, (Q - df) / Q) * 100 else 0

  list(ok = TRUE, k_used = k_used,
       TE = TE_hat, lo = lo, hi = hi,
       tau2 = 0, I2 = I2, Q = Q)
}

# --- Random-effects pooling via meta::metagen (robust wrapper) ----------------
metagen_safe_RE <- function(te, se, studlab,
                            hakn_try = TRUE,
                            method_tau_primary = "REML",
                            method_tau_fallback = "DL") {
  ok <- is.finite(te) & is.finite(se) & se > 0 & !is.na(studlab)
  te <- te[ok]; se <- se[ok]; studlab <- studlab[ok]
  k_used <- length(te)

  null_result <- list(ok = FALSE, k_used = k_used,
                      TE = NA_real_, lo = NA_real_, hi = NA_real_,
                      tau2 = NA_real_, I2 = NA_real_, Q = NA_real_)
  if (k_used < 2) return(null_result)

  fmls <- names(formals(meta::metagen))

  # Base arguments
  args_base <- list(
    TE = te, seTE = se, studlab = studlab,
    sm = "MD", prediction = FALSE
  )

  # Adapt to meta package version (common/random vs. comb.fixed/comb.random)
  if ("random" %in% fmls) {
    args_base$common <- FALSE
    args_base$random <- TRUE
  } else {
    args_base$comb.fixed  <- FALSE
    args_base$comb.random <- TRUE
  }
  if ("method.predict" %in% fmls) args_base$method.predict <- "V"

  call_one <- function(args) {
    args2 <- args[names(args) %in% fmls]
    do.call(meta::metagen, args2)
  }

  run_re <- function(method_tau, hakn_flag) {
    args <- args_base
    if ("method.tau" %in% fmls) args$method.tau <- method_tau
    if ("hakn" %in% fmls)       args$hakn <- isTRUE(hakn_flag)
    tryCatch(call_one(args), error = function(e) e)
  }

  # Try primary tau method with HAKN
  fit <- run_re(method_tau_primary, isTRUE(hakn_try))
  if (inherits(fit, "error") && isTRUE(hakn_try)) {
    fit2 <- run_re(method_tau_primary, FALSE)
    if (!inherits(fit2, "error")) fit <- fit2
  }

  # Fallback tau method
  if (inherits(fit, "error")) {
    fit3 <- run_re(method_tau_fallback, FALSE)
    if (inherits(fit3, "error")) return(null_result)
    fit <- fit3
  }

  # Extract random-effects estimates
  TE_r <- suppressWarnings(as.numeric(fit$TE.random))
  lo_r <- suppressWarnings(as.numeric(fit$lower.random))
  hi_r <- suppressWarnings(as.numeric(fit$upper.random))

  if (!is.finite(TE_r) || !is.finite(lo_r) || !is.finite(hi_r)) {
    tau2 <- if (!is.null(fit$tau2)) as.numeric(fit$tau2) else NA_real_
    I2   <- if (!is.null(fit$I2))   as.numeric(fit$I2)   else NA_real_
    Q    <- if (!is.null(fit$Q))    as.numeric(fit$Q)    else NA_real_
    return(list(ok = FALSE, k_used = k_used,
                TE = NA_real_, lo = NA_real_, hi = NA_real_,
                tau2 = tau2, I2 = I2, Q = Q))
  }

  tau2 <- if (!is.null(fit$tau2)) as.numeric(fit$tau2) else NA_real_
  I2   <- if (!is.null(fit$I2))   as.numeric(fit$I2)   else NA_real_
  Q    <- if (!is.null(fit$Q))    as.numeric(fit$Q)    else NA_real_

  list(ok = TRUE, k_used = k_used,
       TE = TE_r, lo = lo_r, hi = hi_r,
       tau2 = tau2, I2 = I2, Q = Q)
}
