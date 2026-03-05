# =============================================================================
# 08_fh_bounds.R — Fréchet–Hoeffding Feasibility Bounds for Binary Covariates
# =============================================================================
# ADEMP Component: Methods (M) — Ensures correlation matrices are feasible
#                   under Bernoulli marginals (Nelsen et al., 2004)
# =============================================================================

# --- FH correlation bounds for a pair of Bernoulli variables ------------------
fh_bounds_cor <- function(p, q, eps = 1e-12) {
  p <- pmin(pmax(as.numeric(p), eps), 1 - eps)
  q <- pmin(pmax(as.numeric(q), eps), 1 - eps)
  sd <- sqrt(p * (1 - p) * q * (1 - q))
  if (!is.finite(sd) || sd <= 0) return(c(rmin = NA_real_, rmax = NA_real_))
  lo11 <- max(0, p + q - 1)
  hi11 <- min(p, q)
  rmin <- (lo11 - p * q) / sd
  rmax <- (hi11 - p * q) / sd
  c(rmin = rmin, rmax = rmax)
}

# --- Check all pairwise correlations against FH bounds ------------------------
fh_check_cor_matrix <- function(cor_mat, target_p, vars = COVARS, eps = 1e-8) {
  cor_mat  <- as.matrix(cor_mat)
  target_p <- as.numeric(target_p[vars])
  names(target_p) <- vars

  pairs <- utils::combn(vars, 2, simplify = FALSE)
  out <- lapply(pairs, function(pr) {
    v1 <- pr[1]; v2 <- pr[2]
    r  <- suppressWarnings(as.numeric(cor_mat[v1, v2]))
    b  <- fh_bounds_cor(target_p[v1], target_p[v2])
    ok <- is.finite(r) && is.finite(b["rmin"]) &&
          (r >= b["rmin"] - eps) && (r <= b["rmax"] + eps)
    data.table::data.table(
      var1 = v1, var2 = v2, r = r,
      rmin = b["rmin"], rmax = b["rmax"],
      ok = as.integer(ok)
    )
  })
  dt <- data.table::rbindlist(out)
  list(ok = all(dt$ok == 1L), report = dt)
}

# --- Clip correlation matrix to FH feasible bounds ----------------------------
fh_clip_cor_to_bounds <- function(cor_mat, target_p, vars = COVARS, eps = 1e-8) {
  C <- as.matrix(cor_mat)
  diag(C) <- 1
  target_p <- as.numeric(target_p[vars])
  names(target_p) <- vars

  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (j <= i) next
      v1 <- vars[i]; v2 <- vars[j]
      r  <- suppressWarnings(as.numeric(C[v1, v2]))
      b  <- fh_bounds_cor(target_p[v1], target_p[v2])
      if (is.finite(r) && is.finite(b["rmin"]) && is.finite(b["rmax"])) {
        r2 <- min(max(r, b["rmin"] + eps), b["rmax"] - eps)
        C[v1, v2] <- r2
        C[v2, v1] <- r2
      }
    }
  }
  C
}

# --- Rebuild covariance matrix from clipped correlation and target marginals --
sigma_from_cor_and_targetp <- function(cor_mat, target_p, vars = COVARS,
                                       eig_min = JITTER_EIG_MIN) {
  target_p <- as.numeric(target_p[vars])
  names(target_p) <- vars

  v <- target_p * (1 - target_p)
  v[!is.finite(v)] <- eig_min
  v <- pmax(v, eig_min)

  D <- diag(sqrt(v), length(vars))
  dimnames(D) <- list(vars, vars)

  C <- as.matrix(cor_mat)
  C[!is.finite(C)] <- 0
  diag(C) <- 1
  dimnames(C) <- list(vars, vars)

  S <- D %*% C %*% D
  dimnames(S) <- list(vars, vars)
  make_psd(S, eig_min = eig_min)
}

# --- Combined FH check-and-fix for a covariance matrix ------------------------
fh_fix_sigma_to_target <- function(Sigma, target_p, vars = COVARS,
                                   fh_eps = 1e-8, eig_min = JITTER_EIG_MIN) {
  Sigma <- as.matrix(Sigma)
  Sigma <- (Sigma + t(Sigma)) / 2
  Sigma <- make_psd(Sigma, eig_min = eig_min)

  cor0 <- suppressWarnings(stats::cov2cor(Sigma))
  if (any(!is.finite(cor0))) {
    cor0[!is.finite(cor0)] <- 0
    diag(cor0) <- 1
  }

  chk <- fh_check_cor_matrix(cor0, target_p = target_p,
                               vars = vars, eps = fh_eps)
  if (isTRUE(chk$ok)) {
    return(list(Sigma_fix = Sigma, fh_ok = TRUE,
                action = "pass", report = chk$report))
  }

  cor_clip  <- fh_clip_cor_to_bounds(cor0, target_p = target_p,
                                      vars = vars, eps = fh_eps)
  Sigma_fix <- sigma_from_cor_and_targetp(cor_clip, target_p = target_p,
                                           vars = vars, eig_min = eig_min)
  chk2 <- fh_check_cor_matrix(
    suppressWarnings(stats::cov2cor(Sigma_fix)),
    target_p = target_p, vars = vars, eps = fh_eps
  )

  list(Sigma_fix = Sigma_fix, fh_ok = isTRUE(chk2$ok),
       action = "clip_fix", report = chk2$report)
}
