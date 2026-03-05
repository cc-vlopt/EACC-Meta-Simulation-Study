# =============================================================================
# 01_utilities.R — General-Purpose Helper Functions
# =============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_mkdir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

round_half_up <- function(x) floor(x + 0.5)

inv_probit <- function(z) pnorm(z)

# --- Truncated normal draw (single value) ------------------------------------
rtruncnorm_1 <- function(mu, sigma, a, b) {
  pa <- pnorm((a - mu) / sigma)
  pb <- pnorm((b - mu) / sigma)
  u  <- runif(1, pa, pb)
  mu + sigma * qnorm(u)
}

# --- Flexible file reader (CSV / XLSX / XLS) ---------------------------------
read_any <- function(base, rwd_dir) {
  cand <- file.path(rwd_dir, paste0(base, c(".csv", ".xlsx", ".xls")))
  if (file.exists(cand[1])) return(readr::read_csv(cand[1], show_col_types = FALSE))
  if (file.exists(cand[2])) return(readxl::read_xlsx(cand[2]))
  if (file.exists(cand[3])) return(readxl::read_xls(cand[3]))
  stop("File not found: ", base, " (.csv/.xlsx/.xls) in ", rwd_dir)
}

# --- Column name fuzzy matching -----------------------------------------------
col_has <- function(df, patt) {
  nm  <- names(df)
  hit <- nm[stringr::str_detect(tolower(nm), patt)]
  if (length(hit)) hit[1] else NA_character_
}

# --- Standardize binary variables ---------------------------------------------
std_bin <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    return(as.integer(tolower(x) %in% c("1", "y", "yes", "true", "male", "m", "current", "ever")))
  }
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    if (all(na.omit(x) %in% c(0, 1))) return(as.integer(x))
    if (is.finite(min(x, na.rm = TRUE)) && is.finite(max(x, na.rm = TRUE)) &&
        min(x, na.rm = TRUE) >= 0 && max(x, na.rm = TRUE) <= 1)
      return(as.integer(x >= 0.5))
    if (all(na.omit(x) %in% c(1, 2))) return(as.integer(x == 1))
    return(as.integer(x == 1))
  }
  stop("Cannot standardize variable to binary")
}

# --- Extract and standardize 4 binary covariates from a data frame -----------
prep_covars_4 <- function(df) {
  # Age >= 65
  nm_age65 <- col_has(df, "age\\s*65|age65|65\\+|65plus")
  nm_age   <- col_has(df, "^age$|age_year|agey")
  age65 <- if (!is.na(nm_age65)) {
    std_bin(df[[nm_age65]])
  } else if (!is.na(nm_age) && is.numeric(df[[nm_age]])) {
    as.integer(df[[nm_age]] >= 65)
  } else {
    rbinom(nrow(df), 1, 0.5)
  }

  # Hypertension
  nm_htn <- col_has(df, "htn|hypertens")
  htn <- if (!is.na(nm_htn)) std_bin(df[[nm_htn]]) else rbinom(nrow(df), 1, 0.5)

  # Male sex
  nm_sex <- col_has(df, "^male$|sex|gender")
  male <- if (!is.na(nm_sex)) {
    v <- df[[nm_sex]]
    if (is.character(v) || is.factor(v)) {
      as.integer(tolower(v) %in% c("male", "m", "1"))
    } else {
      std_bin(v)
    }
  } else {
    rbinom(nrow(df), 1, 0.5)
  }

  # Smoking status
  nm_smok <- col_has(df, "smok")
  smoke <- if (!is.na(nm_smok)) std_bin(df[[nm_smok]]) else rbinom(nrow(df), 1, 0.5)

  tibble::tibble(
    age65 = as.integer(age65),
    htn   = as.integer(htn),
    male  = as.integer(male),
    smoke = as.integer(smoke)
  )
}

# --- Probit-scale treatment effect from 2x2 counts ---------------------------
te_probit_from_counts <- function(r1, n1, r2, n2, cc = 0.5, eps = 1e-6) {
  if (any(!is.finite(c(r1, n1, r2, n2))) || any(c(n1, n2) <= 0)) {
    return(list(te = NA_real_, se = NA_real_))
  }
  p1 <- (r1 + cc) / (n1 + 2 * cc)
  p2 <- (r2 + cc) / (n2 + 2 * cc)
  p1 <- pmin(pmax(p1, eps), 1 - eps)
  p2 <- pmin(pmax(p2, eps), 1 - eps)
  z1 <- qnorm(p1); z2 <- qnorm(p2); te <- z1 - z2
  phi1 <- dnorm(z1); phi2 <- dnorm(z2)
  var1 <- p1 * (1 - p1) / (n1 * phi1^2)
  var2 <- p2 * (1 - p2) / (n2 * phi2^2)
  se   <- sqrt(var1 + var2)
  list(te = as.numeric(te), se = as.numeric(se))
}

# --- Enforce positive semi-definiteness ---------------------------------------
make_psd <- function(S, eig_min = JITTER_EIG_MIN) {
  S  <- as.matrix(S)
  dn <- dimnames(S)
  S  <- (S + t(S)) / 2
  S[!is.finite(S)] <- 0
  d <- diag(S)
  d[!is.finite(d)] <- 0
  d <- pmax(d, eig_min)
  diag(S) <- d
  eg   <- eigen(S, symmetric = TRUE)
  vals <- pmax(eg$values, eig_min)
  vals[!is.finite(vals)] <- eig_min
  S2 <- eg$vectors %*% diag(vals, length(vals)) %*% t(eg$vectors)
  S2 <- (S2 + t(S2)) / 2
  dimnames(S2) <- dn
  S2
}

# --- Output directory path helpers --------------------------------------------
fmt_scene_tag <- function(scene_id, pL, pV) {
  sprintf("scene%02d_pL%03d_pV%03d",
          as.integer(scene_id),
          as.integer(round(100 * pL)),
          as.integer(round(100 * pV)))
}

fmt_k_tag   <- function(K)   sprintf("K%02d", as.integer(K))
fmt_deg_tag <- function(deg)  sprintf("deg%02d", as.integer(round(10 * deg)))

rep_dir   <- function(rep_id)
  file.path(OUT_ROOT_RUN, sprintf("rep%03d", as.integer(rep_id)))
scene_dir <- function(rep_id, scene_id, pL, pV)
  file.path(rep_dir(rep_id), fmt_scene_tag(scene_id, pL, pV))
k_dir     <- function(rep_id, scene_id, pL, pV, K)
  file.path(scene_dir(rep_id, scene_id, pL, pV), fmt_k_tag(K))
case_dir  <- function(rep_id, scene_id, pL, pV, K, pool, varset)
  file.path(k_dir(rep_id, scene_id, pL, pV, K),
            sprintf("pool_%s", pool),
            sprintf("varset_%s", varset))
