# =============================================================================
# SESSION_INFO.R — Capture R Session Information for Reproducibility
# =============================================================================
# Run this script after completing the simulation to record
# the exact R version and package versions used.
#
# Usage: source("SESSION_INFO.R")
# =============================================================================

sink("SESSION_INFO.txt")
cat("EACC-Meta Simulation Study — Session Information\n")
cat("Generated:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
cat("=================================================\n\n")
sessionInfo()
cat("\n\n--- Installed Package Versions ---\n")
pkgs <- c("dplyr", "tibble", "readr", "purrr", "stringr", "tidyr",
           "MASS", "survey", "readxl", "data.table", "meta",
           "future", "future.apply", "progressr")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "NOT INSTALLED")
  cat(sprintf("  %-20s %s\n", p, v))
}
sink()
message("Session info saved to SESSION_INFO.txt")
