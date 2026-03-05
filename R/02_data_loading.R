# =============================================================================
# 02_data_loading.R — Read and Standardize Real-World Data Pools
# =============================================================================
# ADEMP Component: Data-generating mechanisms (D) — external data sources
# =============================================================================

load_rwd_pools <- function(rwd_dir, match_pools, gen_pool) {
  rwd_map <- list()

  for (pool_name in c("CVD_all", "CVD_and_DM")) {
    base_name <- paste0("final_analysis_dataset_", pool_name)
    found <- any(file.exists(
      file.path(rwd_dir, paste0(base_name, c(".csv", ".xlsx", ".xls")))
    ))
    if (found) {
      rwd_map[[pool_name]] <- read_any(base_name, rwd_dir)
    }
  }

  if (!length(rwd_map)) {
    stop("No RWD files found in: ", rwd_dir)
  }

  # Standardize covariates
  rwd_std <- lapply(rwd_map, prep_covars_4)

  # Validate requested pools
  valid_pools <- intersect(match_pools, names(rwd_std))
  if (!length(valid_pools)) {
    stop("No valid match pools found. Available: ",
         paste(names(rwd_std), collapse = ", "))
  }

  if (!gen_pool %in% names(rwd_std)) {
    stop("Generation pool not found: ", gen_pool,
         ". Available: ", paste(names(rwd_std), collapse = ", "))
  }

  list(
    rwd_std     = rwd_std,
    match_pools = valid_pools,
    pool_df_gen = rwd_std[[gen_pool]]
  )
}
