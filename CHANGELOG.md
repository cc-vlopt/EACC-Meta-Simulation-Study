# CHANGELOG — Original → Submission-Ready Code

## Summary of Refactoring

The original single-file script (`重复显示进程.R`, 1,977 lines) was refactored into a modular, well-documented project structure suitable for public repository submission.

---

## Critical Bug Fixes

### 1. Duplicate `match_external_base_4()` definition (SEVERE)
**Original:** The function was defined twice — a basic version at line 513 and an FH-enhanced version at line 1872. R would only use the second definition, but the second definition was placed **after** the main simulation loop, meaning the basic (non-FH) version was the one actually executed during simulation.

**Fix:** Merged into a single, unified version in `R/04_raking_matching.R` that includes FH bounds checking, tolerance relaxation, and clip-fix — ensuring the intended FH logic is always active.

### 2. Orphan code block at lines 1964–1976 (SEVERE)
**Original:** Lines 1964–1976 contained bare statements (`rowc <- cov_obs_k %>% ...`, `corS <- ...`, etc.) outside any function. These would execute at source time against undefined variables, causing errors or silent failures.

**Fix:** This FH-check-per-draw logic was integrated into the bootstrap draw loop inside `R/09_fullboot_pipeline.R`, where `cov_obs_k` and `Sig_b` are properly in scope.

### 3. `gen_one_trial_AB_fast()` used `<<-` for `agd_list` (FRAGILE)
**Original:** The `add_cov_rows()` inner function used `<<-` to modify the parent scope's `agd_list`, which is fragile and can cause subtle bugs in parallel contexts.

**Fix:** Restructured `build_agd9()` in `R/03_data_generation.R` as a self-contained function using standard `list()` accumulation — no `<<-` needed.

---

## Structural Changes

| Aspect | Original | Refactored |
|--------|----------|------------|
| Files | 1 monolithic file (1,977 lines) | 14 files across clear modules |
| ADEMP mapping | Not explicit | Each file header states its ADEMP component |
| Entry point | `source("重复显示进程.R")` | `source("main.R")` |
| Parameters | Mixed into code body | Isolated in `R/00_config.R` |
| FH bounds | Appended after main loop | Properly integrated in pipeline |
| Working directory | Hardcoded `C:/Users/vlpot/Desktop/meta插补` | Uses `getwd()`, configurable in `00_config.R` |
| RWD folder name | `真实世界数据` (Chinese) | `real_world_data` (English) |
| Language | Mixed Chinese/English comments | All English (submission-ready) |

---

## File Mapping (Original Line Ranges → New Files)

| Original Lines | Content | New File |
|---------------|---------|----------|
| 1–128 | Parameters, packages, parallel setup | `R/00_config.R` + `main.R` |
| 130–145 | Path helpers, format tags | `R/01_utilities.R` |
| 148–225 | Utility functions | `R/01_utilities.R` |
| 228–251 | RWD loading | `R/02_data_loading.R` |
| 254–500 | DGP, trial generation, targets | `R/03_data_generation.R` |
| 270–326 | Raking sampling | `R/04_raking_matching.R` |
| 504–572 | `choose_N_match`, `match_external_base_4`, `make_psd` | `R/04_raking_matching.R` + `R/01_utilities.R` |
| 575–659 | BLUP imputation | `R/05_blup_imputation.R` |
| 661–696 | WLS fitting and prediction | `R/06_wls_extrapolation.R` |
| 698–834 | FE and RE pooling | `R/07_meta_pooling.R` |
| 838–1292 | FULLBOOT pipeline | `R/09_fullboot_pipeline.R` |
| 1294–1456 | `run_one_rep` | `R/10_run_simulation.R` |
| 1458–1548 | Main loop and summary | `main.R` |
| 1550–1757 | Rubin post-processing | `R/11_postprocessing.R` |
| 1765–1868 | FH bounds functions | `R/08_fh_bounds.R` |
| 1872–1962 | Duplicate `match_external_base_4` | Merged into `R/04_raking_matching.R` |
| 1964–1976 | Orphan FH check code | Integrated into `R/09_fullboot_pipeline.R` |

---

## Removed Items

- **Debug prints:** All `cat(sprintf("Trial %s: agd_list 长度 = %d ...", ...))` and `print(table(...))` removed
- **Redundant `library()` calls:** Second `library(data.table)` at line 1551 removed
- **Duplicate `progressr::handlers()`:** Called twice (lines 103 and 112); consolidated to once in `main.R`
- **Warning emoji:** `⚠️` characters in `warning()` calls replaced with plain text
- **`SAVE_IPD` default:** Remains `FALSE` (200 reps × 15 trials × ~10K rows = massive files)

---

## Added Items

- `README.md` — Project overview with ADEMP table, requirements, usage instructions
- `.gitignore` — Excludes outputs/ and real_world_data/
- `SESSION_INFO.R` — Captures exact R and package versions for reproducibility
- `CHANGELOG.md` — This file
- ADEMP component labels in every file header
- FH bounds check inside bootstrap draw loop (was orphan code, now properly placed)
- `fh_ok`, `fh_action`, `fh_relax_mult` fields in base match summary output

---

## Action Items Before Submission

1. **Rename RWD folder:** `真实世界数据/` → `real_world_data/`
2. **Test run:** Execute `main.R` with `N_REP <- 2L` and `B_BOOT <- 5L` to verify
3. **Run `SESSION_INFO.R`** and commit the output
4. **Upload to GitHub** and create a Zenodo DOI
5. **Update `00_config.R`:** Set final `N_REP` and `B_BOOT` values
