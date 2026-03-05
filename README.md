# EACC-Meta Simulation Study



**Authors:** Lingyao Sun
 

---

## Study Overview (ADEMP Framework)

This simulation study was designed and reported following the **ADEMP framework** (Morris et al., 2019).

| ADEMP Component        | Description |
|------------------------|-------------|
| **A** – Aims           | Evaluate the statistical performance of EACC-Meta vs. conventional meta-analysis for estimating target-population treatment effects under incomplete aggregate data |
| **D** – Data-generating | Probit model with 4 binary covariates, 2 effect modifiers; 1,944 scenarios across transport intensity, covariate availability, and external data pools |
| **E** – Estimand       | Average treatment effect (ATE) in the pre-specified target population (probit scale) |
| **M** – Methods        | EACC-FE, EACC-RE vs. Conventional FE, Conventional RE (4 methods) |
| **P** – Performance    | Bias, RMSE, empirical coverage probability (95% CI), CI width; Monte Carlo SEs reported |

---

## Repository Structure

```
EACC-Meta-Simulation/
├── README.md                          # This file
├── main.R                             # Master run script (entry point)
├── R/
│   ├── 00_config.R                    # [A] All user-configurable parameters
│   ├── 01_utilities.R                 # General-purpose helper functions
│   ├── 02_data_loading.R             # [D] Read and standardize RWD pools
│   ├── 03_data_generation.R          # [D] DGP: probit model, trial simulation
│   ├── 04_raking_matching.R          # [M] Raking-based sampling and external matching
│   ├── 05_blup_imputation.R          # [M] BLUP-based covariate completion
│   ├── 06_wls_extrapolation.R        # [M] Within-trial WLS regression and extrapolation
│   ├── 07_meta_pooling.R             # [M] Fixed-effect and random-effects pooling
│   ├── 08_fh_bounds.R                # [M] Fréchet–Hoeffding feasibility bounds
│   ├── 09_fullboot_pipeline.R        # [M] Full-process bootstrap (FULLBOOT)
│   ├── 10_run_simulation.R           # Main simulation loop across reps and scenarios
│   └── 11_postprocessing.R           # [P] Rubin's rules and summary statistics
└── seeds/
    └── (seed manifests generated at runtime)
```

---

## Requirements

- **R** >= 4.1.0
- **Packages:** dplyr, tibble, readr, purrr, stringr, tidyr, MASS, survey, readxl, data.table, meta, future, future.apply, progressr

All packages are auto-installed on first run if missing.

---

## Data availability and runnable demo

**Restricted external data.** The full analysis in the manuscript uses external real-world cohort data to construct matching pools. These datasets are subject to data-provider agreements and cannot be shared publicly by the authors.

**Public code and workflow.** All simulation and analysis code in this repository is provided for transparency and reproducibility.

**Runnable demo without restricted data.** To enable users to run the workflow end-to-end without the restricted datasets, this repository includes a minimal worked example that automatically generates a synthetic external pool when the required real-world datasets are not available.

**Important note.** Results from the synthetic demo are illustrative and will not reproduce the manuscript’s numerical results. However, the demo reproduces the computational workflow and generates the same types of outputs (tables and figures).

## How to Run

1. **Prepare real-world data:** Place the external data files (`final_analysis_dataset_CVD_all.csv` and/or `final_analysis_dataset_CVD_and_DM.csv`) in a subdirectory named `real_world_data/` within the working directory.

2. **Configure parameters:** Edit `R/00_config.R` to adjust:
   - `N_REP`: Number of simulation repetitions (default: 200)
   - `B_BOOT`: Bootstrap draws per case (default: 200)
   - `WORK_DIR`: Working directory path
   - Other scenario parameters as needed

3. **Run the simulation:**
   ```r
   source("main.R")
   ```

4. **Outputs** are saved to `outputs_meta_reps/RUN_<timestamp>/`, including:
   - `REP_LEVEL_ALL.csv` – Rep-level estimates
   - `SUMMARY_200reps.csv` – Aggregated performance (Bias, RMSE, Coverage, etc.)
   - `SUMMARY_200reps_RUBIN.csv` – Rubin's rules combined estimates
   - Per-rep diagnostic files (WLS coefficients, imputation diagnostics, etc.)

---

## Reproducibility

- Random seeds are centrally managed via three independent RNG streams:
  - `SEED_GEN` for data generation
  - `SEED_EXT` for external matching
  - `SEED_BOOT` for bootstrap draws
- Trial-level intermediate outputs are exported by scenario.
- Seed manifests are saved per repetition in `seeds/`.

---

## Key References

- Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. *Statistics in Medicine*. 2019;38(11):2074-2102. doi:10.1002/sim.8086
- Siepe BS, Bartoš F, Morris TP, et al. Simulation studies for methodological research in psychology: A standardized template for planning, preregistration, and reporting. *Psychological Methods*. 2024. doi:10.1037/met0000695

---

## License

This code is provided for academic reproducibility purposes accompanying the submitted manuscript. Please cite the associated publication if reusing any part of this code.
