# Analysis Overview

This repo supports two reproducibility paths:

1. **Deposited-table reproducibility (default):**
   regenerate figure assets from `3_results/tables/fig3_table2_single_experiment.csv`.
2. **Raw-data recomputation (HPC):**
   regenerate `fig3_summary_stats.csv` and upstream outputs from POD5/BAM.

The deposited-table path is deterministic and intended for strict reproducibility
without requiring Turbo/HPC access.
