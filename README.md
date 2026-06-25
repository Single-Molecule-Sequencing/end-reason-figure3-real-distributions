# End-Reason Figure 3 — POD5-Derived Read-Length and Q-Score Distributions

<!-- LAB:DASHBOARD-BADGE BEGIN -->
📊 **[Live dashboard](https://single-molecule-sequencing.github.io/end-reason-figure3-real-distributions/)**
<!-- LAB:DASHBOARD-BADGE END -->

**Paper:** End-reason filtering in single-molecule sequencing of native DNA (Oxford Nanopore)  
**Figure:** Figure 3 (single-experiment view)  
**Run ID:** `20250519_1041_MN48328_AYJ384_c3faa658`

This repository is the source-of-record analysis package for Figure 3 real-distribution outputs and Table 2 single-experiment summary.

## Structure

- [`1_experiment/`](1_experiment/README.md): run metadata and data-location pointers
- [`2_analysis/`](2_analysis/README.md): reproducibility scripts, exact commands, environment
- [`3_results/`](3_results/README.md): figure and table artifacts (CSV/TEX/lineage)
- [`HAPPY_PATH.md`](HAPPY_PATH.md): end-to-end reproducibility workflows
- [`analysis.yaml`](analysis.yaml): analysis manifest
- [`docs/index.html`](docs/index.html): published analysis dashboard
- [`provenance/`](provenance/README.md): append-only run logs

## Reproduce quickly (from deposited table; no HPC)

```bash
git clone https://github.com/Single-Molecule-Sequencing/end-reason-figure3-real-distributions
cd end-reason-figure3-real-distributions

python 2_analysis/scripts/fig3_from_deposited_table.py \
  --table 3_results/tables/fig3_table2_single_experiment.csv \
  --out-prefix 3_results/figures/fig3_real_distributions_from_deposited_table
```

## Full recompute (from POD5/BAM on Turbo)

See [`2_analysis/commands.md`](2_analysis/commands.md) and [`HAPPY_PATH.md`](HAPPY_PATH.md).
