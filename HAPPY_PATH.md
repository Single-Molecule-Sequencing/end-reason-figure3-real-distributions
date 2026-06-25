# Happy Path — End-to-End Reproducibility

## Path A (recommended): reproduce from deposited table artifacts

No HPC access required.

```bash
git clone https://github.com/Single-Molecule-Sequencing/end-reason-figure3-real-distributions
cd end-reason-figure3-real-distributions
python 2_analysis/scripts/fig3_from_deposited_table.py \
  --table 3_results/tables/fig3_table2_single_experiment.csv \
  --out-prefix 3_results/figures/fig3_real_distributions_from_deposited_table
```

Expected outputs:

- `3_results/figures/fig3_real_distributions_from_deposited_table.{png,pdf,svg}`
- `3_results/figures/fig3_real_distributions_from_deposited_table.lineage.json`

## Path B: recompute summary from raw run data (Turbo/HPC)

1. Allocate Great Lakes compute and activate environment.
2. Run one of the figure scripts in `2_analysis/scripts/` (POD5+BAM or BAM `er:Z:` path).
3. Re-render Table 2:

```bash
python 2_analysis/scripts/render_table2_single_experiment.py
```

4. Optionally regenerate deposited-table reproducibility figure:

```bash
python 2_analysis/scripts/fig3_from_deposited_table.py
```

## Verification checklist

- `3_results/tables/fig3_table2_single_experiment.csv` exists
- `3_results/tables/fig3_table2_single_experiment.tex` exists
- `3_results/tables/fig3_table2_single_experiment.lineage.json` exists
- `3_results/figures/fig3_real_distributions.{png,pdf,svg}` exists
- `3_results/figures/fig3_real_distributions_from_deposited_table.{png,pdf,svg}` exists
