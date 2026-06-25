# 3. Results

Final source artifacts for Figure 3 single-experiment analysis.

## Figure assets

- `3_results/figures/fig3_real_distributions.{png,pdf,svg}` (primary figure assets)
- `3_results/figures/fig3_real_distributions_from_deposited_table.{png,pdf,svg}` (table-derived reproducibility figure)
- `3_results/figures/fig3_real_distributions_from_deposited_table.lineage.json`

## Table assets

- `3_results/tables/fig3_table2_single_experiment.csv`
- `3_results/tables/fig3_table2_single_experiment.tex`
- `3_results/tables/fig3_table2_single_experiment.lineage.json`

Regenerate table assets:

```bash
python 2_analysis/scripts/render_table2_single_experiment.py
```

Regenerate deposited-table reproducibility figure:

```bash
python 2_analysis/scripts/fig3_from_deposited_table.py
```
