# 3. Results

This section contains the **final, polished outputs** of the analysis: the figure
as it will appear in the paper, and the figure legend text.

---

## Contents

| File / Folder | Description |
|---|---|
| [`figures/`](figures/) | Polished figure files (PDF for print, SVG for editing, PNG for web) |
| [`tables/`](tables/) | Source table artifacts harvested by the paper repo |
| [`figure_legends.md`](figure_legends.md) | Figure legend text as it will appear in the paper |

## What Figure 3 shows

Figure 3 validates that POD5-derived `end_reason` annotations reliably separate
physically concordant from discordant reads. Using the Cutting-resistant E
Regular run (`20250519_1041_MN48328_AYJ384_c3faa658`), `signal_positive` reads
produce a read-length distribution that matches the known physical fragment sizes,
while non-`signal_positive` reads do not. Crucially, Q-score filtering at Q10
does not rescue this: non-`signal_positive` reads survive quality filtering,
demonstrating that end-reason filtering is necessary and complementary to quality
filtering.

## Table 2 source artifact

The paper's Figure 3-associated Table 2 is sourced from this repo, not authored
directly in the paper repo. Regenerate it with:

```bash
python 2_analysis/scripts/render_table2_single_experiment.py
```

Outputs:

- `3_results/tables/fig3_table2_single_experiment.csv`
- `3_results/tables/fig3_table2_single_experiment.tex`

The table intentionally uses only the Cutting-resistant E Regular run shown in
Figure 3, merging the previous paper-level eight-run summary and per-end-reason
summary table into one single-experiment validation table.

## Draft status

- **Polished draft date:** —
- **Illustrator refinements:** Not yet applied
- **Pending:** Run script on all three cohorts; review raw output; Illustrator pass
