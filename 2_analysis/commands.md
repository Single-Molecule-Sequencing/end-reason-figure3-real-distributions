# Commands — How to Run the Analysis

Run these in order to reproduce the figure from scratch.

> ⚠️ **HPC required.** The raw POD5 and sequencing_summary inputs live on Turbo
> and are only accessible from a Great Lakes compute node. Run from a GL allocation:
> `srun --partition=standard --cpus-per-task=4 --mem=32G --pty bash`

---

## 0. Prerequisites

```bash
# SSH to Great Lakes
ssh <uniqname>@greatlakes.arc-ts.umich.edu

# Get an interactive allocation (if not already in one)
srun --partition=standard --cpus-per-task=4 --mem=32G --pty bash

# Activate the lab environment
conda activate atheylab

# Verify key packages are importable
python -c "import pod5, pandas, numpy, scipy, matplotlib; print('All imports OK')"
```

---

## 1. Run the figure script

This reads end reasons from the POD5 files for run
`20250519_1041_MN48328_AYJ384_c3faa658`, joins the sequencing_summary data, and
produces read-length KDE + Q-score KDE panels.

```bash
cd /path/to/end-reason-figure3-real-distributions

python 2_analysis/scripts/fig3_real_distributions.py
```

**Output:** KDE figure files written to `3_results/figures/raw_output/`  
**Runtime:** ~5–10 minutes per cohort depending on POD5 size and allocation CPUs

---

## 2. Record the provenance run

After executing, stamp the run so this figure has a traceable lineage:

```bash
lab-analysis record-run \
    --figure fig3_real_distributions \
    --command "python 2_analysis/scripts/fig3_real_distributions.py" \
    --output 3_results/figures/Figure_3_final.pdf \
    --output 3_results/figures/Figure_3_final.png
```

---

## 3. (Manual) Illustrator refinements

Open `3_results/figures/raw_output/` in Adobe Illustrator, apply final
typographic/layout refinements, and export as:
- `3_results/figures/Figure_3_final.pdf`
- `3_results/figures/Figure_3_final.svg`
- `3_results/figures/Figure_3_final@4x.png`

Commit the exported files.

---

## 4. Verify the happy path

```bash
# From repo root
lab-analysis verify-happy-path
```

See [`HAPPY_PATH.md`](../HAPPY_PATH.md) for the full end-to-end narrative.
