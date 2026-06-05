# Happy Path — End-to-End Reproduce Guide

This document is your **single source of truth** for reproducing Figure 3 from
scratch. Follow these steps in order. Each step links to the relevant files.

> ⚠️ **HPC required for steps 1–3.** The POD5 files live on Turbo.
> Steps 4–5 can be done anywhere once you have the output files.

---

## Step 0 — Prerequisites

1. Access to Great Lakes HPC with an active `umms-atheylab` allocation
2. `conda activate atheylab` (or equivalent with the packages in [`2_analysis/overview.md`](2_analysis/overview.md))
3. This repo cloned locally (or on Turbo)

---

## Step 1 — Confirm the raw data is accessible

```bash
# Should list the POD5 files for the run
ls /nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular/20250519_1041_MN48328_AYJ384_c3faa658/
```

📂 More details: [`1_experiment/raw_data/README.md`](1_experiment/raw_data/README.md)

---

## Step 2 — Confirm basecalled data is accessible

```bash
# Should list sequencing_summary files for each cohort
ls <path/to/basecalled/>   # see 1_experiment/basecalled_data/README.md
```

📂 More details: [`1_experiment/basecalled_data/README.md`](1_experiment/basecalled_data/README.md)  
🔧 How it was produced: [`1_experiment/dorado_commands.md`](1_experiment/dorado_commands.md)

---

## Step 3 — Run the analysis script

```bash
conda activate atheylab

python 2_analysis/scripts/fig3_real_distributions.py \
    --run-dir /nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular/20250519_1041_MN48328_AYJ384_c3faa658 \
    --out-dir 3_results/figures/raw_output \
    --peak-bp <EXPECTED_BP>   # fill in expected fragment size(s)
```

📖 What the script does (plain English): [`2_analysis/overview.md`](2_analysis/overview.md)  
🔧 Full command details: [`2_analysis/commands.md`](2_analysis/commands.md)

**Expected output:** KDE figure files written to `3_results/figures/raw_output/`

---

## Step 4 — Apply Illustrator refinements (manual)

Open the raw output in Adobe Illustrator, apply typographic and layout
refinements, export as:

```
3_results/figures/Figure_3_final.pdf
3_results/figures/Figure_3_final.svg
3_results/figures/Figure_3_final@4x.png
```

---

## Step 5 — Verify and record

```bash
# Record the provenance run (stamps the figure atom with lineage)
lab-analysis record-run \
    --figure fig3_real_distributions \
    --command "python 2_analysis/scripts/fig3_real_distributions.py" \
    --output 3_results/figures/Figure_3_final.pdf

# Check current maturity status
lab-analysis status
```

📋 Results: [`3_results/README.md`](3_results/README.md)  
📖 Figure legend: [`3_results/figure_legends.md`](3_results/figure_legends.md)

---

## Provenance trail

Every artifact-producing run is logged in [`provenance/runs.jsonl`](provenance/runs.jsonl).
This file is committed to git and is the **durable source of truth** for this
figure's lineage. Inspect it with:

```bash
lab-analysis status --provenance
```
