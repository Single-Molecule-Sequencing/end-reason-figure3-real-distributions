# Commands — How to Run the Analysis

Run these in order to reproduce the figure from scratch.

> ⚠️ **HPC required.** The raw POD5 and BAM inputs live on Turbo
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
```

Choose the right import check for your script (see Step 1):

```bash
# POD5 script (dorado < v1.3.1)
python -c "import pod5, pysam, pandas, numpy, matplotlib, seaborn; print('All imports OK')"

# BAM er:Z: script (dorado v1.3.1+)
python -c "import pysam, pandas, numpy, matplotlib, seaborn; print('All imports OK')"
```

---

## 1. Run the figure script

Choose the script that matches the dorado version used to basecall the experiment.

### Option A — dorado < v1.3.1 (end reasons from POD5)

Reads end reasons from POD5 acquisition files. Requires `pod5/` and `bam_pass/`
subdirectories inside `--run-dir`.

```bash
cd /path/to/end-reason-figure3-real-distributions

python 2_analysis/scripts/fig3_real_distributions_pod5_end_reasons.py \
    --run-dir /nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular/20250519_1041_MN48328_AYJ384_c3faa658 \
    --out-dir 3_results/figures/raw_output \
    --peak-bp 5500   # adjust to expected fragment size(s)
```

### Option B — dorado v1.3.1+ (end reasons from BAM er:Z: tag)

Reads end reasons directly from the BAM `er:Z:` aux tag. Only requires
`bam_pass/` inside `--run-dir` — no POD5 files needed.

> Note: the `er:Z:` tag was first written by dorado in **v1.3.1** (DOR-1307
> backport, 2026-01-12). It is absent in v1.3.0 and earlier.

```bash
cd /path/to/end-reason-figure3-real-distributions

python 2_analysis/scripts/fig3_real_distributions_bam_er_tag.py \
    --run-dir /nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular/20250519_1041_MN48328_AYJ384_c3faa658 \
    --out-dir 3_results/figures/raw_output \
    --peak-bp 5500   # adjust to expected fragment size(s)
```

Both scripts produce identical output:
**`3_results/figures/raw_output/fig3_real_distributions.{pdf,png}`** and `lineage.json`  
**Runtime:** ~5–15 minutes depending on BAM size and allocation CPUs

---

## 2. Record the provenance run

After executing, stamp the run so this figure has a traceable lineage:

```bash
lab-analysis record-run \
    --figure fig3_real_distributions \
    --command "python 2_analysis/scripts/fig3_real_distributions_bam_er_tag.py" \
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
