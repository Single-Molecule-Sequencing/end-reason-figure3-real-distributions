# End-Reason Figure 3 — POD5-Derived Read-Length and Q-Score Distributions

<!-- LAB:DASHBOARD-BADGE BEGIN -->
📊 **[Live dashboard](https://single-molecule-sequencing.github.io/end-reason-figure3-real-distributions/)**
<!-- LAB:DASHBOARD-BADGE END -->


**Paper:** End-reason filtering in single-molecule sequencing of native DNA (Oxford Nanopore)  
**Figure:** Figure 3  
**What this shows:** Read-length kernel density estimates (KDEs) and Q-score KDEs,
stratified by POD5-derived `end_reason`, from the Cutting-resistant E Regular
sequencing run (`20250519_1041_MN48328_AYJ384_c3faa658`). Validates that
`signal_positive` reads best recapitulate the known physical molecule sizes, and
demonstrates that quality-score thresholds alone do not remove physically
discordant non-`signal_positive` populations.

---

## How this repo is organized

| Folder | What's inside |
|---|---|
| [`1_experiment/`](1_experiment/README.md) | Sequencing run details, dorado basecalling commands, pointers to raw and basecalled data |
| [`2_analysis/`](2_analysis/README.md) | Scripts, commands to run them, plain-English explanation of what the code does |
| [`3_results/`](3_results/README.md) | Polished figures and figure legends |
| [`HAPPY_PATH.md`](HAPPY_PATH.md) | Full end-to-end reproduce walkthrough (experiment → analysis → figure) |
| [`provenance/runs.jsonl`](provenance/runs.jsonl) | Append-only log of every run that produced an artifact in this repo |
| [`analysis.yaml`](analysis.yaml) | Lab system manifest (links this repo into the lab paper pipeline) |

## Quick start — reproduce the figure

> ⚠️ The raw POD5 input files live on Athey Lab Turbo (Great Lakes HPC). You
> need an active `umms-atheylab` allocation to reproduce from scratch.

```bash
# 1. Clone this repo
git clone https://github.com/Single-Molecule-Sequencing/end-reason-figure3-real-distributions
cd end-reason-figure3-real-distributions

# 2. Activate the lab Python environment
conda activate atheylab

# 3. Run the figure script
python 2_analysis/scripts/fig3_real_distributions.py
```

See [`HAPPY_PATH.md`](HAPPY_PATH.md) for the full annotated walkthrough.

## Current status

- **Draft date:** —
- **Status:** In progress — data sources confirmed; script being ported
- **Next step:** Run script on all three cohorts; review output
