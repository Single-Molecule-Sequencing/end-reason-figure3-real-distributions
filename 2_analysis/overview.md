# Analysis Overview — Plain English

This document explains **what the analysis code does and why**, without requiring
you to read the script. Think of it as the methods section in plain language.

---

## Goal

Produce Figure 3: read-length KDE and Q-score KDE panels stratified by
`end_reason` for the Cutting-resistant E Regular run
(`20250519_1041_MN48328_AYJ384_c3faa658`):

| Panel | Content |
|---|---|
| Left | Read-length KDE by `end_reason` (log-scaled bp x-axis; dashed vertical lines for expected physical peak centers) |
| Middle | Q-score KDE by `end_reason` (Q10 threshold reference line) |
| Right | Expected-physical-window percentage table per `end_reason` class |

The key result: `signal_positive` reads have a read-length distribution that
recapitulates the known physical fragment sizes. Non-`signal_positive` reads do
not. A Q10 quality filter does not fix this — discordant reads survive quality
filtering, which is the main motivation for end-reason filtering.

---

## Step-by-step logic

### Step 1 — Read end reasons from POD5

The `pod5` Python API exposes an `end_reason` field per read in the POD5
acquisition file. We read these directly from the Cutting-resistant E Regular
run directory using `pod5.Reader` → `read_table`. This gives us the
ground-truth sequencer-assigned end reason without any downstream annotation.

### Step 2 — Join read length and Q-score from sequencing_summary

The `sequencing_summary.txt` produced by dorado contains per-read metrics
including `sequence_length_template` and `mean_qscore_template`. We join this
table to the POD5 end-reason table on `read_id` using pandas.

### Step 3 — Compute expected physical peak centers

For each cohort, expected fragment sizes (in bp) are known from the restriction
enzyme cut-site positions. We infer the peak center from the `signal_positive`
read-length distribution using a trough-bounded peak workflow (finding the
dominant peak in the KDE of `signal_positive` reads).

### Step 4 — Plot read-length KDEs (left column)

For each cohort, plot one KDE curve per `end_reason` class on a log-scaled
base-pair x-axis. Overlay dashed vertical lines at the expected physical peak
centers.

### Step 5 — Plot Q-score KDEs (middle column)

For each cohort, plot one KDE curve per `end_reason` class on a linear Q-score
x-axis. Overlay a dashed vertical reference line at Q10.

### Step 6 — Compute physical-window percentages (right column)

For each cohort and each `end_reason` class, compute the fraction of reads
falling within the expected physical fragment-size window. Display as a small
table or bar chart per cohort.

---

## Dependencies

| Library | Purpose |
|---|---|
| `pod5` | Read `end_reason` from POD5 acquisition files |
| `pandas` | `read_id` join between POD5 and sequencing_summary |
| `numpy` | KDE computation and trough-bounded peak workflow |
| `scipy` | Gaussian KDE |
| `matplotlib` | Plotting the KDE panels and summary table |
