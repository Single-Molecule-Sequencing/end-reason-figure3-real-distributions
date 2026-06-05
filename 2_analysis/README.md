# 2. Analysis

This section contains the **code that turns raw data into figures**, along with
plain-English explanations of what the code does and how to run it.

---

## Contents

| File / Folder | Description |
|---|---|
| [`overview.md`](overview.md) | Plain-English walkthrough of the analysis logic (no code) |
| [`commands.md`](commands.md) | Exact commands to run each script, in order |
| [`scripts/`](scripts/) | The actual code (`fig3_real_distributions.py`) |

## One-line summary

The analysis reads end reasons from POD5 acquisition files, joins read lengths
and Q-scores from `sequencing_summary` files by `read_id`, then plots read-length
KDEs and Q-score KDEs stratified by `end_reason` for each of the three
proof-of-principle cohorts, with dashed vertical lines marking expected physical
fragment-size peaks inferred from the `signal_positive` trough-bounded peak workflow.
