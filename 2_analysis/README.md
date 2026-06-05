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

The analysis reads end reasons from the POD5 acquisition files for the
Cutting-resistant E Regular run, joins read lengths and Q-scores from the
`sequencing_summary` file by `read_id`, then plots read-length KDEs and
Q-score KDEs stratified by `end_reason`, with dashed vertical lines marking
expected physical fragment-size peaks inferred from the `signal_positive`
trough-bounded peak workflow.
