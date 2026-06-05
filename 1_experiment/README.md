# 1. Experiment

This section documents everything about **how the sequencing data was generated**:
the physical sequencing runs, the exact basecalling commands used to convert raw
signal to reads, and where the resulting data files live.

---

## Contents

| File / Folder | Description |
|---|---|
| [`run_details.md`](run_details.md) | Sequencing run metadata for each cohort |
| [`dorado_commands.md`](dorado_commands.md) | Exact dorado basecalling commands, in the order they were run |
| [`raw_data/`](raw_data/README.md) | Pointers to raw POD5 files on Turbo |
| [`basecalled_data/`](basecalled_data/README.md) | Pointers to sequencing_summary files on Turbo |

## What was the experiment?

Figure 3 draws from **three proof-of-principle restriction-digest cohorts**,
each sequenced on an Oxford Nanopore device:

| Cohort | Description | Expected fragment sizes |
|---|---|---|
| **Eco53KI** | Single-species digest | Known sharp peak(s) |
| **RC Eco53KI+PvuII** | Double-enzyme digest | Multiple defined peaks |
| **Cutting-resistant E** | Cutting-resistant control | Longer, uncut molecules dominate |

End reasons are read directly from **POD5 acquisition files** using the pod5
Python API. Read lengths and Q-scores are joined from **basecaller
`sequencing_summary` files** by `read_id`. This self-contained data path (POD5 +
sequencing_summary) requires no additional annotation step.
