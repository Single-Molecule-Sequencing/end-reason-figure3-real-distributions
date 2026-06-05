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

Figure 3 uses a **single sequencing run** from the SMS POP data collection:
the Cutting-resistant E (Regular) run on a MinION (run ID
`20250519_1041_MN48328_AYJ384_c3faa658`).

End reasons are read directly from **POD5 acquisition files** using the pod5
Python API. Read lengths and Q-scores are joined from the **basecaller
`sequencing_summary` file** by `read_id`. This self-contained data path (POD5 +
sequencing_summary) requires no additional annotation step.
