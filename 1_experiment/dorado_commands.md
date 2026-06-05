# Dorado Basecalling Commands

These are the **exact commands** used to basecall the raw signal data into reads
for each cohort. Copy-paste these to reproduce the `sequencing_summary` files
from the raw POD5 inputs.

---

## Prerequisites

```bash
# Activate the lab environment
conda activate atheylab

# Confirm dorado version (should match what was used originally)
dorado --version
```

---

## Cohort 1 — Eco53KI

<!-- TODO: Fill in the actual dorado command. -->

```bash
# [PLACEHOLDER — replace with the actual command]
dorado basecaller \
    <model> \
    <path/to/eco53ki/pod5/> \
    --reference <reference.fa> \
    > eco53ki_basecalled.bam
```

**Model used:** [?]  
**Reference:** [?]  
**Output:** see [`basecalled_data/`](basecalled_data/README.md)

---

## Cohort 2 — RC Eco53KI+PvuII

<!-- TODO: Fill in the actual dorado command. -->

```bash
# [PLACEHOLDER — replace with the actual command]
dorado basecaller \
    <model> \
    <path/to/rc_eco53ki_pvuii/pod5/> \
    --reference <reference.fa> \
    > rc_eco53ki_pvuii_basecalled.bam
```

**Model used:** [?]  
**Reference:** [?]

---

## Cohort 3 — Cutting-resistant E

<!-- TODO: Fill in the actual dorado command. -->

```bash
# [PLACEHOLDER — replace with the actual command]
dorado basecaller \
    <model> \
    <path/to/cutting_res_e/pod5/> \
    --reference <reference.fa> \
    > cutting_res_e_basecalled.bam
```

**Model used:** [?]  
**Reference:** [?]

---

## Notes

<!-- Any flags, filters, or post-processing steps applied after the main basecall commands -->
