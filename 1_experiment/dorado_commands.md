# Dorado Basecalling Commands

These are the **exact commands** used to basecall the raw signal data into reads.
Copy-paste these to reproduce the `sequencing_summary` file from the raw POD5 input.

---

## Prerequisites

```bash
# Activate the lab environment
conda activate atheylab

# Confirm dorado version (should match what was used originally)
dorado --version
```

---

## Basecalling command

<!-- TODO: Fill in the actual dorado command. -->

```bash
# [PLACEHOLDER — replace with the actual command]
dorado basecaller \
    <model> \
    /nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular/20250519_1041_MN48328_AYJ384_c3faa658/ \
    --reference <reference.fa> \
    > basecalled.bam
```

**Model used:** [?] (e.g., `dna_r10.4.1_e8.2_400bps_sup@v4.3.0`)  
**Reference:** [?]  
**Output:** see [`basecalled_data/`](basecalled_data/README.md)

---

## Notes

<!-- Any flags, filters, or post-processing steps applied after the main basecall command -->
