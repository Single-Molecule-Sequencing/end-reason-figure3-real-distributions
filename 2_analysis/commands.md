# Commands — Exact Reproducibility Runs

## Environment

```bash
conda env create -f 2_analysis/scripts/environment.yaml -n end-reason-fig3 || true
conda activate end-reason-fig3
python -c "import pandas, numpy, matplotlib, seaborn; print('python deps OK')"
```

## A) Regenerate figure from deposited table (no HPC)

```bash
python 2_analysis/scripts/fig3_from_deposited_table.py \
  --table 3_results/tables/fig3_table2_single_experiment.csv \
  --out-prefix 3_results/figures/fig3_real_distributions_from_deposited_table
```

## B) Re-render Table 2 artifacts from summary stats

```bash
python 2_analysis/scripts/render_table2_single_experiment.py
```

## C) Recompute summary stats from raw data (HPC required)

```bash
# POD5 + custom BAM join
python 2_analysis/scripts/fig3_real_distributions_pod5_er_custom_bam.py \
  --pod5-dir /nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular/20250519_1041_MN48328_AYJ384_c3faa658/pod5 \
  --bam /nfs/turbo/umms-atheylab/hrli/Code/dorado-run/Output/20250519_1041_MN48328_AYJ384_c3faa658_sup_v5.2.0_trim1_10.bam \
  --out-dir 3_results/figures

# or BAM er:Z: path (dorado v1.3.1+)
python 2_analysis/scripts/fig3_real_distributions_bam_er_tag.py \
  --bam /path/to/basecalled_v1.3.1_or_later.bam \
  --out-dir 3_results/figures
```
