# Dorado Basecalling Commands

Exact original dorado invocation for this run is not present in deposited artifacts.
Known reproducible command pattern and known output path are below.

```bash
conda activate atheylab
dorado --version

# Template: fill model/reference to match the run environment.
dorado basecaller \
  <model> \
  /nfs/turbo/umms-atheylab/gregfar/SMS/SMS_POP_data/Single_Molecule_Seqeuncing_Cutting_Res_E/Regular/20250519_1041_MN48328_AYJ384_c3faa658/pod5 \
  --emit-sam \
  > /nfs/turbo/umms-atheylab/hrli/Code/dorado-run/Output/20250519_1041_MN48328_AYJ384_c3faa658_sup_v5.2.0_trim1_10.bam
```

See `1_experiment/unresolved.json` for provenance-gapped command parameters.
