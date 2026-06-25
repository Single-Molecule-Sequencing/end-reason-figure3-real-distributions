# Figure 3 Legend


**Figure 3. POD5-derived read-length and Q-score distributions by end-reason class.**
Read-length kernel density estimates (KDEs) and Q-score KDEs stratified by
POD5-derived `end_reason` for the Cutting-resistant E Regular run
(`20250519_1041_MN48328_AYJ384_c3faa658`). Left panel: read-length KDE on
log-scaled base-pair coordinates; dashed vertical lines mark expected physical
fragment-size peak centers inferred from the `signal_positive` trough-bounded
peak workflow. Middle panel: per-read Q-score KDE with a dashed reference line
at Q10. Right panel: fraction of reads per `end_reason` class falling within the
expected physical fragment-size window. End reasons are read directly from POD5
acquisition files; read lengths and Q-scores are joined from the basecaller
`sequencing_summary` file by `read_id`. `signal_positive` reads recapitulate the
known physical molecule sizes; non-`signal_positive`
reads do not, and quality-score thresholds alone do not remove the physically
discordant populations.
