# Figure 3 Legend

<!-- TODO: Finalize the figure legend text for submission. The draft below
     reflects the fig.quality atom caption + legend — update once co-authors
     have reviewed the polished figure. -->

**Figure 3. POD5-derived read-length and Q-score distributions by end-reason class.**
Read-length kernel density estimates (KDEs) stratified by POD5-derived `end_reason`
for three proof-of-principle restriction-digest cohorts: single-species Eco53KI
(top row), RC Eco53KI+PvuII double-enzyme digest (middle row), and cutting-resistant
E control (bottom row). Left panels: read-length KDEs on log-scaled base-pair
coordinates; dashed vertical lines mark expected physical fragment-size peak centers
inferred from the `signal_positive` trough-bounded peak workflow. Middle panels:
per-read Q-score KDEs with a dashed reference line at Q10. Right panels: fraction
of reads per `end_reason` class falling within the expected physical fragment-size
window. End reasons are read directly from POD5 acquisition files; read lengths
and Q-scores are joined from basecaller `sequencing_summary` files by `read_id`.
`signal_positive` reads recapitulate the known physical molecule sizes; non-`signal_positive`
reads do not, and quality-score thresholds alone do not remove the physically
discordant populations.

---

## Notes for revision

- [ ] Confirm the exact cohort identities and display names
- [ ] Confirm expected fragment sizes (bp) for each cohort / enzyme combination
- [ ] Confirm Q10 as the reference threshold
- [ ] Add panel labels (A–I or row/column labels) if the journal requires them
- [ ] Verify that all end-reason class names match the canonical taxonomy
