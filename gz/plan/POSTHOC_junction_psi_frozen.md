# Post hoc plan — junction-level splicing in GSE280522 (frozen 2026-09-20, before any junction is counted)

POST HOC. Written after the Gate 4 verdict and after the cross-fitting, and labelled as such throughout. Nothing here
changes a gate verdict; the frozen Gate 4 result stands as reported.

## Why

Limitation 3 of the manuscript states that splicing was quantified from pseudoalignment-based transcript estimates on
annotated events, not from junction reads. That limitation can be tested directly: junction reads are an independent
measurement of the same quantity, obtained from a different aligner and a different statistic. This analysis asks
whether the Gate 4 readout survives a change of quantification method, and produces the read-level view (sashimi) that
transcript-level quantification cannot give.

## Data

GSE280522 only (P1): the dataset with three arms, the rescue arm and the cross-fitting, 23 libraries after the
prespecified quality control (the same 23 scored in Gate 3 and Gate 4). Raw reads are already on disk and were
verified against the archive MD5 checksums. The other two primary datasets are not re-aligned: P3 has 77 control
events, too few for a junction-level comparison, and P2 has no rescue arm.

## Reference and alignment

Ensembl GRCm39 release 112 primary assembly and its GTF — the same release and annotation used for the kallisto index
and for SUPPA2 event generation, so the event coordinates are identical. STAR 2.7 (the installed build), two-pass per
sample, default parameters except `--outSAMtype BAM SortedByCoordinate`, `--outSJfilterReads Unique` and
`--quantMode GeneCounts`. Splice junctions are taken from `SJ.out.tab` (uniquely mapped reads only).

## Readout

For every skipped-exon (SE) event in the Gate 4 filtered set of GSE280522, junction PSI is

    PSI_junction = (I1 + I2) / 2 / ((I1 + I2) / 2 + S)

where I1 and I2 are the uniquely mapped reads spanning the two inclusion junctions and S the reads spanning the
skipping junction, with the coordinates taken from the SUPPA2 event identifier. An event is scored in a library only
when I1 + I2 + S >= 10 uniquely mapped reads; a library with fewer than 20 scorable events is not used.

SE events only: the junction formula above is unambiguous for skipped exons, and the other six event classes need
either transcript-end evidence (AF, AL) or a different formula (RI, MX, A5, A3) that mixes junction and coverage
counts.

## Questions and reading rules, fixed before any junction is counted

1. **Agreement of the two quantifications.** Pearson and Spearman correlation between junction PSI and
   transcript-TPM PSI over all library-event pairs scored by both. Read as supporting the Gate 4 quantification if
   Spearman rho >= 0.6; as partially supporting if 0.4 <= rho < 0.6; as not supporting if rho < 0.4.
2. **Does the Gate 4 result reproduce?** The progression score of Section 4.10 is recomputed on junction PSI over the
   control-defined SE events, with the scale set on the control arm exactly as before. The interaction
   I = progression(A485) - progression(control) and the rescue difference R = progression(A485+DUX) -
   progression(A485) are reported with the same bootstrap (2,000 replicates, seed 20260921). Read as reproducing if
   I < 0 and R > 0; the magnitude is read only against the cross-fitted magnitude already reported, not against the
   in-sample one.
3. **Read-level illustration.** Sashimi coverage and junction counts for the three control-defined SE events with the
   largest |dPSI| in the control arm that also pass the junction depth rule in every arm, one panel per arm. Chosen by
   that rule, not by eye.

Outcome 1 or 2 failing is reported as it stands. No result here is used to change a gate verdict, an interaction
already reported, or the manuscript's conclusions; it is reported as a post-hoc test of the quantification method.

## Outputs

    results/posthoc_junction_psi.tsv        per library and event: junction counts, junction PSI, TPM PSI
    results/posthoc_junction_psi.md         the three readings above
    figures/source_data/figJ_junction_psi.tsv
    figures/source_data/figJ_sashimi.tsv    coverage and junction counts for the three illustrated events

## Seeds and software

Bootstrap seed 20260921. STAR 2.7 (system build), samtools, bedtools, Python 3.11 with the versions already recorded
in Section 4.11.
