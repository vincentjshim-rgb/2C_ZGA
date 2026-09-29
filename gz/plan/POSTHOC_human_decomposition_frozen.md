# Post hoc plan — does the decomposition hold in human embryos? (frozen 2026-09-27, before any human clock score)

POST HOC. Written after every gate verdict and labelled as such throughout. Nothing here changes a gate verdict,
a reported interaction, or any conclusion already drawn in mouse; it asks one new question of a new species.

## Why

The manuscript's central result is that a scalar clock value in the mouse two-cell embryo is a small residual of
large opposing per-gene contributions, so the value can stay flat while those contributions are removed and
restored. That is a claim about how this clock reads early embryos. If it is a property of the clock rather than a
peculiarity of one mouse dataset, the same structure should appear in human embryos, where the clock's decrease has
been reported at a different interval. The Discussion already names this as the first open question; this plan turns
it into a result, whichever way it falls.

## Data

GSE36552 (Yan et al., 2013, *Nat Struct Mol Biol*), human preimplantation single cells, 124 samples. Deposited
per-sample expression files carry `Uniq_reads_num` — unique read counts per gene symbol — alongside RPKM; the counts
are used, the RPKM is not. Each file lists only the genes the sample expresses, so the matrix is the union of gene
symbols with absent entries set to zero.

Chosen over E-MTAB-3929 (Petropoulos et al., 2016) because that deposit was not reachable from this machine at the
time of writing; if it becomes reachable it may be added as a second human series, under this same plan and with the
same rules, and the attempt recorded either way.

**Exclusions, fixed here.** Embryonic stem cell samples (34) are not embryos and are dropped. Blastocyst lineage
samples are kept for the trajectory but are not part of the prespecified interval.

**Unit.** The embryo, not the cell: cells carrying the same embryo number in the deposited title are summed to an
embryo pseudobulk, as was done for GSE45719 in Gate 2b. Cells of one embryo are repeated measures, not replicates.

## Readout

The same clock and the same preprocessing as everywhere else in this study (EN_Chronoage_Multispecies_Multitissue
scaled-difference, Python port of tAge v1.1.0), scored against the oocytes of this dataset, which is the convention
used for GSE225056. Human is one of the clock's four training species, so no ortholog mapping is needed; gene
symbols are mapped to Entrez identifiers.

**Prespecified interval: 8-cell → morula.** This is the interval at which the human transcriptomic-age decrease was
reported (Zakar-Polyák et al., 2024). It is not a within-stage interval — it spans a cleavage — and that asymmetry
with the mouse two-cell window is a limitation of the comparison, not a result of it.

## Questions and reading rules, fixed before any human score

1. **Is there a decrease across the prespecified interval?** Report the clock change with a 95% percentile bootstrap
   interval over embryo pseudobulks within stage (2,000 replicates, seed 20260927). Read as present if the change is
   negative; the interval is reported whatever it does.

2. **Is it a residual of opposing contributions?** Decompose the interval with the exact identity used in Section 4.8
   (contribution = clock coefficient × change in the preprocessed feature; the terms sum to the reported difference).
   Define

       residual fraction = |D| / max(|sum of negative contributions|, |sum of positive contributions|)

   The mouse comparators, already computed, are 0.17 (GSE280522, within-two-cell) and 0.12 (GSE225056 mouse,
   within-two-cell). Reading rule: **the human decrease has the same structure if D < 0 and the residual fraction is
   at most 0.5**, i.e. the one-sided sums are at least twice the net. The fraction is reported whatever its value.

3. **How much of the clock does this dataset carry?** Report the number of the 1,839 weighted clock genes detected,
   the quantity that mouse GSE225056 turned out to report at only 912. A dataset carrying fewer than 300 of them is
   reported as too thin to decompose and question 2 is not read.

Outcome 1 or 2 failing is reported as it stands, in the manuscript, with the same prominence as a positive result. No
outcome licenses a change to any mouse conclusion.

## Outputs

    results/posthoc_human_decomposition.tsv    per-gene contributions across the prespecified interval
    results/posthoc_human_decomposition.md     the three readings above
    results/human_embryo_tage.tsv              clock value of every embryo pseudobulk
    figures/source_data/figH_human_*.tsv

## Seeds and software

Bootstrap seed 20260927. Software versions as recorded in Section 4.11.
