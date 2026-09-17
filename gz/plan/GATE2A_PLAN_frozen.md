# Gate 2a plan — is the transcriptomic-age drop a compositional artefact of maternal clearance / zygotic activation?
(frozen 2026-09-14, before any clock score under the variants below was computed)

## Status
Robustness test of the Gate 1 / Gate 1-alt observation (itself post hoc for pig/rabbit). The test is prespecified
here; its outcome decides whether the "reset" reading survives.

## Threat
tAge normalises the whole transcriptome (RLE → log → per-sample scaling → YuGene → reference centring). When
maternal mRNA that makes up a large share of the library is cleared (mouse −36%, cow −28% share in the ZGA
interval), every remaining gene's relative level shifts, which could move the clock without any change in
age-related programmes.

## Data
GSE225056 5' gene counts, mouse and cow (the two species where the clock drop coincided with the data-driven ZGA /
clearance interval), plus pig and rabbit as secondary. Same QC and identifier mapping as Gate 1 (+ addendum 1).
Clock: EN_Chronoage_Multispecies_Multitissue_scaleddiff, oocyte reference, Python port.

## Gene sets (recomputed per species WITHOUT excluding clock genes, same thresholds as Gate 1-alt)
Maternal genes: CPM >= 10 in >= 80% of oocytes and zygotes, and median CPM in some later stage <= 1/8 of the oocyte median.
Zygotic genes: CPM < 1 in >= 80% of oocytes and zygotes, and CPM >= 10 in >= 50% of libraries of some later stage.
Dynamic set D = maternal ∪ zygotic.

## Variants
- V0 original (reproduces Gate 1).
- V1 remove D genes that are NOT clock features before preprocessing (normalisation no longer sees the bulk of
  maternal/zygotic mass; clock features untouched).
- V2 remove all D genes before preprocessing, including D genes that are clock features (those features become
  missing and are mean-imputed by the model's own imputer) — removes both the compositional and the direct
  contribution of maternal/zygotic genes.
- V3 composition-only simulation (positive control for the artefact): take the real oocyte libraries of the species,
  multiply the counts of maternal genes by the species' observed ratio (mean maternal share in the later stage of
  the data-driven ZGA interval / mean share in the earlier stage), leave every other gene unchanged, and score these
  pseudo-libraries together with the real oocytes (oocyte reference). Δ_sim = mean(pseudo) − mean(oocytes).

## Primary quantity
Clock Δ in the data-driven ZGA interval (Gate 1-alt: mouse Early-2-cell → Late-2-cell; cow 8-cell → 16-cell), with 95%
bootstrap CI (embryos within stage, 2000 replicates, seed 20260916). Retention = Δ_V2 / Δ_V0.

## Decision rules (mouse and cow)
- RESET-ROBUST: in both species Δ_V2 < 0 with CI excluding 0 AND retention >= 0.5, AND |Δ_sim| < 0.5 × |Δ_V0|.
- COMPOSITIONAL: in both species either the V2 CI includes 0 or retention < 0.25, OR Δ_sim reproduces >= 0.75 of Δ_V0 in both.
- INTERMEDIATE: anything else (reported per species).
Secondary: V1 retention; all intervals and pig/rabbit reported descriptively.

## Outputs
results/gate2a_intervals.tsv, results/gate2a_simulation.tsv, results/gate2a_decision.md
