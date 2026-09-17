# Gate 2b plan — independent replication in mouse: does transcriptomic age drop within the 2-cell stage?
(frozen 2026-09-14, before any clock score on GSE45719 or GSE66582)

## Hypothesis (from GSE225056 mouse, Gate 1)
The largest pre-blastocyst decrease in transcriptomic age occurs within the 2-cell stage (early → late 2-cell), the
interval of major ZGA and bulk maternal clearance, with no cell division.

## Dataset R1 (primary): GSE45719 (Deng et al. 2014, Science), single-cell Smart-seq, CAST × B6 F1
Per-cell files `*_expression.txt.gz` (gene symbol, reads). Unit = embryo: reads of all cells from the same embryo are
summed (pseudobulk). Embryo ID = the token before the final "-cell" index in the sample name (e.g. early2cell_0r-1 and
early2cell_0r-2 → embryo early2cell_0r; zy1 … zy4 are four zygotes).
Included stages and names: zygote (zy1–4), early 2-cell (early2cell_*), mid 2-cell (mid2cell_*), late 2-cell
(late2cell_*), 4-cell (4cell_*), 8-cell (8cell_1/2/5/8), 16-cell (16cell_1/4/5/6), early/mid/late blastocyst
(earlyblast_*, midblast_*, lateblast_*).
Excluded before analysis: *split*, *pooled*, *smartseq2* (technical / different protocol), C57twocell_* (different cross,
unspecified timing), fibroblasts, liver RNA/cells.
Genes: symbol → mouse Ensembl via the tAge gene table (unique symbols only); summed if several symbols map to one ID.
QC: exclude an embryo if log10(detected genes, reads >= 1) < median − 3 × 1.4826 × MAD.
Clock: EN_Chronoage_Multispecies_Multitissue_scaleddiff (primary), reference = zygotes (no oocytes in R1).

Primary statistics (bootstrap embryos within stage, 2000 replicates, seed 20260917):
- S1 = mean(late 2-cell) − mean(early 2-cell), 95% CI.
- Adjacent intervals among zygote, early 2C, mid 2C, late 2C, 4C, 8C, 16C: Δ and P_boot(most negative).
Decision R1:
- REPLICATED: S1 < 0 with CI excluding 0 AND the most negative pre-blastocyst adjacent interval is early→mid 2C or
  mid→late 2C with combined P_boot >= 0.6.
- NOT REPLICATED: S1 CI includes 0 or S1 > 0.
- PARTIAL: S1 < 0 with CI excluding 0 but the most negative adjacent interval lies elsewhere.
Blastocyst intervals are descriptive.

## Dataset R2 (secondary): GSE66582 (Wu et al. 2016, Nature), bulk embryos, B6N × DBA2N
Runs: MII oocyte, zygote, early 2-cell, 2-cell, 4-cell, 8-cell (2 replicates each; ICM descriptive; knock-down and
mESC runs excluded). Raw FASTQ re-quantified with the project kallisto index (GRCm39 r112 cDNA+ncRNA) and summed
to Ensembl genes. Reference = MII oocytes.
With n = 2 per stage no interval CI is interpretable; the prespecified readout is direction only:
- CONSISTENT: mean(2-cell) − mean(early 2-cell) < 0 and it is the most negative adjacent interval among
  MII→zygote, zygote→early 2C, early 2C→2C, 2C→4C, 4C→8C.
- INCONSISTENT otherwise.

## Secondary (both datasets)
Repeat with the Gate 2a V2 variant (all maternal/zygotic dynamic genes removed before preprocessing; gene sets defined
within each dataset with the Gate 2a thresholds, using the earliest stage available as "oocyte" and the next as
"zygote": R1 zygote + early 2C; R2 MII + zygote).

## Outputs
results/gate2b_R1_embryo_tage.tsv, results/gate2b_R1_intervals.tsv, results/gate2b_R2_*.tsv, results/gate2b_decision.md
