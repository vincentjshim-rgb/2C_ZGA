# Gate 1-alt plan — does the transcriptomic-age drop align with the ZGA timing measured in the same embryos?
(frozen 2026-09-14, after Gate 1 = MIXED and before any zygotic-gene score was computed)

## Status
POST HOC / EXPLORATORY. Motivated by Gate 1: the literature-fixed ZGA intervals may not match when ZGA occurs in
GSE225056. The clock values are already known (Gate 1); only the ZGA timing is new. A positive result here cannot
be quoted as confirmatory; it requires replication in independent data (e.g. mouse GSE45719/GSE66582, another cow
series) under a new frozen plan.

## Data and preprocessing
Same libraries, QC and identifier mapping as Gate 1 (plan + addendum 1): GSE225056 5' gene counts, mouse / cow /
pig / rabbit, genes as mouse one-to-one orthologs (mouse genes as is). CPM = count / library total over mapped
genes × 1e6. Rhesus excluded.
Clock-feature exclusion: genes that are features of EN_Chronoage_Multispecies_Multitissue_scaleddiff (mouse Entrez →
mouse Ensembl via the tAge gene table) are removed before defining the gene sets below, so the ZGA and clock readouts
share no genes.

## Data-driven ZGA readout (primary, A)
Zygotic genes (per species): CPM < 1 in >= 80% of oocytes AND >= 80% of zygotes, and CPM >= 10 in >= 50% of the
libraries of at least one later stage.
Embryo score A = fraction of CPM from zygotic genes.
Data-driven major ZGA interval = stage-to-stage interval with the largest positive change in mean A.

## Secondary readout (B, maternal clearance)
Maternal genes: CPM >= 10 in >= 80% of oocytes and >= 80% of zygotes, and median CPM in at least one later stage
<= 1/8 of the oocyte median. Score B = fraction of CPM from maternal genes; clearance interval = largest negative change.

## Clock readout (fixed, from Gate 1)
Primary-analysis Δ per interval from `results/gate1_intervals.tsv` (not recomputed).

## Statistics
Bootstrap embryos within stage, 2000 replicates, seed 20260915: 95% CI for every Δ of A and B, and
P_boot(interval = argmax ΔA). Pooled alignment: Spearman correlation between ΔA and −Δclock across all intervals of
the four species, with 10,000 permutations of interval order within species (seed 20260915).

## Decision rules
- Timing informative: data-driven ZGA intervals take >= 2 distinct cell-number labels across the four species
  (otherwise ZGA and division count cannot be separated in these data → NOT SEPARABLE).
- ALIGNED: the most negative clock interval equals the data-driven ZGA interval in >= 3 of 4 species, AND pooled
  Spearman > 0 with permutation p < 0.05.
- PARTIAL: exactly one of the two ALIGNED conditions holds.
- NOT ALIGNED: neither holds.
Secondary (reported only): same comparison with the maternal-clearance interval (B); agreement of data-driven vs
literature ZGA intervals.

## Outputs
results/gate1alt_stage_scores.tsv, results/gate1alt_intervals.tsv, results/gate1alt_gene_sets.tsv,
results/gate1alt_decision.md, figures/FigGZ_2_zga_vs_clock.
