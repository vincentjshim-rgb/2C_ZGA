# Gate 1 plan — does the transcriptomic age drop track ZGA or cell division? (frozen 2026-09-14, before any clock score on embryo data)

Study question (umbrella): is the embryonic molecular-age reset tied to reprogramming transitions
(ZGA, totipotency exit) rather than to developmental progression / cell division per se?
Gate 1 is the cross-species timing test. Exploratory; no causal language.

## Data
GSE225056 (Smart-seq+5', single oocytes/embryos, one lab, same protocol). Deposited per-species gene
count matrices `GSE225056_<Species>_count_5prime.txt.gz` (Ensembl gene rows; TE rows dropped).
Primary species: Mouse (76), Cow (82), Pig (81), Rabbit (82; "Zika" = rabbit hybrid line, not infection).
Rhesus (19; 1–6 per stage) descriptive only. ESC libraries excluded.
Unit: embryo within stage (no donor/female IDs deposited — limitation reported).
Gate 1a uses these deposited 5' counts. Any positive result must be confirmed in Gate 1b by full-length
re-quantification (kallisto) of mouse and cow raw reads before it is quoted.

## Clock
tAge Python port `src/tage_py.py`, validated against TACO golden values (`results/00_tage_py_validation.tsv`, PASS).
Primary model: EN_Chronoage_Multispecies_Multitissue_scaleddiff (normalised output, no lifespan rescaling).
Preprocessing run separately per species on all embryos of that species; reference samples = oocytes
(so tAge is relative to the species' oocyte). Non-mouse genes mapped to mouse via Ensembl one-to-one orthologs.
Model SHA256: tools_tAge_models/SHA256SUMS.

## QC (before scoring)
Exclude a library if log10(number of genes with count >= 1) is below species median − 3×1.4826×MAD.

## Stage order and intervals
Mouse: Oocyte, Zygote, Early-2-cell, Late-2-cell, 4-cell, 8-cell, 16-cell
Cow, Rabbit: Oocyte, Zygote, 2-cell, 4-cell, 8-cell, 16-cell, Morula
Pig: Oocyte, Zygote, 2-cell, 4-cell, Day2 (5–8 cell), Day3 (8–16 cell), Morula
Cell-number label of each interval: fertilisation (Oocyte→Zygote), 1→2, (mouse only: within-2-cell, no division),
2→4, 4→8, 8→16, 16→morula.

## Major ZGA interval (fixed from literature, not from these data)
Mouse: Early-2-cell → Late-2-cell (no cell division in this interval)
Pig: 2-cell → 4-cell
Cow: 8-cell → 16-cell
Rabbit: 8-cell → 16-cell

## Statistics
tAge per embryo; stage mean; Δ_k = mean(stage k+1) − mean(stage k); negative Δ = younger.
Bootstrap: resample embryos within each stage, 2000 replicates, seed 20260914 → 95% CI for every Δ and
P_boot(interval with most negative Δ = ZGA interval).

## Decision rules
- INFORMATIVE: in >= 3 of 4 primary species at least one Δ has a 95% CI excluding 0. Otherwise UNINFORMATIVE (stop; clock does not move).
- ZGA-LINKED (pass): the most negative Δ falls in the species' ZGA interval with its CI < 0 and P_boot >= 0.6,
  in >= 3 of 4 species, including >= 1 early-ZGA species (mouse or pig) and >= 1 late-ZGA species (cow or rabbit).
- DIVISION-LINKED: the most negative Δ falls in the same cell-number interval in >= 3 species.
- Otherwise MIXED: report, do not adopt the reprogramming framing without redesign.

## Sensitivity (reported, not used for the decision)
yugenediff version; EN_Mortality_Multispecies; mouse-specific EN_Chronoage_Mouse in mouse; tAge residualised on
log10 total gene counts and log10 detected genes within species; interval with the largest absolute change;
clock feature coverage per species; rhesus descriptive panel.

## Outputs
results/gate1_embryo_tage.tsv, results/gate1_intervals.tsv, results/gate1_decision.md, figure later.
