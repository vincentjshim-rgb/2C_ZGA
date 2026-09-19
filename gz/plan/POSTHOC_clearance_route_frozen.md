# Post-hoc note — through which genes does simulated maternal clearance move the clock?
(frozen 2026-09-19, before the number was computed; display-level follow-up of Gate 2a, no new data, no test)

## Why
Section 2.2 reports that scaling maternal transcripts in oocytes by the observed clearance factor reproduces 0.70 of the
mouse two-cell decrease, while Section 2.4 reports that 99.5% of the decrease in GSE280522 comes from clock genes that
meet neither the maternal nor the zygotic definition. Read side by side the two statements look contradictory. They are
not, if the clearance acts on the clock through the preprocessing (relative-log-expression normalisation and per-library
scaling shift the features of every other gene when abundant maternal transcripts are removed) rather than through the
coefficients of maternal clock genes. This note quantifies that route with the exact decomposition already used.

## Computation
Mouse GSE225056, exactly as in `src/05_gate2a_composition.py` (same QC, gene sets, clearance ratio r and oocyte
reference). For the simulated oocytes against the real oocytes, the clock change Δ_sim = Σ_g β_g Δx_g is split into
(i) the sum over clock genes that belong to the maternal set and (ii) the sum over all other clock genes. Δ_sim must
equal the stored `gate2a_simulation.tsv` value for mouse (tolerance 10⁻⁹) before anything is written. Also reported:
the number of clock genes in the maternal set.

## Reading (decided in advance)
Whatever the split, it is reported as one sentence in Section 2.4 (or the Discussion) and one line in the audit log; no
figure, no table, no gate. If route (ii) carries most of Δ_sim, the sentence says that maternal clearance reaches the
clock mainly through the normalisation of the remaining genes; if route (i) does, the sentence says the two results
overlap in the maternal clock genes.

Script `src/23_posthoc_clearance_route.py`; output `results/posthoc_clearance_route.tsv`.
