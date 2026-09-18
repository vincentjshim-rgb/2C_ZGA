# Post-hoc plan — is the within-2-cell clock decrease explained by global transcriptome remodelling?
(frozen 2026-09-18, before any number of this analysis was computed; extends `POSTHOC_gate3_contribution_frozen.md`)

## Why
Early embryogenesis is accompanied by a genome-wide fall and then rise of mRNA abundance. A clock is a weighted sum, so
a decrease could in principle arise from every clock gene moving by roughly the same amount (a global shift that the
preprocessing did not remove) rather than from specific genes. The source study of the clock ran the analogous control
for its embryonic trajectory by recomputing partial tAge from only up- or only down-regulated genes
(Tyshkovskiy et al., 2026, Extended Data Fig. 9). Reviewers will ask the same of the decomposition reported here.

## Data and inputs
The V0 preprocessed feature matrices (`scaled_diff`) and quality-control set of Gate 3, loaded exactly as in
`src/13_posthoc_gate3_contribution.py`: P1 GSE280522 (control, A485, A485+DUX) and P3 GSE300734 (control, Brg1 matKO).
Coefficients: the 1,839 non-zero coefficients of the clock. Genes absent from a dataset keep Δx = 0 and therefore
contribute zero, as in the existing decomposition.

## Readouts (no new model, no new data)
For an arm, with Δx_g the change in the preprocessed feature between late and early two-cell libraries:

1. **Exact split of the drop into a global and a specific part**
   D = Σ_g β_g Δx_g = G + S, where G = (Σ_g β_g) · mean_g(Δx_g) is what the drop would be if every clock gene moved by
   the average amount, and S = Σ_g β_g (Δx_g − mean_g Δx_g) is the gene-specific remainder.
2. **Permutation null**: 2,000 random permutations of the assignment of the coefficient vector to genes, holding the
   observed Δx vector fixed (seed 20260920). Report the null mean (which must equal G to numerical tolerance), the
   2.5th and 97.5th percentiles, and the observed D's position in the null.
3. **Direction-restricted partial drops**: D recomputed on genes with Δx > 0 only and with Δx < 0 only, and on genes
   whose control expression rises or falls across the window (the categories of Figure 3a, recomputed here for
   completeness).
4. The same three readouts for the interaction I = D(A485) − D(control) in P1, with Δx replaced by the difference of the
   two arms' Δx.

## Prespecified reading rule
For the control arm, the decrease is

- **NOT EXPLAINED by global remodelling** if, in both P1 and P3, |S| > |G| and the observed D lies below the 2.5th
  percentile of its permutation null;
- **PARTLY EXPLAINED** if exactly one dataset meets both conditions;
- **EXPLAINED** if |G| ≥ |S| in both datasets.

For the P1 interaction the same rule is applied with the observed I above the 97.5th percentile of its null.

## Decided before seeing the result (to keep the paper from growing)
- If the verdict is NOT EXPLAINED: add at most two sentences to Results 2.4 and one clause to Methods 4.8; no new
  figure, no new table.
- If the verdict is PARTLY EXPLAINED or EXPLAINED: do not add a result to the paper; add one limitation sentence to the
  Discussion instead, and report the numbers in the audit log only.
- Either way this stays post hoc and descriptive: no gate verdict changes, and the permutation percentile is reported as
  a position within a re-assignment null, not as a test of a biological hypothesis.

## Outputs
`results/posthoc_global_remodelling.tsv` (long format: dataset, arm, quantity, value) and
`results/posthoc_global_remodelling.md` (short report). Script `src/21_posthoc_global_remodelling.py`. Seed 20260920.
