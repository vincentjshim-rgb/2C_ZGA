# POST HOC — which genes produce the within-2-cell clock decrease, and which contributions A485 removes
(frozen 2026-09-17, before running `src/13_posthoc_gate3_contribution.py`; written after the Gate 3 verdict was seen)

**Status.** POST HOC and descriptive. It is labelled post hoc wherever it is reported. It does not change Gate 3, does
not re-open any gate, adds no dataset, and produces no p-value and no new decision rule. Its purpose is to say what the
clock decrease consists of, because a single score moving down is not yet a biological statement.

## Why this is exact, not an approximation

The clock is a scikit-learn pipeline: median imputation → mean-centring (`StandardScaler(with_std=False)`) →
`SelectKBest` (10,487 features) → `ElasticNet` (1,839 non-zero coefficients). With `species=None`, as used in Gate 3,
no species factor is applied, so for a sample

```text
tAge = intercept + Σ_g coef_g · (x_g − centre_g)
```

and for a difference between two groups of libraries the intercept and the centres cancel:

```text
D = mean tAge(L2C) − mean tAge(E2C) = Σ_g coef_g · ( mean x_g(L2C) − mean x_g(E2C) )
```

where `x_g` is the `scaled_diff` feature of gene g (Entrez space) produced by the same preprocessing as Gate 3.
Genes absent from the matrix are imputed with one constant and cancel in the difference, contributing exactly 0.
The per-gene terms therefore sum to the drop reported in `results/gate3_arm_drops.tsv`; the script asserts this
(tolerance 1e-9).

## Input (identical to the Gate 3 run, no re-quantification)

P1 GSE280522, V0 gene space, the QC-passing libraries used in Gate 3 (SRR30827777 excluded by the frozen rule),
preprocessing reference = control E2C libraries. P3 GSE300734 is decomposed the same way as a direction check.
P2 GSE221985 is not decomposed (n = 2, direction not stable across V0/V2).

## Quantities

Per gene g, within the control arm and within the A485 arm:

```text
c_control_g = coef_g · ( mean x_g(control L2C) − mean x_g(control E2C) )     Σ = D(control)
c_A485_g    = coef_g · ( mean x_g(A485 L2C)    − mean x_g(A485 E2C)    )     Σ = D(A485)
Δ_g         = c_A485_g − c_control_g                                          Σ = I
```

Reported alongside each gene: Entrez ID, symbol, clock coefficient, log2 fold change of CPM (control L2C vs control
E2C), and the maternal / zygotic / neither label.

## Prespecified outputs

1. `results/posthoc_gate3_contribution_genes.tsv` — every clock feature with a non-zero contribution in P1, sorted by
   |c_control|: the columns above plus the P3 equivalents.
2. `results/posthoc_gate3_contribution_summary.tsv` —
   - concentration: share of D(control) from the top 1, 5, 10, 25, 50, 100 genes; number of genes reaching 50% and 80%;
   - direction split: D(control) partitioned into four cells by (expression rises or falls from E2C to L2C) ×
     (coefficient positive or negative), i.e. how much of the decrease comes from genes being switched on versus genes
     being cleared;
   - programme split: share of D(control) from genes labelled maternal, zygotic, or neither;
   - the same three splits for Δ (what A485 removes).
3. `results/posthoc_gate3_contribution.md` — the numbers above in prose, marked POST HOC.

## Gene labels

The maternal and zygotic labels use the Gate 3 V2 definitions, applied to the control arm of the same dataset
(maternal: CPM ≥ 10 in ≥ 80% of control E2C and control L2C median ≤ 1/8 of control E2C median; zygotic: CPM < 1 in
≥ 80% of control E2C and CPM ≥ 10 in ≥ 50% of control L2C). Deviation recorded here: they are applied in Entrez space,
the space the clock features live in, rather than in Ensembl space as in `11_gate3_analysis.py`; gene counts in the two
spaces therefore differ slightly.

## Limits fixed before the numbers are seen

- Coefficients come from an adult multi-species multi-tissue clock. A gene carrying a large contribution means the
  clock weights it and it moves in this window; it is not evidence that the gene has a role in the embryo.
- No causal verb. "Contribution to the score", not "driver of rejuvenation".
- Ranked lists are descriptive; no enrichment test is run here, and no gene is promoted to a claim on this basis.
- If the decrease turns out to be spread thinly over many genes, that is the result and is reported as such; no
  threshold is tuned afterwards to produce a shortlist.
