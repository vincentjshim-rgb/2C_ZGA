# Post hoc: does the decomposition hold in human embryos? (plan: plan/POSTHOC_human_decomposition_frozen.md)

GSE36552, 20 embryo pseudobulks from 90 cells; embryonic stem cell samples excluded by the plan.

## 3. How much of the clock does this dataset carry?

Weighted clock genes detected: **1,370 of 1,839** (0.74). Mouse comparators: 912 (GSE225056), 1,355 (GSE280522), 1,509 (GSE300734).

## 1. Is there a decrease across the prespecified interval?

8-cell -> Morula: **D = -0.098** [-0.232, +0.036], 3 and 2 embryo pseudobulks.

Clock value by stage (mean over embryo pseudobulks, relative to oocytes):

- Oocyte: +0.026 (n = 3)
- Zygote: +0.031 (n = 3)
- 2-cell: +0.010 (n = 3)
- 4-cell: +0.146 (n = 3)
- 8-cell: +0.105 (n = 3)
- Morula: +0.007 (n = 2)
- Late blastocyst: -0.039 (n = 3)

## 2. Is it a residual of opposing contributions?

Downward contributions sum to **-1.243**, upward to **+1.145**, and they add to the reported change of -0.098.
Residual fraction = |D| / max(|down|, |up|) = **0.08** (mouse: 0.17 in GSE280522, 0.12 in GSE225056; plan threshold 0.50).

**Prespecified reading: SAME STRUCTURE AS MOUSE — the decrease is a small residual of opposing contributions**

## Follow-up (computed after the planned items, not in the frozen plan)

Is it the same genes? Over the 1,098 clock genes with a non-zero contribution in both the human 8-cell to morula interval and the mouse GSE280522 control window, the per-gene contributions are essentially uncorrelated: Pearson r = 0.072, Spearman rho = 0.112, and 10 of the 50 largest downward contributors are shared.

For comparison, two mouse studies measuring the same interval in the same species agreed at r = 0.74. The near-cancellation therefore reproduces across species and intervals while the gene identity behind it does not, which is consistent with reading the largest contributors as a property of this clock in a given window rather than as a ranking of genes.
