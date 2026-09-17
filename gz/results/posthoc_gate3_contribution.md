# POST HOC — what the within-2-cell clock decrease is made of
(plan `plan/POSTHOC_gate3_contribution_frozen.md`; script `src/13_posthoc_gate3_contribution.py`; run 2026-09-17.
Descriptive. No p-value, no gate, no change to Gate 3, which stands at MIXED.)

The decomposition is exact: per-gene terms sum to the drops in `results/gate3_arm_drops.tsv` (asserted, 1e-9).
P1 GSE280522: D(control) = −0.2375, D(A485) = −0.0566, I = +0.1809. P3 GSE300734: D(control) = −0.1714,
D(Brg1 matKO) = −0.0623, I = +0.1091. 1,355 clock genes carry a non-zero contribution in P1, 1,509 in P3.

## 1. The decrease is a small imbalance between two large opposing flows

In P1 the negative contributions sum to −1.4030 and the positive ones to +1.1655; the observed drop, −0.2375, is the
17% residual. Reaching half of the total absolute movement takes 161 genes. Statements of the form "the top 10 genes
explain 77% of the drop" are arithmetically true against the net value but misleading: those same ten genes are a small
part of the movement on either side.

## 2. Activation and clearance contribute in similar measure

Splitting P1's D by the direction of expression change (control E2C → L2C) and the sign of the clock coefficient:

| cell | contribution | reading |
|---|---|---|
| expression rises × negative coefficient | −0.746 | genes switched on in the window push the score down |
| expression falls × positive coefficient | −0.650 | genes cleared in the window push the score down |
| expression rises × positive coefficient | +0.531 | opposing |
| expression falls × negative coefficient | +0.626 | opposing |

So the decrease is not attributable to activation alone or to maternal clearance alone; the two are of comparable size
here (0.75 versus 0.65). P3 gives the same picture (+4.28 / +4.10 of D against −3.42 / −3.97, D = −0.1714).

## 3. It is mostly not the canonical maternal or zygotic genes

Using the Gate 3 V2 definitions applied in the clock's Entrez space, 99.5% of P1's D comes from genes labelled
*neither* maternal nor zygotic (1,346 genes); 8 maternal genes contribute 0.2% and the single zygotic clock gene 0.3%.
P3: 92.2% neither, 8.0% maternal, −0.2% zygotic. This is consistent with V2 barely changing the interaction
(I = +0.1809 → +0.1810 in P1), and it means the clock is not simply reading the maternal-to-zygotic switch.

## 4. The same genes carry the decrease in two independent datasets

Per-gene contributions in P1 and P3 correlate (Pearson 0.739, Spearman 0.606 over the 1,839 clock genes); 28 of the 50
largest down-contributors are shared. The two datasets are different laboratories, different perturbations (A485 versus
maternal Brg1 loss) and different mice, so the composition of the decrease is reproducible at the gene level even though
the Gate 3 verdict is MIXED.

## 5. What A485 removes is concentrated on those same genes

Contribution and change are anticorrelated (r = −0.565 between `c_control` and `delta`): the genes that push the score
down in control are the genes that stop doing so under A485. The ten largest control down-contributors account for 63%
of the whole interaction I.

Largest down-contributors in P1 (full table: `results/posthoc_gate3_contribution_genes.tsv`):

| symbol | coefficient | c control | c A485 | Δ | log2FC control L2C/E2C |
|---|---|---|---|---|---|
| Klf9 | −0.0199 | −0.0327 | −0.0173 | +0.0154 | +4.48 |
| Neto2 | −0.0227 | −0.0237 | −0.0028 | +0.0209 | +2.92 |
| Pi4k2a | −0.0159 | −0.0182 | −0.0032 | +0.0150 | +3.03 |
| Smyd2 | +0.0318 | −0.0180 | −0.0139 | +0.0041 | −0.87 |
| Gpnmb | +0.0208 | −0.0175 | −0.0056 | +0.0119 | −1.57 |
| Gpatch4 | −0.0135 | −0.0161 | −0.0036 | +0.0125 | +2.89 |
| S100a6 | +0.0244 | −0.0158 | +0.0048 | +0.0206 | −0.72 |

Genes with a negative coefficient that are strongly induced in the window (Klf9, Neto2, Pi4k2a, Gpatch4, Psmb5;
log2FC +2.9 to +4.5) lose almost all of their contribution when minor ZGA is blocked. Genes with a positive coefficient
that fall in the window (Smyd2, Gpnmb, S100a6) lose part of theirs.

## 6. Limits

The coefficients belong to an adult multi-species multi-tissue clock. A large contribution means the clock weights that
gene and the gene moves in this window; it is not evidence of a role in the embryo, and none of these genes is promoted
to a claim here. No enrichment test was run. Section 4 is an observation about reproducibility of composition, not a
statistical test of it.
