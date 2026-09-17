# POST HOC — the A485+Dux arm did restore zygotic transcription; the clock score still did not move
(plan `plan/POSTHOC_rescue_frozen.md`; script `src/14_posthoc_rescue.py`; run 2026-09-17. Descriptive. Gate 3 stands
at MIXED and its rescue term R = −0.0071 [−0.0629, 0.0572] is unchanged.)

The plan named two separable explanations: (a) the rescue did not happen biologically, (b) it happened but not on what
the clock reads. The measurements answer (b), and then show why the score stayed flat, which neither option anticipated.

## 1. Dux co-expression restored the zygotic programme

Zygotic set (13 genes, control-arm definition), mean log2(CPM+1) at late 2-cell:

| arm | score | difference [95% CI] |
|---|---|---|
| control E2C | 0.49 | — |
| control L2C | 3.57 | — |
| A485 L2C | 1.10 | vs control −2.469 [−2.723, −2.214] |
| A485+Dux L2C | 2.76 | vs control −0.802 [−1.335, −0.358]; vs A485 **+1.667 [1.139, 2.037]** |

The maternal set moves in the opposite, smaller way (A485 +0.688 [0.227, 1.228] above control; A485+Dux +0.386
[−0.011, 0.792]). Markers confirm it: Zscan4f 6.71 → 9.15, Zscan4c 5.28 → 8.16, Zscan4d 5.25 → 7.70 from A485 to
A485+Dux, in each case above the control L2C value (7.45, 6.39, 6.19); Tdpoz and Tcstv genes return to control levels.
So the rescue arm is not a failed manipulation — explanation (a) is excluded.

## 2. The clock's own down-contributors were also largely restored

Summed contributions of the 666 genes that push the score down in control: control −1.4030, A485 −0.7592,
**A485+Dux −1.1867** (85% of control). Over the 20 largest down-contributors: −0.2955, −0.1304, −0.2210. Per-gene
contribution profiles correlate with control at 0.647 (A485) and 0.833 (A485+Dux). Klf9 recovers from −0.0173 to
−0.0283 against −0.0327 in control; Pi4k2a, Gpatch4, Psmb5, Chn1 recover similarly.

## 3. Why the score did not move: the opposing flow grew too

Per-arm flows (unplanned follow-up computed after §2, descriptive):

| arm | negative flow | positive flow | net D |
|---|---|---|---|
| control | −1.4030 | +1.1655 | −0.2375 |
| A485 | −1.1529 | +1.0964 | −0.0566 |
| A485+Dux | −1.3803 | +1.3166 | −0.0637 |

A485 damps both flows. A485+Dux restores the downward flow to 98% of control but the upward flow reaches 113% of it, so
the residual — which is what the clock reports — stays where A485 left it. The genes gaining the most upward
contribution in this arm relative to control are Cd74 (+0.0241), Parp3 (+0.0233), Cst7 (+0.0137), Icam1 (+0.0119),
Upp1 (+0.0116), Gpc4 (+0.0115), Cdkn1a (+0.0111), Lgals3 (+0.0089).

That list reads as a stress- and inflammation-associated set, and Dux overexpression is reported elsewhere to impose
cellular stress, but this analysis tests nothing of the kind: the genes were selected by their clock contribution, no
gene set was tested, n = 4 libraries per arm, and the labels come from an adult multi-tissue clock. It is recorded as an
observation worth a prespecified test, not as a finding.

## 4. What this changes

The Gate 3 rescue term is not evidence that the decrease is independent of zygotic transcription. It is evidence that a
single scalar clock difference can stay flat while both of its component flows are restored or exceeded. Any statement
about the rescue arm should be made on the decomposition, not on R alone. Gate 3 remains MIXED; no rule was changed and
no number in `results/gate3_arm_drops.tsv` was recomputed.
