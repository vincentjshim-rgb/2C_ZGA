# POST HOC — cross-fitted check of the Gate 4 splicing progression score (P1 only)
(frozen 2026-09-17, before running `src/16_posthoc_gate4_crossfit.py`; written after the Gate 4 verdict was seen)

**Status.** POST HOC. It does not change the Gate 4 verdict or rule. It measures how much of the Gate 4 interaction
comes from scoring the control arm on the same libraries used to select and scale the events.

## The problem it addresses

In Gate 4, control ZSA events were selected (|ΔPSI| ≥ 0.10, diffSplice p < 0.05) and the progression scale was set
(control E2C → 0, control L2C → 1) on the control libraries, and the control arm was then scored on those same
libraries. Control P = 1 is therefore in-sample; every perturbed arm is out-of-sample. Selection on noisy ΔPSI makes
an in-sample value larger than an out-of-sample value of the same biology, so the interaction
`P(perturbed) − P(control)` is biased downward by an unknown amount.

The P1 rescue term R = P(A485+Dux) − P(A485) compares two out-of-sample arms and is not affected by this bias; it is
recomputed here only for completeness.

## Scope

P1 GSE280522 only: its control arm has 4 E2C and 4 L2C libraries. P2 (2 + 2) and P3 (control E2C 2 after QC) cannot be
split while keeping ≥ 2 libraries per selection group. Their Gate 4 interactions remain "not cross-checked" and are
reported that way.

## Design (deterministic; no random seed needed)

Folds: every choice of one held-out control E2C library and one held-out control L2C library, 4 × 4 = 16 folds.
In each fold, with the remaining 3 + 3 control libraries as the **selection set**:

1. PSI values are the per-library values already computed in Gate 4 (`data/gate4/GSE280522.psi`); PSI of one library
   does not depend on other libraries, so `psiPerEvent` is not re-run.
2. Filter: host-gene TPM ≥ 5 in ≥ 80% of selection-set E2C **and** selection-set L2C libraries; PSI non-missing in
   ≥ 80% of selection-set E2C, selection-set L2C, and each perturbed arm × stage group. Held-out libraries do not enter
   the filter.
3. `diffSplice -m empirical -gc` on the selection set only, SUPPA defaults otherwise (as in Gate 4).
4. Fold ZSA events: filtered events with |ΔPSI_sel| ≥ 0.10 and p < 0.05, ΔPSI from selection-set means.
5. Scale from the selection set only: base = mean selection-set E2C PSI, s_e and |ΔPSI_e| from the selection set.
6. Score with that fold's events and scale: the two held-out control libraries, and all A485 and A485+Dux libraries.
7. Fold quantities:
   - `P_ctrl_out` = score(held-out L2C) − score(held-out E2C)
   - `P_ctrl_in` = mean score(selection L2C) − mean score(selection E2C) (= 1 by construction; reported as a check)
   - `P_A485`, `P_Dux` = arm L2C mean − arm E2C mean
   - `I_A485 = P_A485 − P_ctrl_out`, `I_Dux = P_Dux − P_ctrl_out`, `R = P_Dux − P_A485`

A reproduction check runs first: with no library held out (4 + 4 selection, control scored in-sample) the code must
return the Gate 4 values (P control 1.000, P A485 0.5568, P A485+Dux 0.7869) to 1e-6; otherwise the run stops.

## Summaries

Mean over the 16 folds of each quantity; the 2.5th–97.5th percentile range across folds is reported as **fold spread**,
explicitly not a confidence interval (folds share libraries). Inflation = 1 − mean `P_ctrl_out`. The number of fold
events and their overlap with the Gate 4 event set are reported.

## Reading, fixed now

On `I_cf` = mean `I_A485` over folds:

- **SURVIVES**: `I_cf ≤ −0.20` and `I_A485 < 0` in at least 15 of 16 folds. The Gate 4 splicing attenuation under A485
  may then be reported, with the cross-fitted value, as the effect size.
- **ATTENUATED BUT INFLATED**: `I_cf < 0` but either `I_cf > −0.20` or fewer than 15 of 16 folds negative. The direction
  may be reported; the Gate 4 magnitude may not.
- **CONSTRUCTION**: `I_cf ≥ 0`. The Gate 4 interaction is attributed to the selection construction and is not reported
  as a splicing effect.

The −0.20 threshold is about half of the Gate 4 value (−0.443) and is fixed here before the run. Whatever the outcome,
P2 and P3 stay "not cross-checked", and no change is made to the Gate 4 verdict line in `results/gate4_decision.md`;
the reading is added beside it.

## Outputs

`results/posthoc_gate4_crossfit_folds.tsv` (one row per fold), `results/posthoc_gate4_crossfit.md`, log
`logs/16_posthoc_gate4_crossfit.log`, working files under `data/gate4/crossfit/`.
