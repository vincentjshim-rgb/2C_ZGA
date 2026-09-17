# POST HOC — cross-fitted check of the Gate 4 progression score (P1 GSE280522)

Plan `plan/POSTHOC_gate4_crossfit_frozen.md`. Reproduction check (no hold-out) matched Gate 4 to 1e-6 with an identical event set.

READING: **ATTENUATED BUT INFLATED** — mean I_A485 -0.2733, negative in 14/16 folds; in-sample inflation of the control arm 0.2097.

| quantity | mean over folds | fold spread (2.5–97.5%) |
|---|---|---|
| n_events | 571.1250 | 429.0000 – 808.0000 |
| overlap_with_gate4 | 0.7632 | 0.6012 – 0.8673 |
| P_ctrl_in | 1.0000 | 1.0000 – 1.0000 |
| P_ctrl_out | 0.7903 | 0.3800 – 1.0716 |
| P_A485 | 0.5170 | 0.4679 – 0.5746 |
| P_Dux | 0.7340 | 0.6260 – 0.8217 |
| I_A485 | -0.2733 | -0.5461 – 0.0998 |
| I_Dux | -0.0563 | -0.3075 – 0.2460 |
| R | 0.2170 | 0.1462 – 0.2512 |

Fold spread is not a confidence interval: folds share libraries. P2 and P3 are not cross-checked.
