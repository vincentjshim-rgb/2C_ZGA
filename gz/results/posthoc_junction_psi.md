# Post hoc: junction-level PSI in GSE280522 (plan: plan/POSTHOC_junction_psi_frozen.md)

Skipped-exon events in the Gate 4 filtered set: 4076 (155 control-activated).
Libraries aligned and scored: 23 of 23. Library-event pairs at >= 10 unique junction reads: 59,482.
Junction coordinate offset chosen from the data: +1.

## 1. Agreement of the two quantifications

Pearson r = 0.779, Spearman rho = 0.783 over 59,482 library-event pairs scored by both.
Prespecified reading: **supports the Gate 4 quantification**.

## 2. Does the Gate 4 result reproduce on junction PSI?

Progression (control = 1 by construction): control 1.00, A485 0.39, A485+DUX 0.90 over the 31 control-activated SE events that are scored in every library and whose junction control change is at least 0.10 (of 78 scored in every library).
Junction depth differs by arm (1,499 to 3,491 events scored per library), so the score is computed on complete cases only; scoring each library on whatever it covers would average different event sets and is reported here as the reason for that restriction.
I = -0.61, R = +0.51. Prespecified reading: **reproduces the Gate 4 direction**.
Compare with the cross-fitted transcript-level values (I = -0.27, R = +0.22); the in-sample transcript-level magnitude is not the comparator.

## 3. Events selected by rule for the read-level illustration

- `ENSMUSG00000028693;SE:4:116469278-116471527:116471607-116476067:-` dPSI_control = -0.57
- `ENSMUSG00000071533;SE:16:55842492-55844587:55844749-55849966:-` dPSI_control = +0.51
- `ENSMUSG00000026511;SE:1:181952540-181956249:181956317-181958901:+` dPSI_control = -0.50
