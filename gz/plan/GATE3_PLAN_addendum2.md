# Gate 3 plan — addendum 2 (2026-09-16, in effect)

**Status.** Written and approved on 2026-09-16, before any Gate 3 score: no result from any Gate 3 dataset had been
computed or seen, and `src/11_gate3_analysis.py` had not been run. Approved decisions: P4 GSE162345 moves from primary
to secondary; its two collection windows are reported separately (`T45`, `T54`); the frozen decision rules are applied
unchanged, which with three primary datasets means ZGA-DEPENDENT requires `I > 0` in all three when all are informative.
Recorded in `analysis/logs/AUDIT_LOG.md` (2026-09-16).

**Subject.** P4 GSE162345 — the repository metadata do not match the reading of this dataset in the frozen plan.

---

## 1. Why the frozen reading fails

The frozen plan (2026-09-15) lists P4 as "SCNT E2C, L2C × ± α-amanitin; PE36", i.e. an early→late 2-cell sampling under
a transcription block, from which the drop `D` and the interaction `I` are computed like P1–P3.

GEO (series record and GSM5610457–68) describes a different experiment: C2C12 cells were transplanted (fused) into
2-cell embryos either at 21 hpi ("early") or at 30 hpi ("late"), cultured for 24 h in 0.1 µg/ml demecolcine with or
without α-amanitin, and collected at 45 hpi (`early2cell_NT`) or 54 hpi (`late2cell_NT`); pools of 5 embryos; strains
C57BL/6, DBA/2, C3H; enucleation is not stated.

Consequences:

- "early" and "late" name the recipient stage at transfer, not two samplings of one embryo population across the
  2-cell window. `D = mean tAge(L2C) − mean tAge(E2C)` therefore does not measure the early→late 2-cell change that
  Gate 3 is about, and `I` does not measure its attenuation.
- The two groups differ in three ways at once (recipient stage, collection time 45 vs 54 hpi, and time of α-amanitin
  onset relative to major ZGA). They are not a stage pair.
- Both arms within a window receive demecolcine, so that drug cancels within a window-matched contrast.

Verified on the downloaded data (2026-09-16, all 12 runs quantified): 3 libraries per group; 13.2–21.3 M read pairs;
pseudoalignment 77.0–81.1% (the highest of the six Gate 3 datasets). R1 is 51 nt and R2 is 24 nt; R2 is shorter than
k = 31, so kallisto's paired mode uses R1 only and reports `estimated average fragment length: 0.0` for all 12 libraries
(identical treatment across the dataset; effective lengths are uncorrected).

## 2. Options considered

**A. Drop P4 from Gate 3 entirely.** Primary becomes P1–P3; nothing from GSE162345 is reported. Loses the only
pharmacological transcription block in the set (α-amanitin), which is the most direct ZGA block available here.

**B. Keep computing `D` and `I` for P4 and label them descriptive.** Rejected: the numbers would be printed in the same
columns as P1–P3 and read as the same quantity, which is exactly the confusion this addendum exists to prevent.

**C. Re-cast P4 as a window-matched, single-stage transcription-block contrast (recommended).** Treat GSE162345 like the
secondary datasets S1/S2: compare α-amanitin against control *within each collection window*, never across windows.

## 3. What is adopted (option C)

P4 GSE162345 moves from primary to **secondary**. Its role becomes a transcription-block contrast at a fixed collection
time in a nuclear-transfer reprogramming context.

Arms (3 libraries each): `control_45hpi`, `alpha-amanitin_45hpi`, `control_54hpi`, `alpha-amanitin_54hpi`.

Readout, per window, on V0 (V2 is not defined for single-stage datasets, as in the frozen plan):

```text
T45 = mean tAge(alpha-amanitin_45hpi) − mean tAge(control_45hpi)
T54 = mean tAge(alpha-amanitin_54hpi) − mean tAge(control_54hpi)
```

Reference for preprocessing: the control libraries of this dataset (both windows), as for S1/S2. Uncertainty: the same
bootstrap as the other secondaries (2000 replicates, seed 20260918); with n = 3 per arm the CI is reported, direction is
what is read.

Prespecified direction: if the within-2-cell tAge decrease requires zygotic transcription, blocking transcription leaves
embryos at a higher tAge at the same collection time, i.e. `T > 0`. This is an association in a nuclear-transfer
context; no causal verb, and no claim about normal fertilised embryos.

**Explicitly not computed for P4:** any 45 hpi → 54 hpi difference, in either arm. The windows differ in recipient stage
as well as time, so such a difference is not `D` and is not reported.

## 4. Effect on the Gate 3 decision rules

Primary datasets become **P1 GSE280522, P2 GSE221985, P3 GSE300734** (three, not four). The frozen rule already covers
this case and is applied unchanged:

- ZGA-DEPENDENT: `I > 0` in **all** informative primary datasets, with at least 2 informative, **and** `R < 0` in P1.
- NOT ZGA-DEPENDENT: `I ≤ 0` in ≥ 2 informative primary datasets.
- MIXED otherwise.
- Unchanged: the per-dataset QC precondition (control arm shows `D < 0`), the requirement that V0 and V2 agree before
  any conclusion, and the library QC (pseudoalignment ≥ 30%; detected genes ≥ median − 3×1.4826×MAD).

P4 contributes to none of these counts. It is reported beside S1 and S2.

## 5. Code changes required before the analysis runs

In `src/11_gate3_analysis.py` (to be made, then re-checksummed in the audit log):

1. `annotate()`, GSE162345 branch (line 28–30): return the window and the window-tagged arm
   (`early2cell_NT*` → `('NT45', 'control_45hpi' | 'alpha-amanitin_45hpi')`, `late2cell_NT*` → `NT54` equivalents)
   instead of mapping `early`/`late` to `E2C`/`L2C`.
2. Dataset list (line 83): `('GSE162345', 'P')` → `('GSE162345', 'S')`.
3. Primary loop (line 144): remove `'GSE162345'` from the primary list.
4. Secondary reference mapping (line 123): map `alpha-amanitin_45hpi` → `control_45hpi` and
   `alpha-amanitin_54hpi` → `control_54hpi`; keep both control arms out of the reported contrasts.
5. Output: P4 rows appear in `results/gate3_secondary_L2C.tsv`; rename that file to
   `results/gate3_secondary_single_stage.tsv`, since its rows are no longer all late-2-cell.

Applied 2026-09-16. While making change 4 a pre-existing fault was found in the same expression and fixed before any
score: every arm starting with `SCNT` was pointed at reference arm `control_ICSI`, which exists only in GSE235547, so
all three S1 GSE248499 contrasts would have used an empty reference and returned NaN. The reference is now written per
dataset — S1: IVF (`control`) for `SCNT_control`, `SCNT_Kdm4d` and `SCNT_Kdm3a`, so a rescue arm shows as a smaller
difference, as in S2; S2 unchanged (`control_ICSI`, and `control_siCtrl` for `siObox3`); P4 as above — and an empty
reference now raises instead of producing NaN. Contrasts produced: S1 3, S2 3, P4 2.

## 6. Quantification of P4 (kept as is)

The 12 libraries stay as quantified (paired call, R1 effectively used, fragment length 0.0 for all 12). The bias is
uniform within the dataset and the contrasts are within-dataset, so it does not distort `T45`/`T54`. If any reported
statement turns on P4, a sensitivity re-quantification in single-end mode (`--single -l 200 -s 30`, R1 only) is run and
reported next to it; that re-quantification is not run by default.

## 7. Open item, deferred to the Gate 4 plan

Gate 4 (splicing) currently expects the same FASTQ. With 51 nt single-usable reads and no fragment-length estimate,
transcript-level TPM for P4 is weak, so P4 may need to be excluded from Gate 4 as well. That belongs to the Gate 4 plan
and is not decided here; it must be settled before any PSI is computed.
