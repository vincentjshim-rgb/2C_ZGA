# Gate 3 plan — does the within-2-cell transcriptomic-age decrease require zygotic genome activation / reprogramming?
(frozen 2026-09-15, before any clock score on the datasets below)

## Background (established so far)
Within-2-cell decrease of the tAge clock replicated in three wild-type mouse series (GSE225056, GSE45719, GSE66582);
a variable part (0–60%) is compositional (maternal clearance). Gate 3 asks whether the decrease depends on ZGA.

## Prediction
If the decrease is tied to reprogramming/ZGA: embryos in which ZGA is blocked show an attenuated early→late 2-cell
decrease relative to controls of the same study (interaction > 0), and restoring ZGA (rescue arm) restores the decrease.
Interpretation limit: attenuation shows dependence on zygotic transcription in this window, not a causal "rejuvenation"
mechanism; no causal verbs.

## Datasets (raw FASTQ re-quantified with kallisto 0.51.1, GRCm39 r112 cDNA+ncRNA, summed to Ensembl genes)
Primary (early and late 2-cell within one study):
- P1 GSE280522: E2C, L2C × DMSO (control) / A485 (p300/CBP inhibitor, minor-ZGA block) / A485+Dux (rescue); n = 4.
- P2 GSE221985: E2C, L2C × Ctrl / maternal Tardbp (TDP-43) KO; n = 2. (MII, 1C used only for gene-set definition checks.)
- P3 GSE300734: E2C, L2C × WT / maternal Smarca4 (Brg1) KO; n = 3.
- P4 GSE162345: SCNT E2C, L2C × no α-amanitin (control) / α-amanitin (transcription block); n = 3; PE36 reads.
Secondary (late 2-cell only; single-stage contrasts):
- S1 GSE248499: IVF vs SCNT control vs SCNT+Kdm4d (and other SCNT treatments, descriptive).
- S2 GSE235547: ICSI vs SCNT-EGFP vs SCNT-Obox3; siCtrl vs siObox3 (late 2-cell).
Descriptive: GSE195760 SCNT time course (no IVF arm).

## Readouts
Clock EN_Chronoage_Multispecies_Multitissue_scaleddiff (Python tAge port). Preprocessing per dataset on all included
libraries; reference = control E2C libraries (P1–P4) or control late 2-cell libraries (S1, S2).
- V0: all genes.
- V2 (composition-robust): remove maternal and zygotic dynamic genes defined in the CONTROL arm only:
  maternal = CPM >= 10 in >= 80% control E2C and control L2C median CPM <= 1/8 of control E2C median;
  zygotic = CPM < 1 in >= 80% control E2C and CPM >= 10 in >= 50% control L2C.
Per arm: drop D = mean tAge(L2C) − mean tAge(E2C). Interaction I = D(perturbed) − D(control). Rescue restoration
R = D(A485+Dux) − D(A485) (P1 only).
Uncertainty: bootstrap of libraries within arm × stage (2000 replicates, seed 20260918); with n = 2 (P2) the CI is
reported but not used; direction only.

## Decision rules (evaluated on V0 and V2 separately; the conclusion requires both)
- QC precondition per dataset: the control arm shows D < 0 (otherwise the dataset is uninformative for Gate 3).
- ZGA-DEPENDENT: among informative primary datasets, I > 0 in >= 3 of 4 (or all, if fewer than 4 are informative, with
  >= 2 informative), AND in P1 R < 0 (rescue restores part of the decrease).
- NOT ZGA-DEPENDENT: I <= 0 in >= 2 informative primary datasets (drop unchanged or larger when ZGA is blocked).
- MIXED otherwise.
Secondary single-stage prediction: tAge(L2C perturbed/SCNT) − tAge(L2C control) > 0 and rescue arms lower it; reported
only.

## Also prespecified
Library QC: exclude libraries with kallisto pseudoalignment < 30% or detected genes (count >= 1) below median −
3×1.4826×MAD within dataset. FASTQ retained for Gate 4 (splicing; separate plan, frozen before junction counting).

## Outputs
gz/results/gate3_<dataset>_library_tage.tsv, gz/results/gate3_arm_drops.tsv, gz/results/gate3_decision.md
