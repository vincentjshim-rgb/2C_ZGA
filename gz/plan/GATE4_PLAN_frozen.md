# Gate 4 plan — zygotic splicing activation (ZSA): is it ZGA-dependent, and does it track the clock decrease?
(frozen 2026-09-15, before any PSI value was computed)

## Question
Within the 2-cell stage, does alternative splicing change (zygotic splicing activation, reported in mouse/human),
is that change attenuated when ZGA is blocked, and does splicing progression co-vary with the tAge decrease measured on
the same libraries (Gate 3)?

## Data
Gate 3 primary datasets with reads >= 95 nt: P1 GSE280522 (DMSO / A485 / A485+Dux, E2C & L2C, n = 4), P2 GSE221985
(Ctrl / Tardbp matKO, n = 2), P3 GSE300734 (WT / Brg1 matKO, n = 3). P4 GSE162345 (PE36) excluded from splicing.
Same libraries and QC as Gate 3.

## Method (SUPPA2 v2.3)
- Events: `generateEvents` on Ensembl GRCm39 r112 GTF, local event types SE, RI, A5, A3, MX, AF, AL, boundary "S".
- Transcript TPM: kallisto abundance.tsv (same quantification as Gate 3), transcript IDs matched to the GTF (version stripped).
- PSI per library: `psiPerEvent`.
- Event filter per dataset: host-gene TPM (sum of transcripts) >= 5 in >= 80% of libraries of control E2C and of control L2C;
  PSI non-missing in >= 80% of libraries in every arm × stage group used.
- Control ZSA events: events with |mean PSI(control L2C) − mean PSI(control E2C)| >= 0.10 and SUPPA `diffSplice`
  (empirical method, gene-correction on) p < 0.05, computed on the control arm only.

## Readouts
1. ZSA magnitude (control): number and fraction of filtered events that are control ZSA events, by event type.
2. Splicing progression score per library: mean over control ZSA events of s_e × (PSI_e − mean control E2C PSI_e) / |ΔPSI_e|,
   where s_e = sign of the control ΔPSI; control E2C ≈ 0, control L2C ≈ 1 by construction.
   Arm drop in progression: P = mean score(L2C) − mean score(E2C); interaction = P(perturbed) − P(control)
   (prediction: < 0, i.e. less progression when ZGA is blocked); rescue = P(A485+Dux) − P(A485) (prediction > 0).
3. Coupling to the clock (descriptive): across arm × stage group means in P1–P3, Spearman correlation between the
   progression score and Gate 3 V0 tAge; and across arms, correlation of interaction terms (progression vs clock).

Bootstrap: libraries within arm × stage, 2000 replicates, seed 20260919.

## Decision rules
- ZSA PRESENT: control ZSA events >= 100 in >= 2 of 3 datasets.
- ZSA ZGA-DEPENDENT: ZSA PRESENT, interaction < 0 in >= 2 of 3 datasets, and P1 rescue > 0.
- ZSA NOT ZGA-DEPENDENT: interaction >= 0 in >= 2 of 3 datasets.
- Otherwise MIXED. The clock coupling (readout 3) is reported only; with <= 6 group means per dataset no inferential claim.

## Caveats (stated in advance)
PSI from pseudoalignment-based transcript estimates, not junction reads; ZGA block changes expression of zygotic genes,
so the expression filter is applied in both control stages; annotated events only (no novel junctions).

## Outputs
gz/results/gate4_events_filtered_<gse>.tsv, gz/results/gate4_zsa_summary.tsv, gz/results/gate4_progression.tsv,
gz/results/gate4_decision.md
