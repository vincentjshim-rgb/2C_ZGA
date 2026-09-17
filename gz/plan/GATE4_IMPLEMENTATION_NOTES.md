# Gate 4 — implementation notes (written 2026-09-17, before any PSI value was computed)

These notes fill in details that `GATE4_PLAN_frozen.md` leaves open. They do not change any question, dataset, filter,
threshold or decision rule of the frozen plan. Where the plan is silent, SUPPA2 defaults are used and named here.

## Corrections to the running record (not to the plan)

- The frozen plan already excludes P4 GSE162345 from Gate 4 ("PE36, excluded from splicing"). The statement in
  `results/SUMMARY_gz_*.md` and the review page that P4's Gate 4 status still had to be settled by addendum was wrong;
  no addendum is needed and none is written.
- Read-length condition (≥ 95 nt) checked on the downloaded data: P1 GSE280522 mixed 72/150 nt reads, mean ≈ 110 nt per
  read; P2 GSE221985 150 nt; P3 GSE300734 150 nt. All three meet the condition as stated.

## Inputs

- Events: `data/annotation/suppa/mm112_{SE,A5,A3,MX,RI,AF,AL}_strict.ioe`, generated 2026-09-17 with
  `generateEvents -f ioe -e SE SS MX RI FL -b S` from `Mus_musculus.GRCm39.112.gtf.gz`
  (sha256 a1ad4101…223c). Event counts equal the Windows run of 2026-09-15 (SE 20,851; A5 9,399; A3 10,913;
  MX 2,239; RI 4,731; AF 32,964; AL 7,429). The seven files are concatenated into one ioe (one header).
- Libraries: the Gate 3 QC-passing libraries of P1–P3, with the same QC rule (so SRR30827777 and SRR34172919 are
  excluded as in Gate 3). Group labels come from the same title parser as `11_gate3_analysis.py`.
- Transcript TPM: the `tpm` column of each kallisto `abundance.tsv`; transcript IDs with the version suffix removed.
  Written in SUPPA's expression format (header = sample names only, one field fewer than the data rows).

## SUPPA calls

- `psiPerEvent -i <merged ioe> -e <dataset TPM> -o <out>` with the default total filter (`-f 0`). The plan's own
  expression filter is applied afterwards.
- `diffSplice -m empirical -gc -i <merged ioe> -p <control E2C psi> <control L2C psi> -e <control E2C tpm>
  <control L2C tpm>`, all other parameters SUPPA defaults: area 1000, lower bound 0, alpha 0.05, TPM threshold 1.0,
  NaN threshold 0, unpaired, mean ΔPSI. The p-value used is the `p-val` column of the resulting `.dpsi` file
  (gene-corrected, as `-gc` requests). Only the control arm enters `diffSplice`, as the plan states.

## Computation order

1. PSI for every event in every QC-passing library of the dataset.
2. Plan filter: host-gene TPM (sum of the event gene's transcripts) ≥ 5 in ≥ 80% of control E2C libraries **and**
   ≥ 80% of control L2C libraries; PSI non-missing in ≥ 80% of libraries in every arm × stage group of the dataset.
3. Control ZSA events: filtered events with |mean PSI(control L2C) − mean PSI(control E2C)| ≥ 0.10 **and** diffSplice
   `p-val` < 0.05. The ΔPSI used for the threshold and for `s_e` is computed from the filtered PSI table with means over
   available libraries, not taken from `diffSplice`.
4. Progression score per library: mean over control ZSA events (NaN-skipping) of
   `s_e · (PSI_e − mean control E2C PSI_e) / |ΔPSI_e|`.
5. `P` per arm, interaction `P(perturbed) − P(control)`, rescue `P(A485+Dux) − P(A485)`; bootstrap of libraries within
   arm × stage, 2000 replicates, `numpy.random.default_rng(20260919)`.
6. Decision exactly as the plan, on point estimates; the bootstrap intervals are reported beside them. In P1 the
   interaction entering the rule is the A485 arm (A485+Dux is the rescue), as in Gate 3. The plan's three bullets are
   applied literally in order; if ZSA PRESENT fails, the verdict line still follows the bullets and the absence of ZSA
   is stated next to it, so that a "NOT ZGA-DEPENDENT" without ZSA cannot be read as a finding about ZSA.
7. Clock coupling (descriptive): arm × stage group means of progression against Gate 3 V0 tAge group means from
   `results/gate3_library_tage.tsv`, Spearman per dataset and pooled; and Spearman between the four progression
   interactions (P1 A485, P1 A485+Dux, P2, P3) and the matching Gate 3 V0 interactions. No inference is drawn.

## Known limits carried from the plan

PSI from pseudoalignment transcript estimates, not junction reads; annotated events only; P2 has n = 2 per group and
P3 control E2C has n = 2 after QC, so `diffSplice`'s empirical distribution is built from very few replicates in those
datasets.
