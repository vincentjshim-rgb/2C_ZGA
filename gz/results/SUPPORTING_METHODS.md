## Supporting Methods

These paragraphs give the detail removed from the main-text Methods for length. Section numbers refer to the Methods of
the main text.

**S4.4 Clock preprocessing.** Genes with at least 10 counts in at least 20% of the libraries of a dataset were kept,
mapped to Entrez identifiers, normalised by relative log expression, log10-transformed and scaled per library, then
aligned to the clock gene list and centred on the per-gene median of the reference libraries named in Section 4.4.

**S4.5 Gate 1 quality control and permutation.** A GSE225056 library was excluded when log10 of its detected genes fell
below the species median minus 3 × 1.4826 × MAD; clock feature coverage per species was 0.56–0.69. The within-species
permutation of the clock changes used 10,000 permutations (seed 20260915).

**S4.6 Gate 2b inclusion rules.** In GSE45719, split, pooled, Smart-seq2, C57 two-cell, fibroblast and liver samples
were excluded and the 8-cell and 16-cell stages were restricted to four embryos each; the cells of each retained embryo
were summed to 39 pseudobulks, none of which failed the detected-gene rule. In GSE45719 and GSE66582 the two earliest
stages available took the roles of oocyte and zygote in the maternal and zygotic gene definitions.

**S4.7 Gate 3 arm selection and seeds.** V2 used the control early two-cell libraries as the early stage and the control
late two-cell libraries as the later stage. Bootstraps used 2,000 replicates with seed 20260918. In GSE248499 the SCNT
arm is the deposit's "SCNT, cont" libraries (late two-cell), which are the arms matching the Kdm4d and Kdm3a injections;
in GSE235547 the late two-cell total RNA-seq libraries were used.

**S4.8 Reassignment null.** The null that assigns the coefficient vector to genes at random used 2,000 draws with seed
20260920.

**S4.9 Marker list.** The first run of the rescue analysis omitted Zfp352 from the prespecified marker list through a
scripting error; this was corrected before submission and recorded in the audit log.

**S4.10 SUPPA2 parameters.** Events were generated with generateEvents using strict boundaries (one import was updated
for the installed statsmodels version), and PSI was computed with psiPerEvent under the default filter after removing
transcript version suffixes. An event was kept when its host gene had TPM ≥ 5 in ≥ 80% of the control early libraries
and, separately, in ≥ 80% of the control late libraries. Control splicing-activation events were called with diffSplice
using the empirical method with gene correction. Arm progression, interactions and the rescue difference used seed
20260919.

**S4.11 Seeds.** Bootstraps and permutations used fixed seeds of numpy.random.default_rng: Gate 1 20260914, the
zygotic-timing analysis 20260915, Gate 2a 20260916, Gate 2b R1 20260917, Gate 3 and the rescue analysis 20260918,
Gate 4 20260919 and the reassignment null 20260920. The GSE66582 readout, the exact decomposition and the cross-fitting
involve no random draw.

**S4.12 Figure details.** Bootstrap intervals drawn in Figure 1 were recomputed with the seeds of the original analyses
and equal the reported intervals. Expression panels show
log2(CPM + 1) over all Ensembl genes.


**S4.2 Audit log.** Plans, scripts, outputs, exclusions and failed runs are recorded with checksums in a dated audit
log in the public repository.

**S4.4 Reference libraries.** Preprocessed features were centred on the per-gene median of the oocytes of the same
species (GSE225056), the zygotes (GSE45719), the MII oocytes (GSE66582), the control early two-cell libraries (primary
perturbation series) or the control libraries (secondary series).

**S4.6 Dynamic-gene thresholds.** Maternal genes were those with counts per million (CPM) ≥ 10 in ≥ 80% of oocytes and
of zygotes and a median ≤ 1/8 of the oocyte median in at least one later stage; zygotic genes those with CPM < 1 in
≥ 80% of oocytes and of zygotes and CPM ≥ 10 in ≥ 50% of the libraries of at least one later stage.

**S4.7 Full Gate 3 rule.** The primary datasets were classified as ZGA-dependent when I > 0 in all informative datasets
(at least two) and R < 0, as not ZGA-dependent when I ≤ 0 in at least two informative datasets, and as mixed otherwise,
separately in V0 and V2; an overall conclusion required the two variants to agree.

**S4.8 Reassignment null, in full.** The null assigns the coefficient vector to genes at random while holding the
observed feature changes fixed (2,000 draws, seed 20260920). The frozen reading rule asked for |S| > |G| and the
observed drop below the 2.5th percentile of the null, in both datasets; |S| exceeded |G| by factors of 34 and 15, but
the drop lay at the 2.6th (GSE280522) and 8.1st (GSE300734) percentile, so the rule did not classify the outcome.

**S4.10 Event filters.** An event was kept when the host gene had TPM ≥ 5 in ≥ 80% of the control early libraries and,
separately, in ≥ 80% of the control late libraries, and PSI was available in ≥ 80% of the libraries of every
arm-by-stage group.

**S4.12 Figure specification.** Every figure is 167 mm wide, rendered at 600 dpi, with no text below 6 pt at final size.

**S4.4 Reference libraries, in full.** Features were centred on the per-gene median of the oocytes of the same species
(GSE225056), the zygotes (GSE45719), the MII oocytes (GSE66582), the control early two-cell libraries (primary
perturbation series) or the control libraries (secondary series).

**S4.5 Gate 1 library counts.** Libraries passing quality control: mouse 75 of 76, cow 82 of 82, pig 78 of 81, rabbit
81 of 82.

**S4.7 Secondary-arm selection.** In GSE248499 the SCNT arm is the deposit's "SCNT, cont" libraries, which match the
Kdm4d and Kdm3a injections; in GSE235547 the late two-cell total RNA-seq libraries were used. Each GSE162345
collection window was compared with the control of that window.

**S4.8 Reassignment control, in full.** |S| exceeded |G| by factors of 34 and 15, but the drop lay at the 2.6th
(GSE280522) and 8.1st (GSE300734) percentile of the null, so the frozen rule, which asked for both, did not classify
the outcome.

**S4.10 Read-length condition and SUPPA2 settings.** The plan restricted the splicing analysis to primary datasets
whose mean read length was at least 95 nt: GSE280522 has 150-nt reads with adapter-trimmed shorter reads, mean about
110 nt, and GSE221985 and GSE300734 have 150 nt; GSE162345 (51 + 25 nt) was excluded. Events were generated with
strict boundaries.

**S4.12 Figure specification.** Every figure is 167 mm wide, rendered at 600 dpi, with no text below 6 pt at final
size. Bootstrap intervals drawn in Figure 1 were recomputed with the seeds of the original analyses and equal the
reported intervals.

**S4.1 Screen of the perturbation series.** Fourteen GEO series of mouse preimplantation RNA-seq with a ZGA
perturbation were retrieved and screened against the frozen requirement for a primary dataset: early *and* late
two-cell libraries in the same study, in a perturbed arm *and* a matched control arm.

| Outcome | Series | Reason |
|---|---|---|
| Primary | GSE280522, GSE221985, GSE300734 | both stages in both arms |
| Secondary, single stage | GSE248499, GSE235547 | perturbation and control at the late two-cell stage only |
| Secondary, moved before scoring | GSE162345 | the deposit is not an early-to-late two-cell pair (addendum 2) |
| Descriptive only | GSE195760 | both stages present, but every two-cell arm is a nuclear transfer, so no control arm exists for an interaction |
| Not eligible | GSE166338, GSE196671, GSE214878, GSE229740, GSE262039, GSE269417, GSE298245 | no early/late two-cell split, so the drop cannot be computed |

**S4.4 Coverage of the weighted clock genes.** Feature coverage as reported for Gate 1 (0.56–0.69) is the fraction of
all clock input features present in a dataset. Over the 1,839 genes that carry a non-zero coefficient — the genes that
determine the value — coverage is 0.50 for mouse GSE225056 (912 genes), 0.74 for GSE280522 and 0.82 for GSE300734.
Genes that a dataset does not detect are median-imputed by the clock pipeline and contribute nothing to any
difference, so they neither add to nor subtract from a reported change; they do mean that the cross-species
comparison rests on half of the weighted clock.
