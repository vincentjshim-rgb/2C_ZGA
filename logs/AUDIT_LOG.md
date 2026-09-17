
## 2026-09-13 — canonical input change: gene matrices rebuilt with ncRNA

**What was wrong.** `metadata/tx2gene.tsv` was built from `mm_cdna.fa.gz` only, but the
kallisto index `data/processed/mm_full.idx` was built from cDNA + ncRNA. At gene
aggregation (`src/14`, `src/27`) all 29,930 ncRNA transcripts were dropped. Consequences:
gene-level TPM columns did not sum to 1e6 (B6N median ~590,000; per-library scale
differed), and ncRNA genes (Xist, Tsix, H19, Meg3, Airn, Kcnq1ot1, rRNA, miRNA, snoRNA)
were absent from every analysis.

**Second defect.** Some libraries were quantified from FASTQ files truncated when a
download or cutadapt run aborted (`results/verify_truncated_libraries.tsv`). Three pass QC:
BlastoB6NY2208, BlastoFVBO3122, BlastoFVBxxY2xx2xx25.

**Found by.** `find-the-real-finding` workflow (run wf_347f7323-083), then verified
directly: ncRNA transcripts present in tx2gene = 0 of 29,930; unmapped TPM in one B6N
library = 43.5%.

**Change.** New files, old files kept for comparison, not overwritten:
- `metadata/tx2gene_full.tsv` (`src/50`) — 145,918 transcripts, 53,954 genes
- `data/processed/{b6n,fvb}_{tpm,counts}_full.parquet` (`src/51`) — gene_id level
- `data/processed/gene_annotation_full.tsv`

**Impact to be assessed.** Every established number (reproduction r, variance
partition, mixed models, leave-one-sire-out, calibration, MDE, sex calls) is re-run on the
full matrices and compared side by side with the old values before any figure or
manuscript number is reused. Truncated libraries are re-downloaded or excluded.

**Also withdrawn on 2026-09-13 (reporting errors, not data changes).**
- "52.6% of read pairs carry TSO" and "pseudoalignment 21.6% → 67.5%": the first counted
  chance 3-nt adapter matches, the second compared different denominators. Paired
  trimmed/untrimmed comparison on 30 B6N embryos: median |dLFC| 0.009.
- "EV11 orientation is inverted": most likely a young-minus-mature reporting convention.
- "388 embryos": 371 deposited libraries.

### 2026-09-13 — impact of the ncRNA fix on established B6N numbers (`src/55`, `results/rerun_old_vs_full.tsv`)

Same 154 embryos (155 QC-pass minus BlastoB6NY2208, truncated); prespecified QC gate
re-applied on the full matrix selects the identical 154. The two zero-TPM libraries
(BlastoB6NO2218, BlastoB6NO5108) were never in the QC-pass set: they were lost to a
duplicate-worker race (log shows OK then QUANT_FAIL two seconds later), not to embryo
quality. They and the 4 truncated libraries are being re-quantified.

| quantity | old matrix | full matrix |
|---|---|---|
| expressed genes | 12,976 | 16,547 |
| published DEGs present | 220 | 222 |
| reproduction Pearson / Spearman / sign | 0.952 / 0.870 / 99.5% | 0.953 / 0.872 / 99.5% |
| mixed model beta, p | +0.425, 0.0015 | +0.424, 0.0015 |
| sire-level exact permutation p | 0.0095 | 0.0095 |
| median age share / df chance | 0.64% / 0.65% | **0.58% / 0.65%** |
| median sire-within-age share / df chance | 7.35% / 5.23% | **6.89% / 5.23%** |
| leave-one-sire-out, SD (95% CI), perm p | −0.26 (−1.05 to +0.53), 0.45 | −0.28 (−1.02 to +0.46), 0.38 |
| minimum detectable \|log2FC\| | 0.364 | 0.380 |

Conclusion: no established conclusion changes. Two framing consequences: the median
paternal-age share is at or below its degrees-of-freedom chance level, and the sire share
exceeds chance by only about 1.7 percentage points, so "sire exceeds age in 96.7% of
genes" is mostly 8 df versus 1 df and is withdrawn as a headline.

### 2026-09-13 19:23 — P0 complete: all 371 libraries intact, final established numbers

**Re-quantification.** Six libraries re-run with verified complete input (cutadapt input
pairs = ENA-declared pairs for every one): BlastoB6NY2108 (6,332,186), BlastoB6NY2208
(13,012,934), BlastoFVBO3122 (14,980,489), BlastoFVBxxY2xx2xx25 (8,468,741), plus the two
duplicate-worker losses BlastoB6NO2218 (9,807,872) and BlastoB6NO5108 (13,653,692).
Truncation audit (`src/52`): 0 of 371 truncated; 0 zero-TPM libraries.

**Process failures during re-quantification (recorded, all resolved).** (1) Job list written
by Windows Python with CRLF, so read-2 URLs ended in "\r"; (2) killing a job's parent left
retry subshells alive, which later deleted files the relaunched job was using and wrote a
stale completion marker; (3) Git-Bash /tmp and Windows-Python /tmp are different
directories, so one relaunch ran on an empty list and wrote a false completion marker.
Fixes: CR-stripped lists in `work/`, PID-level kills, unique completion markers
(`LAST3B_DONE`, `FINAL_REQUANT_DONE`), chain waits for zero live workers.

**QC on the full matrices (`src/58`).** B6N 157/169 pass (was 155); FVB 186/202 pass.

**Final established numbers, B6N, 157 embryos / 10 sires (`results/rerun_old_vs_full.tsv`).**
| quantity | original (155, old matrix) | final (157, full matrix) |
|---|---|---|
| reproduction Pearson / Spearman / sign | 0.951 / 0.869 / 99.5% | 0.953 / 0.871 / 99.5% (222 genes) |
| embryo-level beta, p | +0.448, 2.3e-13 | +0.450, 1.1e-13 |
| mixed model beta (SE), p | +0.424 (0.0895), 0.0015 | +0.425 (0.0919), 0.0017 |
| SE inflation | 1.61 | 1.66 |
| sire-level exact permutation p | 0.0095 | 0.0095 |
| leave-one-sire-out, SD (95% CI), perm p | −0.27 (−1.05 to +0.50), 0.41 | −0.25 (−1.02 to +0.51), 0.44 |
| minimum detectable \|log2FC\| | 0.364 | 0.378 |

**Degrees-of-freedom-adjusted variance partition (`src/56`, `results/df_adjusted_partition.tsv`).**
B6N: age median 0.58% vs analytic chance 0.64% and age-label permutation null 0.47%
(p = 0.29); sire-within-age 6.76% vs analytic chance 5.13% and within-(age x batch)
permutation null 6.09% (p < 1/300); excess over analytic chance 1.63 points (95% CI
1.55–1.72). FVB (age term aliased with flowcell): sire-within-age 4.45% vs chance 2.70%,
permutation null 2.36%.

**Planted-effect calibration on the final set (`src/57`, 8 spike-in gene sets × 9 levels).**
Observed −0.261 SD, upper 95% bound 0.500 SD. Recovered mean: 0.25 log2FC → +0.386 (0/8
replicates reach the bound); 0.30 → +0.684 (8/8). Transferable effects of 0.30 log2FC or
larger are excluded; 0.25 is not. This replaces the interpolated "about 0.25–0.3".

## 2026-09-13 — Figure 4 programme bounds on the final set (`src/64`)

Plan frozen before output: `results/FIG4_PLAN_frozen.md`. B6N final set, 157 embryos / 10 sires (4 mature / 6 young),
sire as unit; exact permutation over 210 labelings, batch-stratified over 36; Welch 95% CI; MDE at 80% power.
Results in `results/fig4_programme_bounds.tsv`; replaces the pre-repair table in MANUSCRIPT_FACTS_v2 §6.
No programme shows a detectable sire-level association (all CIs include 0): residual +0.57 SD (−0.20, +1.34) p 0.18;
tau +0.0018 (−0.0026, +0.0062) p 0.40; TE−ICM +0.23 SD p 0.43; male fraction +0.19 (−0.09, +0.47) p 0.15;
2C-transient calibrated +0.03 log2 p 0.89; PEG−MEG +0.25 log2 (−0.09, +0.60) p 0.16; dispersion ratio 1.00 p 0.96.
Changes vs pre-repair probes: lineage estimator switched from uncalibrated NNLS ICM share to TE−ICM marker z;
2C score now background-calibrated against the maternal-clearance control set; imprinted dosage newly computed.
Sex calls on final set match Figure 1 (94 M / 57 F / 6 ambiguous). Exploratory; no causal language.

## 2026-09-14 — Research direction change: paternal-age reanalysis → embryonic age reset × reprogramming

Decision by the user after the progress review (`results/PROGRESS_MAP_KR.md`):
- The Dura reanalysis / Biology of Reproduction path is paused. Its results (Figures 1–4, MANUSCRIPT_FACTS_v2)
  are kept unchanged as provenance; no claim from it is withdrawn or promoted.
- The original CLAUDE.md estimand (phase-conditioned ZGA-resolution residual × paternal age) was never taken
  through Step 3 locked validation; the Aging Cell go/no-go conditions 2–3 failed and 4 was untested. It is
  recorded as stopped, not as confirmed or refuted.
- New umbrella question: is the embryonic molecular-age reset tied to reprogramming transitions (ZGA,
  totipotency exit, gastrulation re-methylation) rather than to development/cell division per se?
  Layers: mouse perturbations (SCNT ± rescue, ZGA mutants, totipotent cells), cross-species timing test,
  human replication. Proposal: `results/IDEA_ground_zero_reprogramming.md`.
- Work directory: `gz/`. Public clocks are applied, not built (tAge, Tyshkovskiy 2026; Python port of the
  R preprocessing validated against the package's TACO golden values, max abs diff 7e-8 / 3e-5).
- Gate 1 plan frozen before any embryo clock score: `gz/plan/GATE1_PLAN_frozen.md` (+ .sha256).

## 2026-09-14 — Gate 1 (cross-species timing of transcriptomic-age drop vs major ZGA): MIXED

Plan `gz/plan/GATE1_PLAN_frozen.md` + addendum 1 (identifier mapping only: cow RefSeq → Ensembl; rhesus via
M. mulatta orthologs; species scored only with all stages). Script `gz/src/01_gate1_cross_species.py`.
Data GSE225056 deposited 5' gene counts; clock EN_Chronoage_Multispecies scaleddiff (Python tAge port, TACO-validated).
Libraries after QC: mouse 75, cow 82, pig 78, rabbit 81 (rhesus 18, descriptive). Clock feature coverage 0.56–0.69.

Most negative stage-to-stage change (95% bootstrap CI; P that it is the most negative interval):
- Mouse: Early-2-cell → Late-2-cell (major ZGA, no division) −0.19 (−0.25, −0.13), P 1.00 → ZGA criterion met.
- Cow: 16-cell → morula −0.39 (−0.47, −0.31), P 0.85; ZGA interval 8 → 16-cell also negative −0.31 (−0.38, −0.25), P 0.15.
- Pig: Day3 → morula −0.14 (−0.22, −0.05), P 0.84; ZGA interval 2 → 4-cell −0.04 (−0.11, +0.03).
- Rabbit: 4 → 8-cell −0.13 (−0.18, −0.08), P 0.99; ZGA interval 8 → 16-cell +0.02 (−0.04, +0.08).
Fertilisation (oocyte → zygote) is positive in mouse, pig, rabbit (+0.11 to +0.13, CIs exclude 0).
All species informative. ZGA criterion met in 1/4 → VERDICT MIXED; all sensitivity analyses (yugene, mortality clock,
depth/detection-residualised, cow excluded) also MIXED. Per frozen rule, the reprogramming framing is not adopted
from Gate 1; no Gate 1b re-quantification is triggered. Any reinterpretation (e.g. alternative ZGA-stage
definitions, cumulative-drop onset) is post hoc and must be labelled as such.

## 2026-09-14 — Gate 1-alt (POST HOC): clock drop vs ZGA timing measured in the same embryos: PARTIAL

Plan `gz/plan/GATE1ALT_PLAN_frozen.md` (frozen after Gate 1, before any zygotic-gene score; sha256 logged).
Script `gz/src/03_gate1alt_zga_timing.py`. Clock-feature genes excluded from ZGA/maternal gene sets.
Zygotic gene sets: mouse 683, pig 87, cow 81, rabbit 49 (genes considered after clock-gene exclusion: mouse 44,347;
non-mouse ~5,000–5,600 because only one-to-one orthologs remain — asymmetry, limits the non-mouse ZGA readout).
Data-driven major ZGA (largest zygotic-share gain; P_boot): mouse Early→Late 2-cell (1.00, = literature);
pig Day2→Day3 / 8→16 cell (0.66; literature 2→4); cow 8→16 (1.00, = literature); rabbit 16-cell→morula (1.00; literature 8→16).
Clock most-negative interval equals data ZGA in 1/4 (mouse); equals maternal-clearance interval in mouse and rabbit;
in pig and cow the clock drop is one interval after the data ZGA. Pooled Spearman(dZygotic, −dClock) 0.524,
within-species permutation p 0.0034. VERDICT PARTIAL. Exploratory only; hypothesis for independent replication:
the transcriptomic-age drop coincides with or follows (0–1 stage) zygotic activation / maternal clearance.

## 2026-09-14 — Gate 2a (composition artefact test): INTERMEDIATE

Plan `gz/plan/GATE2A_PLAN_frozen.md`; script `gz/src/05_gate2a_composition.py`. Clock change in the data-driven ZGA
interval (95% bootstrap CI):
- Mouse E2C→L2C: V0 −0.191 (−0.251, −0.134); V1 (non-clock maternal/zygotic genes removed) −0.179 (retention 0.94);
  V2 (all 3,556 dynamic genes removed, 1,798 of them clock features) −0.081 (−0.126, −0.032), retention 0.42;
  composition-only simulation on oocytes (maternal genes × 0.356) −0.134 (−0.176, −0.091) = 0.70 of V0.
- Cow 8→16C: V0 −0.314; V1 −0.311 (0.99); V2 −0.128 (−0.171, −0.084), retention 0.41; simulation −0.148 = 0.47 of V0.
- Pig / rabbit (descriptive): small, CIs near 0.
Reading: a large part of the drop (≈0.5–0.7) is reproducible by maternal clearance alone, and most of it runs through
clock features that are themselves maternal/zygotic genes (V1 ≈ V0, V2 ≈ 0.4); a residual decrease (≈40% of the
original) remains in both species after all dynamic genes are removed. Neither RESET-ROBUST nor COMPOSITIONAL.
Gate 2b plan frozen (`gz/plan/GATE2B_PLAN_frozen.md`) before any GSE45719/GSE66582 clock score.

## 2026-09-14 — Gate 2b R1 (independent mouse replication, GSE45719 Deng 2014): REPLICATED

Plan `gz/plan/GATE2B_PLAN_frozen.md`; script `gz/src/06_gate2b_R1_deng.py`. 259 cells → 39 embryo pseudobulks
(zygote 4, early 2C 4, mid 2C 6, late 2C 5, 4C 4, 8C 4, 16C 4, early/mid/late blastocyst 3/3/2); reference = zygotes;
clock feature coverage 0.78.
- S1 (late − early 2-cell) −0.137 (95% CI −0.241, −0.026). Most negative pre-blastocyst interval: mid → late 2-cell
  −0.144 (−0.257, −0.030); P_boot within 2-cell 0.73 → REPLICATED by the frozen rule.
- V2 (3,049 maternal/zygotic dynamic genes removed): S1 −0.177 (−0.247, −0.099) — the within-2-cell drop is not reduced
  in this dataset (contrast with Gate 2a, where V2 retained ~0.4 in GSE225056).
- Descriptive: 16-cell → early blastocyst −0.241 (−0.296, −0.188), the largest change overall (outside the tested window).
Limitations: 4–6 embryos per stage; zygote reference n = 4; F1 CAST×B6 hybrid; mm9/2011 annotation counts.
GSE66582 re-quantification started (R2, direction-only readout).

## 2026-09-15 — Gate 2b R2 (GSE66582 Wu 2016, bulk, direction only): CONSISTENT

Script `gz/src/07_gate2b_R2_quant.sh` (kallisto 0.51.1, GRCm39 r112 cDNA+ncRNA; 15/15 runs, mapping 69–90%;
two lookup bugs fixed before any score: grep -P locale, kallisto path; three interrupted downloads resumed) and
`gz/src/08_gate2b_R2_wu.py`. Reference = MII oocytes; clock coverage 0.85; n = 2 per stage.
Stage-to-stage change (primary): MII→zygote −0.128; zygote→early 2C +0.092; early 2C→2C −0.260 (largest, both 2C
replicates −0.28/−0.30 vs early-2C mean −0.03); 2C→4C +0.159; 4C→8C −0.004. ICM mean −0.38 (descriptive).
V2 (3,286 maternal/zygotic dynamic genes removed): early 2C→2C −0.182 (retention 0.70), still the largest → CONSISTENT.
Note: oocyte→zygote is negative here but positive in GSE225056 (mouse, pig, rabbit) — not concordant across datasets.
Summary across mouse datasets: within-2-cell decrease seen in GSE225056 (Gate 1), GSE45719 (R1, REPLICATED) and
GSE66582 (R2, CONSISTENT); after removing dynamic genes it is retained at 0.42 / >1.0 / 0.70.

## 2026-09-15 — Gate 3 (ZGA dependence of the within-2-cell clock decrease) and Gate 4 (zygotic splicing activation): plans frozen, data transfer started

User choice: (가) perturbation test, with (나) splicing kept. Both use the same raw FASTQ (kept on disk).
- Gate 3 plan `gz/plan/GATE3_PLAN_frozen.md` + addendum 1 (download scope). Primary: GSE280522 (A485 ± Dux rescue),
  GSE221985 (maternal Tardbp KO), GSE300734 (maternal Brg1 KO), GSE162345 (SCNT ± α-amanitin); secondary late-2-cell
  contrasts GSE248499 (IVF / SCNT / Kdm4d / Kdm3a), GSE235547 (ICSI / SCNT / SCNT+Obox3; siObox3). Readout: interaction
  of early→late 2-cell drop (perturbed − control), V0 and V2 (control-defined dynamic genes removed); conclusion needs both.
- Gate 4 plan `gz/plan/GATE4_PLAN_frozen.md`: SUPPA2 PSI from kallisto transcript TPM (no Linux aligner on this machine;
  WSL not installed), control ZSA events, splicing progression score, interaction and rescue; clock coupling descriptive.
- 76 runs (206.6 GB) selected: `gz/metadata/perturb/gate3_selected_runs.tsv`; 3 parallel download/kallisto workers
  (`gz/src/10_gate3_download_quant.sh`). Analysis scripts written before any data: `gz/src/11_gate3_analysis.py`.
- pip package "SUPPA" 2.3 installed only a top-level `lib` module without suppa.py; uninstalled; SUPPA2 used from the
  GitHub clone `gz/tools_SUPPA`.
- SUPPA2 (GitHub clone `gz/tools_SUPPA`) patched for current statsmodels: `lib/diff_tools.py` import
  `statsmodels.sandbox.stats.multicomp.multipletests` → `statsmodels.stats.multitest.multipletests` (same function).
  Events from Ensembl GRCm39 r112 GTF (strict boundaries): SE 20,851; A5 9,399; A3 (see file); MX 2,239; RI 4,731;
  AF 32,964; AL 7,429 (`gz/data/annotation/suppa/`). No PSI computed yet.

## 2026-09-15 — Work paused on Windows for transfer to an Ubuntu machine (disk space)

Gate 3 download workers stopped (2 of 76 runs quantified; partial FASTQ not transferred). No Gate 3/4 score computed.
Transfer kit: `gz/tools/pack_for_transfer.sh`, `gz/tools/migrate_paths_ubuntu.sh` (path-only edits, logs itself),
`gz/requirements_gz.txt` (exact versions), handoff `gz/HANDOFF_UBUNTU_KR.md`; CLAUDE.md status note points to it.

## 2026-09-15 — Migration to Ubuntu (shim): path-only edits
Root /home/shim/Downloads/embryogenesis_aging_260915/embryogenesis_Aging_transfer_20260915/embryogenesis_Aging; kallisto /home/shim/anaconda3/envs/gz/bin/kallisto (kallisto, version 0.51.1).
Script checksums changed only by path substitution (C:/shimlab → /home/shim/Downloads/embryogenesis_aging_260915/embryogenesis_Aging_transfer_20260915/embryogenesis_Aging, kallisto.exe → /home/shim/anaconda3/envs/gz/bin/kallisto).
Frozen plan files are unchanged. Post-migration script checksums:
```
0e0de6d9181964ac4a8467d168ce07e557d01cd3dc0647aa3cffb13be1cbc7af  src/07_gate2b_R2_quant.sh
1b94ed6fef039b932ecc6f57246c2acfe5c7a488fa023545bb797dfb5966f728  src/11_gate3_analysis.py
2f20aa1bf2ca73ae82fc093cfb13b5928fcaf3c1cb8b2b8900da33fdf0eb4eee  src/08_gate2b_R2_wu.py
3e8b1d47952dc1773dab5c1005d17dcb94e4bd2fed1392f5a4238afd289d88a4  src/01_gate1_cross_species.py
4866979fbe255a12880b2dcfcbc21c070140afccadecbcee3b53a4f5e5bfaa06  src/09_fig_gate2.py
548b336295949ae9f346404a145ee2def5cae6f99195d55ec84c1c898c9cdaf9  src/02_fig_gate1_trajectories.py
69f1a7bc48068a83f1ce02244fdf26b152d100c33939ec60ff336d588da7f5fb  src/00_validate_tage_py.py
82c7adc4a9c8979b2b8f3e7eaff6994c9850b6510cf8f3d694df9229e2b606ef  src/10_gate3_download_quant.sh
8b69ced258b4bf9b1135e55edc1fbd08112ab9bdad2b540ea48a8f577040f320  src/tage_py.py
bc6d360cbea4d8c634d37e243a1e0d4953dc0abb6182cfb2cdf33561e5a58d84  src/03_gate1alt_zga_timing.py
c929e81e3b3011b74bfcb3b8852e48ef3d3c5fcfabce6659b6988863f755c29f  src/06_gate2b_R1_deng.py
e86e1ed364c12b2b535eb5a101ad2cea5773ee65a01f59affc80ecc93bf04339  src/05_gate2a_composition.py
f6a1dac8b3456d58945678205600389a70cf956042711a6d117f290bbb9b8d15  src/04_fig_gate1alt_zga_vs_clock.py
```

## 2026-09-15 — Ubuntu setup: environment, kallisto index rebuilt, Gate 3 downloads restarted; P4 metadata discrepancy (P4 on hold)

**Environment.** conda env `gz` (python 3.11.16; `requirements_gz.txt` exact pins installed; kallisto 0.51.1 bioconda).
Host `shim`, Linux 5.15.0-139, 24 cores, 62 GB RAM, 948 GB free. tAge port re-validated on this machine
(`src/00_validate_tage_py.py`, rewrites `results/00_tage_py_validation.tsv`): mortality max diff 7.2e-08, chronoage 3.3e-05 → PASS
(same as Windows).

**Transferred kallisto index unusable on Linux.** `analysis/data/processed/mm_full.idx` (Windows build, sha256 fdfd6a89…25bc)
makes both `kallisto inspect` and `kallisto quant` segfault while loading (exit 139). Kept, renamed
`mm_full_windows_built.idx`. Rebuilt from Ensembl release-112 `Mus_musculus.GRCm39.cdna.all.fa.gz` + `ncrna.fa.gz`
(BSD sums 14397/49588 and 35684/8007 = Ensembl CHECKSUMS), concatenated cDNA then ncRNA as `data/raw/mm_cdna_ncrna.fa.gz`
(sha256 88f4c5c1…a30d). The 145,918 FASTA IDs equal `metadata/tx2gene_full.tsv` in content and order. `kallisto index -t 8`
(log `logs/kallisto_index_full_ubuntu.log`): 750,479 contigs, 117,049,481 k-mers, 145,918 targets = the original Windows build
log. New `mm_full.idx` sha256 22c45310…c29b. Check on 200,000 read pairs of SRR30827777: 87.8% pseudoaligned, fragment length
73.06 (Windows full run: 88.0%, 73.3).

**Windows-quantified Gate 3 runs moved, not deleted.** The handoff said 2 of 76; there were 3 DONE (SRR30827777, SRR16201705,
SRR16201707) plus an unfinished SRR16201702. Moved to `gz/data/perturb/<gse>/quant_windows_index/` so that every Gate 3 library
is quantified with one index binary and its FASTQ is kept for Gate 4. These 3 runs are downloaded first; their Linux-index vs
Windows-index abundances will be compared (also bears on Gate 2b R2, quantified with the Windows index). No score computed.

**Download script (new file, `10_` kept).** `src/10u_gate3_download_quant.sh` + `src/10u_make_worklists.py`: same 76 runs, same
kallisto parameters; each FASTQ checked against ENA size and md5 (`metadata/perturb/gate3_ena_md5.tsv`, fetched 2026-09-15;
layout, read counts and URLs equal the frozen lists; 206.6 GB); files failing md5 are deleted and re-fetched; interrupted
downloads resume over repeated passes. Order (execution only): the 3 re-quantified runs, then P2 GSE221985, P1 GSE280522,
P3 GSE300734, S2 GSE235547, S1 GSE248499, P4 GSE162345. Worker lists `metadata/perturb/gate3_ubuntu_worker{0,1,2}.tsv`.
Network: wireless link only (wired interface unused); ENA, Ensembl and AWS all give about 1.0–1.3 MB/s total, parallel
connections share it → expected transfer time about 44 h.

**P4 GSE162345 — repository metadata differ from the frozen plan (analysis on hold, no score computed).**
The plan reads P4 as "SCNT E2C, L2C × ± α-amanitin; PE36". GEO (series + GSM5610457–68) says: C2C12 cells were transplanted
(fused) into *early (21 hpi) or late (30 hpi) 2-cell embryos*, cultured 24 h in 0.1 µg/ml demecolcine ± α-amanitin, and
collected at 45 hpi (early2cell_NT) or 54 hpi (late2cell_NT); pools of 5 embryos; strains C57BL/6, DBA/2, C3H. So "early" and
"late" name the recipient stage at transfer, not an early→late 2-cell sampling of the same embryo type; α-amanitin starts
before (early) or after (late) major ZGA onset; enucleation is not stated. The early→late drop D and the interaction I defined
in the plan therefore do not measure the intended window in this dataset. Also: deposited FASTQ read lengths are about 51 nt
(R1) and 25 nt (R2); R2 is shorter than k = 31, so paired-mode kallisto uses R1 only and cannot estimate fragment length
(Windows runs reported 0.0). Per CLAUDE.md, P4 is held; its role (drop, keep as descriptive, or a new single-stage contrast)
and the resulting primary decision count must be settled in an addendum before any Gate 3 score. FASTQ are still downloaded
(last in order, 10.3 GB) so either option remains possible. `src/11_gate3_analysis.py` currently maps early2cell_NT → E2C
and must not be run on P4 until that addendum exists.

### 2026-09-15 03:00 — Gate 3 workers running; cross-index concordance; DDBJ not used

Workers (`src/10u_gate3_download_quant.sh`, 3 × `gate3_ubuntu_worker{0,1,2}.tsv`) started 02:42. First 3 runs: all 6 FASTQ
md5 = ENA; quant ok. Measured total throughput over 3 workers 2.62 MB/s (15-min window) and 2.78 MB/s (10-min window) → about
21 h for the remaining ≈201 GB. The earlier single-connection estimate of ≈44 h is superseded.

**Linux-index vs Windows-index quant, same FASTQ (`src/12u_quant_concordance.py`, `results/ubuntu_quant_concordance.tsv`).**
SRR16201705, SRR16201707: transcript est_counts identical (no fragment-length estimate in either run, see P4 note).
SRR30827777: fragment length 73.2 (Linux) vs 73.3 (Windows); transcript max |Δcount| 3.02, max |Δeff_length| 0.12; gene-level
max |log2 Δ| 0.0064 over 10,291 genes with ≥10 counts; 0 genes > 0.1. Conclusion: the rebuilt index reproduces the Windows index;
Gate 2b R2 (Windows index) remains comparable. No clock score involved.

**DDBJ SRA-lite mirror tested and not adopted.** Domestic mirror 10 MB/s and DDBJ about 20 MB/s over 3 connections show that the
international route to ENA is the limit, but DDBJ lacks most recent runs (45 of the 67 probed returned 404, including all
GSE280522), and SRA-lite has simplified quality values and no ENA md5. Mixing sources is not worth the provenance cost; all
Gate 3 FASTQ come from ENA.

### 2026-09-15 20:15 — Gate 3 workers stopped by a system hang at 11:17; state checked, workers restarted

**Cause (kernel journal, previous boot).** 2026-09-13 00:02 kernel oops in memory management (`BUG: unable to handle page
fault`, RIP `__mod_lruvec_page_state`, process `R`); the kernel kept running tainted (D). 2026-09-15 11:17:20 `kswapd0` soft
lockup on a page-table spinlock (RBX = the R08 page address of the 09-13 oops), then curl, gnome-shell and other tasks locked;
last journal line 12:41; reboot 17:29. Not a power cut and not caused by the pipeline. Separately, the Intel Wi-Fi firmware
crashed and restarted 150 times in that boot (slower downloads).

**State at stop (no score computed).** 22/76 runs quantified with the Ubuntu index (GSE280522 9/24, GSE221985 7/8,
GSE300734 1/12, GSE235547 3/10, GSE248499 0/10, GSE162345 2/12); 26 partial, 27 not started; 88.9 of 206.6 GB on disk.
Integrity: all 22 DONE have 145,918-row `abundance.tsv` and a readable `run_info.json`; md5 recomputed for the 8 FASTQ verified
after 10:40 = ENA; no FASTQ has a zero-filled tail; the only non-gzip file was below.

**Two problems found.**
- `GSE300734/fastq/SRR34172935_1.fastq.gz` was a 199-byte "403 Forbidden" HTML page (10:19). `curl` in
  `src/10u_gate3_download_quant.sh` has no `-f`, so an HTTP error body is saved; a resumed download would append the real file
  after it and fail md5 at full size. File deleted; ENA now answers 200 with the listed size. Script not changed (md5 still
  rejects such files; cost is one re-download).
- SRR34172920 (GSE300734): both FASTQ md5 = ENA, but kallisto died at about 6 M reads (kernel: general protection fault in
  kallisto, 09:54; `quant.log` "quant failed"); this was while the kernel was tainted. The worker retries it (both `.ok` present).

**Restart.** 20:13:16, same command as the handoff except worker logs are appended (`>>`) instead of overwritten, so the
worker0 segfault line is kept. Kallisto parameters, index and worklists unchanged. Measured total throughput over 41 s after
restart about 2.35 MB/s (≈118 GB left).

**SRR34172920 re-test after reboot (outside the project, scratchpad only; same index, same `quant -t 4` call).** Exit 0,
4 min 49 s, max RSS 1.3 GB; 29,053,854 pairs processed, 58.8% pseudoaligned (the failed run showed 58.9% at 6 M pairs).
The 09:54 crash did not reproduce on the same md5-verified FASTQ, consistent with the unstable kernel rather than the data.
This output is not used; the project copy is produced by the worker, which re-runs the run itself.

### 2026-09-16 12:21 — Gate 3 downloads: 71/76 quantified; transfer corruption rate and an at-rest integrity check

Worker1 finished all its runs at 09:25 (pass 9). Worker0 and worker2 are on pass 7 with five runs left (SRR26923894,
SRR26923891, SRR26923899 of GSE248499; SRR34172919 of GSE300734; SRR16201701 of GSE162345), about 11 GB.

**Transfer corruption.** 32 downloads reached the declared ENA size but failed md5 and were deleted and re-fetched
(4 on 09-15, 28 on 09-16; 49.4 GB re-downloaded; 144 files passed, so about 18% of completed transfers). Every affected run
passed on a later attempt. Intel Wi-Fi firmware crashed and restarted 136 times in the current boot (150 in the previous one);
the machine has no wired connection in use. No kernel oops, no I/O error and no EDAC event this boot.

**At-rest check (rules out silent memory or disk corruption).** md5 recomputed for all 144 FASTQ that carry an `.ok`
file: 144 checked, 0 mismatches. Stored data therefore match ENA; the corruption happens in transfer, not after writing.
Together with the SRR34172920 re-test this leaves no evidence that the 09-13/09-15 kernel fault affected any quantified run.

### 2026-09-16 17:55 — Gate 3 downloads and quantification complete: 76/76 (no score computed)

All three workers exited with "all runs DONE": worker1 09:25 (pass 9), worker2 17:34 (pass 12), worker0 17:55 (pass 18).
76/76 runs quantified with the Ubuntu-rebuilt `mm_full.idx` (kallisto 0.51.1, GRCm39 r112 cDNA+ncRNA); 152/152 FASTQ carry
an `.ok` file (md5 = ENA), 206.6 GB kept on disk for Gate 4; 753 GB free. No run was dropped.

**Transfer cost.** 43 completed downloads failed md5 and were re-fetched (96.9 GB wasted; 152 passed). Worst:
SRR26923894_2 (4.5 GB) 8 attempts, SRR34172919_1 (5.0 GB) 6 attempts — both large files on a Wi-Fi link whose firmware
restarted 136+ times per boot, so long transfers are interrupted and resumed repeatedly. Data at rest were verified
(previous entry): the corruption is in transfer only, and every file that entered the analysis matches the ENA md5.

**Library QC against the frozen plan (`p_pseudoaligned` < 30% is a prespecified exclusion).**

| dataset | n | reads (M) | pseudoaligned % min/median/max |
|---|---|---|---|
| GSE280522 (P1) | 24 | 10.6–35.9 | 56.7 / 68.5 / 87.8 |
| GSE221985 (P2) | 8 | 8.6–28.4 | 40.2 / 65.4 / 69.9 |
| GSE300734 (P3) | 12 | 25.3–73.2 | 21.4 / 55.2 / 63.2 |
| GSE162345 (P4, on hold) | 12 | 13.2–21.3 | 77.0 / 79.2 / 81.1 |
| GSE248499 (S1) | 10 | 13.7–49.0 | 62.7 / 69.9 / 72.8 |
| GSE235547 (S2) | 10 | 13.3–38.5 | 66.4 / 72.5 / 76.9 |

One library fails the rule: **SRR34172919 (GSE300734 WT_E2C_r2), 21.4% pseudoaligned, 6.7% unique, 73.2 M reads**
(other GSE300734 libraries 36.8–63.2%; same 150 nt read length; fragment length 155.6). Both mates match the ENA md5, so
this is the deposited library, not a download artefact. Under the prespecified QC it is excluded, which leaves P3 control
(WT) early-2-cell at n = 2 of 3 — to be applied, not decided, when `src/11_gate3_analysis.py` runs. The detected-genes
criterion is evaluated at that point as well.

**Still blocking Gate 3 analysis:** the P4 GSE162345 addendum (metadata discrepancy, entry of 2026-09-15) must be written
before `src/11_gate3_analysis.py` is run. Gate 4 needs the SUPPA events regenerated on this machine
(`data/annotation/suppa/` was not transferred; the GTF is in place).

### 2026-09-16 — Gate 3 addendum 2 (P4 GSE162345 primary → secondary) written, approved and applied to the code

Written and approved before any Gate 3 score; `src/11_gate3_analysis.py` had never been run, and no result from any
Gate 3 dataset had been computed or seen. Plan file `gz/plan/GATE3_PLAN_addendum2.md`.

**Decision.** GEO (series, GSM5610457–68) describes transfer of C2C12 cells into 2-cell embryos at 21 hpi ("early") or
30 hpi ("late"), 24 h culture in demecolcine ± α-amanitin, collection at 45 or 54 hpi — so "early"/"late" are recipient
stages at transfer, not a stage pair, and the frozen `D`/`I` do not measure the intended window in this dataset.
P4 therefore moves from primary to secondary and is read as a window-matched transcription-block contrast:
`T45 = tAge(α-amanitin 45 hpi) − tAge(control 45 hpi)` and `T54` likewise, V0 only, 3 libraries per arm, bootstrap as
for the other secondaries, direction prespecified as `T > 0` if the decrease depends on zygotic transcription. Any
45 → 54 hpi difference is explicitly not computed. User approved both the move and keeping the two windows separate.

**Effect on the decision rules (unchanged rules, changed counts).** Primary = P1 GSE280522, P2 GSE221985, P3 GSE300734.
The frozen "fewer than four informative" clause then requires `I > 0` in all informative primary datasets (≥ 2
informative) plus `R < 0` in P1 for ZGA-DEPENDENT — a stricter bar than 3 of 4. NOT-ZGA-DEPENDENT and MIXED unchanged.
P4 contributes to none of these counts.

**Quantification of P4 kept as quantified** (paired call, R2 = 24 nt < k, R1 effectively used, fragment length 0.0 in all
12 libraries, pseudoalignment 77.0–81.1%): the bias is uniform within the dataset and the contrasts are within-dataset.
A single-end sensitivity re-quantification is defined in the addendum but not run.

**Code changes** (`src/11_gate3_analysis.py`, before any score): docstring; `annotate()` GSE162345 branch returns
window `NT45hpi`/`NT54hpi` and window-tagged arms; GSE162345 moved to the secondary list; removed from the primary
verdict loop; secondary rows now carry the real stage; output renamed `results/gate3_secondary_L2C.tsv` →
`results/gate3_secondary_single_stage.tsv`.

**Pre-existing fault found and fixed in the same expression (no score had ever been computed with it).** Reference
selection read `'control_ICSI' if arm.startswith('SCNT') else …`, so all three S1 GSE248499 SCNT arms pointed at
`control_ICSI`, an arm that exists only in GSE235547; the reference would have been empty and the S1 differences NaN.
References are now per dataset — S1: IVF (`control`) for `SCNT_control`, `SCNT_Kdm4d`, `SCNT_Kdm3a`; S2 unchanged;
P4 same-window control — and an empty reference raises instead of returning NaN. Verified without scoring: the arm
assignment gives S1 3 contrasts, S2 3, P4 2, every reference non-empty (S1 IVF n = 4, S2 n = 2, P4 n = 3).

Checksums after the change:
```
fe6e838cfa2400495d97f03fd677c235be11ef8094e52c98945657217ace17f1  plan/GATE3_PLAN_addendum2.md
d2f8d2645f8096378b9c3a8f957066a0b602f964451b5cab55783186a6e7a2c7  src/11_gate3_analysis.py
9a0fca54be18e6ad44fde1b76bb7067b685ee9e51a6e589913758ac84dfab5cc  plan/GATE3_PLAN_frozen.md (unchanged)
60b33216586549b3b833af728dac6b8a4bf7c25ab64d8ec83d50674aac4e4fa3  plan/GATE3_PLAN_addendum1.md (unchanged)
```

### 2026-09-16 23:53 — Gate 3 run: OVERALL **MIXED** (V0 ZGA-DEPENDENT, V2 MIXED); rescue arm does not restore

First and only run of `src/11_gate3_analysis.py` (log `logs/11_gate3_analysis.log`, 9 s). Before it produced any number a
third pre-existing fault was fixed: the verdict block read the column `drop` as an attribute (`ctrl.drop`), which returns
the DataFrame method, raising `AttributeError`; changed to `ctrl['drop']` in three places. Script checksum after the run:
`1168c609244ecb97d53bb531d7429ee6e3d92a736bbcd30f8b3242db7071baac`.

**Library QC (prespecified, applied automatically).** 3 of 76 libraries excluded: SRR34172919 (GSE300734, control E2C,
21.4% pseudoaligned — the run flagged in the 17:55 entry), SRR30827777 (GSE280522, A485 L2C, 22,471 detected genes vs
cutoff 22,524) and SRR26923904 (GSE248499, IVF control, 26,754 vs 27,071). All three exclusions are by the frozen rule,
not by choice. Control E2C in P3 is therefore n = 2.

**Primary (drop D per arm; interaction I = D(perturbed) − D(control); bootstrap 2000, seed 20260918).**

| dataset | variant | D control | D perturbed | I [95% CI] |
|---|---|---|---|---|
| P1 GSE280522 (A485) | V0 | −0.2375 | −0.0566 | **+0.1809 [0.1082, 0.2526]** |
| P1 GSE280522 (A485) | V2 | −0.2240 | −0.0430 | **+0.1810 [0.1117, 0.2510]** |
| P2 GSE221985 (Tardbp matKO, n = 2) | V0 | −0.1003 | −0.0588 | +0.0415 [−0.0534, 0.1364] |
| P2 GSE221985 | V2 | −0.0722 | −0.0871 | −0.0148 [−0.1196, 0.0907] |
| P3 GSE300734 (Brg1 matKO) | V0 | −0.1714 | −0.0623 | +0.1091 [−0.0108, 0.2228] |
| P3 GSE300734 | V2 | −0.1508 | −0.0460 | +0.1048 [−0.0110, 0.2169] |

All three primary datasets are informative (control D < 0 in both variants). V0: I > 0 in 3 of 3 and R < 0 →
ZGA-DEPENDENT. V2: P2 flips sign (I = −0.0148), so 2 of 3 → MIXED. The plan requires V0 and V2 to agree, so
**OVERALL = MIXED**.

**The P1 rescue arm does not restore the decrease.** R = D(A485+Dux) − D(A485) = −0.0071 [−0.0629, 0.0572] (V0) and
−0.0060 [−0.0633, 0.0583] (V2): A485+Dux (−0.0637) is indistinguishable from A485 (−0.0566), both far from control
(−0.2375). The rule `R < 0` is satisfied by sign alone and should not be read as restoration; recorded here so that the
verdict field `rescue_restores: True` is not quoted without this line.

**Secondary single-stage contrasts (V0 only).**

| dataset | arm | reference | diff [95% CI] |
|---|---|---|---|
| P4 GSE162345 | α-amanitin 45 hpi | control 45 hpi | **+0.1550 [0.0674, 0.2426]** |
| P4 GSE162345 | α-amanitin 54 hpi | control 54 hpi | +0.0971 [−0.0527, 0.2190] |
| S2 GSE235547 | siObox3 | siCtrl | **+0.1091 [0.0532, 0.1650]** |
| S2 GSE235547 | SCNT-Obox3 | ICSI | −0.0558 [−0.1211, 0.0095] |
| S2 GSE235547 | SCNT-EGFP | ICSI | −0.0487 [−0.1156, 0.0182] |
| S1 GSE248499 | SCNT control | IVF | −0.1367 [−0.1927, −0.0807] |
| S1 GSE248499 | SCNT + Kdm3a | IVF | −0.1682 [−0.1901, −0.1376] |
| S1 GSE248499 | SCNT + Kdm4d | IVF | −0.2461 [−0.3112, −0.1809] |

P4 (α-amanitin) and S2 (siObox3) go in the prespecified direction: blocking transcription or knocking down a ZGA
activator leaves a higher tAge at the same collection point. S1 goes the other way — SCNT late-2-cell libraries sit
*below* IVF, most strongly with Kdm4d — so the prespecified secondary direction does not hold there; reported, not
explained here.

**What may be said now.** In the datasets analysed, the within-2-cell tAge decrease is smaller when minor ZGA is blocked
pharmacologically (P1, both variants, CI excluding 0) and, in the same direction but with CIs including 0, when Brg1 is
removed maternally (P3). P2 does not hold its direction across variants. Dux co-expression did not restore the decrease.
No causal verb is licensed, and the overall Gate 3 verdict is MIXED, not ZGA-DEPENDENT.

## 2026-09-17 — POST HOC decomposition of the clock decrease (descriptive; Gate 3 unchanged)

Plan frozen before running (`gz/plan/POSTHOC_gate3_contribution_frozen.md`), script
`gz/src/13_posthoc_gate3_contribution.py`, log `logs/13_posthoc_gate3_contribution.log`, prose
`gz/results/posthoc_gate3_contribution.md`. Labelled POST HOC wherever reported; no p-value, no new dataset, no gate
re-opened, no change to the MIXED verdict. Motivation recorded at the time: a single score moving down is not yet a
biological statement.

The clock pipeline is linear after imputation and centring (`ElasticNet`, 1,839 non-zero coefficients), so with
`species=None` the group difference decomposes exactly: `D = Σ_g coef_g · (mean x_g(L2C) − mean x_g(E2C))`. The script
asserts that the per-gene terms reproduce `results/gate3_arm_drops.tsv` to 1e-9; they do (P1 −0.2375, P3 −0.1714).

1. **The decrease is a residual, not a clean programme.** P1 negative contributions sum to −1.4030, positive to
   +1.1655; the reported drop is the 17% difference. 161 genes are needed for half of the total absolute movement.
2. **Activation and clearance are of similar size.** Genes induced in the window with negative coefficients contribute
   −0.746; genes cleared with positive coefficients −0.650; the opposing cells +0.531 and +0.626. P3 agrees.
3. **Canonical maternal/zygotic genes are not what the clock reads here.** 99.5% of P1's D comes from genes labelled
   neither (1,346 genes); 8 maternal (0.2%), 1 zygotic (0.3%). P3: 92.2% / 8.0% / −0.2%. Consistent with V2 leaving the
   interaction essentially unchanged (+0.1809 → +0.1810).
4. **Composition reproduces across two independent datasets.** Per-gene contributions correlate between P1 (A485) and
   P3 (maternal Brg1 loss): Pearson 0.739, Spearman 0.606 over 1,839 genes; 28 of the 50 largest down-contributors are
   shared. Different labs, different perturbations.
5. **What A485 removes is concentrated on the same genes.** corr(c_control, Δ) = −0.565; the ten largest control
   down-contributors carry 63% of the interaction I. Strongly induced negative-coefficient genes (Klf9 log2FC +4.48,
   Neto2 +2.92, Pi4k2a +3.03, Gpatch4 +2.89, Psmb5 +3.81) lose nearly all of their contribution under A485.

Limits stated in the plan and repeated in the output: the coefficients are from an adult multi-species multi-tissue
clock, so a large contribution means the clock weights that gene and the gene moves in this window — not that the gene
has an embryonic role. No enrichment test was run; no gene is promoted to a claim.

```
5dbe30040cebcf228fa3ce33e38aefc12b8f575bda9770a0cae8d90c1b21e098  plan/POSTHOC_gate3_contribution_frozen.md
7a0499f41f2361abeedfbd5c29985137acf0de15719c9fb65b7825fb7385e03a  src/13_posthoc_gate3_contribution.py
c71d8b5c455db5693da990cab9f7f35d8d4262a28e052cb7bf0d8a7129eebce6  results/posthoc_gate3_contribution_genes.tsv
61902eebefda6eab0697a8a95ab20d7f213ab58a75eba42f30c5bbe8abbcc683  results/posthoc_gate3_contribution_genes_P3.tsv
cf2d0fdf3865127fd09963d4f8846452027ad0ef3b700c6c468ceea4759be060  results/posthoc_gate3_contribution_summary.tsv
424c7f485a8305a7638d1593310bf16f60a2ac2978af9b69d34d6c014ea78b4c  results/posthoc_gate3_contribution.md
```

## 2026-09-17 — POST HOC: the A485+Dux rescue arm worked; the scalar clock difference did not follow

Plan frozen before running (`gz/plan/POSTHOC_rescue_frozen.md`, which named two separable explanations), script
`gz/src/14_posthoc_rescue.py`, log `logs/14_posthoc_rescue.log`, prose `gz/results/posthoc_rescue.md`. Descriptive;
Gate 3 stays MIXED and `results/gate3_arm_drops.tsv` was not recomputed.

**Zygotic transcription was restored (explanation (a) excluded).** Zygotic set (13 genes, control-arm definition),
mean log2(CPM+1) at L2C: control 3.565, A485 1.096, A485+Dux 2.763. Differences with bootstrap CI (2000, seed
20260918): A485 − control −2.469 [−2.723, −2.214]; A485+Dux − A485 +1.667 [1.139, 2.037]; A485+Dux − control
−0.802 [−1.335, −0.358]. Maternal set moves oppositely and less. Markers (fixed list written before looking):
Zscan4f 6.71 → 9.15, Zscan4c 5.28 → 8.16, Zscan4d 5.25 → 7.70 from A485 to A485+Dux, each above the control L2C value.

**The clock's down-contributors were largely restored too.** Sum over the 666 genes with negative control contribution:
control −1.4030, A485 −0.7592, A485+Dux −1.1867 (85%). Over the top 20: −0.2955 / −0.1304 / −0.2210. Per-gene profile
correlation with control: 0.647 (A485), 0.833 (A485+Dux).

**Why the score stayed flat (unplanned follow-up, computed after the two planned items, descriptive).** Per-arm flows:
control −1.4030 / +1.1655 = −0.2375; A485 −1.1529 / +1.0964 = −0.0566; A485+Dux −1.3803 / +1.3166 = −0.0637. The
rescue arm restores the downward flow to 98% of control while the upward flow reaches 113%, so the residual reported by
the clock does not move. Largest upward gains vs control: Cd74 +0.0241, Parp3 +0.0233, Cst7 +0.0137, Icam1 +0.0119,
Upp1 +0.0116, Gpc4 +0.0115, Cdkn1a +0.0111, Lgals3 +0.0089. These genes were selected by clock contribution, not by any
gene-set test; the stress/inflammation reading of that list is an observation for a prespecified test, not a finding.

**Consequence for how the rescue is reported.** R < 0 in the frozen rule was satisfied by sign only (previous entry);
this entry shows the arm is not uninformative but that a single scalar difference can stay flat while both component
flows are restored or exceeded. Statements about the rescue arm are to be made on the decomposition, not on R.

```
45ab355b7fedd62db294b06183229b794ff0d530393b809b18e091c545d76cb9  plan/POSTHOC_rescue_frozen.md
9bde8a69186f317feac4e3c9d46a4152f68b12ad39ddb35da3ed56ad7db9da7d  src/14_posthoc_rescue.py
7a4017583330897f2eb0bfefde3e33e255a15ebdfad9fe697932948804c868cb  results/posthoc_rescue_contributions.tsv
73d7a5d11c5a6e0e2f4e0afc02e979035ee9870fe61b77dc0db0db517ebccfc9  results/posthoc_rescue_markers.tsv
387feeb68cae7ccf8e24052500706d1faa78af61ba37fd1f31e196a17aa0bbdc  results/posthoc_rescue_programme_scores.tsv
6b00eafd125508df3ee3884050d8590e020cc945232a188dd8affea30aab5b68  results/posthoc_rescue.md
```

### 2026-09-17 — Provenance check: the Gate 3 data are GEO series; ENA was only the download route (NCBI cross-check passes)

Question raised while preparing the manuscript: were these GEO datasets, given that the FASTQ came from ENA. They are.
Every Gate 3 run belongs to a GEO series (GSE) whose samples (GSM) are SRA experiments (SRX) with runs (SRR); ENA mirrors
the same runs under the same accessions. Verified rather than assumed, by NCBI E-utilities `runinfo` for all 76 runs
(fetched 2026-09-17, four batches by run accession; GSE numbers are not indexed in SRA for five of the six series, so
the query was by SRR):

- NCBI `spots` = ENA `read_count` for **76/76** runs.
- NCBI `spots` = kallisto `n_processed` for **76/76** runs — every submitted read pair entered the quantification.
- `LibraryLayout` agrees for 76/76.
- 76 runs map 1:1 to 76 GSM samples (SRX → GSM taken from the stored GEO SOFT records; NCBI's `SampleName` holds the
  library name rather than the GSM for GSE280522, so SOFT is the authority for that field).

New file `gz/metadata/perturb/gate3_accessions_geo_sra.tsv` (sha256 `1ad17717c320d40ab7d1dbdedf55368aa7ac45203c0694ba2f5b533082176654`):
gse, gsm, gsm_title, srx, srr, biosample, layout, spots_ncbi, avg_read_len — the table for the data-availability section.

**One difference worth a methods line.** NCBI's default delivery object for all 76 runs is SRA-Lite
(`sra-pub-zq-*`, `.lite.1`), which carries simplified quality strings; the ENA FASTQ we used carry the submitted
quality strings and were each verified against ENA's md5. Read identity and counts are the same either way, and
kallisto pseudoalignment does not read quality values, so no result depends on this — but anyone re-running from
`fasterq-dump` defaults will get the Lite qualities, and the methods should say which copy was used.

### 2026-09-17 — Gate 4 set-up before any PSI: events regenerated, implementation notes frozen

SUPPA2 events regenerated on this machine from `Mus_musculus.GRCm39.112.gtf.gz` (sha256 a1ad41014c12…223c) with
`generateEvents -f ioe -e SE SS MX RI FL -b S` (log `logs/gate4_generateEvents.log`, 9 s). Counts equal the Windows run of
2026-09-15: SE 20,851; A5 9,399; A3 10,913; MX 2,239; RI 4,731; AF 32,964; AL 7,429. The statsmodels import patch in
`gz/tools_SUPPA/lib/diff_tools.py` is present in the transferred clone.

`gz/plan/GATE4_IMPLEMENTATION_NOTES.md` written before any PSI: SUPPA defaults named, computation order, merged ioe,
TPM format, literal application of the plan's decision bullets. No plan change. It also corrects the running record:
the frozen Gate 4 plan already excludes P4 GSE162345, so the note in `SUMMARY_gz_*.md` and the review page that P4's
Gate 4 status was open was wrong. Read-length condition (≥ 95 nt) verified on the data: P1 mixed 72/150 nt (mean ≈ 110),
P2 and P3 150 nt.

```
0f51aae2fe70455076bb8199b3cac109105244eb75f02fd897bc7d027440b64d  plan/GATE4_PLAN_frozen.md (unchanged since 2026-09-15)
7a0bff2b32a3baefe7f3541aedccacb3e66a53d56097940df851d57288daeb9b  plan/GATE4_IMPLEMENTATION_NOTES.md
fc12a1a756b427333bc85372328d0bc111bee540a3c98bcb1a0ba7d870bfc677  src/15_gate4_analysis.py
```

### 2026-09-17 13:10 — Gate 4 run: **ZSA ZGA-DEPENDENT** by the frozen rule; a selection bias in the score is flagged

First and only run of `src/15_gate4_analysis.py` (log `logs/15_gate4_analysis.log`, 2 min 32 s; script unchanged since the
pre-run checksum fc12a1a7…c677). Libraries = Gate 3 QC set: P1 23, P2 8, P3 11. Working files `gz/data/gate4/` (201 MB).

**ZSA magnitude (control arm).** Filtered events / control ZSA events (|ΔPSI| ≥ 0.10 and diffSplice p < 0.05):
P1 GSE280522 13,715 / 511 (3.7%); P2 GSE221985 13,053 / 1,376 (10.5%); P3 GSE300734 5,732 / 77 (1.3%).
≥ 100 in 2 of 3 → ZSA PRESENT. All seven event types contribute; AF and SE carry the most events.

**Progression (control ≡ 1 by construction; bootstrap 2000, seed 20260919).**

| dataset | P control | P perturbed | interaction [95% CI] |
|---|---|---|---|
| P1 A485 | 1.000 | 0.557 | **−0.443 [−0.669, −0.220]** |
| P1 A485+Dux | 1.000 | 0.787 | −0.213 [−0.504, 0.072] |
| P2 Tardbp matKO | 1.000 | 0.579 | **−0.421 [−0.488, −0.355]** |
| P3 Brg1 matKO | 1.000 | 0.549 | **−0.451 [−0.660, −0.232]** |

P1 rescue R = P(A485+Dux) − P(A485) = +0.230 [−0.016, 0.473]. Interaction < 0 in 3 of 3 and R > 0 → by the plan's
bullets, **ZSA ZGA-DEPENDENT**. The rescue CI touches 0; the rule uses the point estimate.

**Clock coupling (descriptive, plan readout 3).** Spearman between arm × stage group means of progression and Gate 3 V0
tAge: P1 −0.943 (6 groups), P2 −0.400 (4), P3 −1.000 (4), pooled −0.603 (14) — more splicing progression goes with a
lower clock value. Across the four perturbation interactions, Spearman(progression I, clock I) = 0.000: P2 shows as
much splicing attenuation (−0.421) as P1 and P3 while its clock interaction is +0.042, so the size of the two effects does
not track across perturbations.

**Threat to interpretation, not anticipated in the frozen plan (recorded before any follow-up is run).** ZSA events are
selected on the control libraries (|ΔPSI| ≥ 0.10, p < 0.05) and the score is scaled on the same libraries, so the
control arm's P = 1 is in-sample while every perturbed arm is scored out-of-sample. Selection on noisy control ΔPSI
inflates the control progression relative to any other arm (winner's curse): an untreated replicate group scored the
same way would also fall below 1. Part of each negative interaction may therefore be produced by the construction. The
A485+Dux arm (0.787) bounds the size of this effect only loosely. The verdict stands as the rule defines it, but it is
not to be reported as a splicing finding until a held-out control estimate exists. A cross-fitted check (select events
on part of the control libraries, score the rest) is feasible in P1 (4 + 4 control libraries) and not in P2 (2 + 2);
it is post hoc and needs its own frozen plan.

```
8cbca1bffb28cbee9d9ae44db25a98065104e1184cecd0fb6aa66d2e8a68ce6a  results/gate4_arm_progression.tsv
27770daec99bf778e3179e374af66c30ce96795f8cea9c18ac308a5323514c9e  results/gate4_events_filtered_GSE280522.tsv
c8c6db136e0a03f068157e911872f24485254e53255a2c280b2f2b230272e5b6  results/gate4_events_filtered_GSE221985.tsv
ed00a0758b49c472f94ba0aa8fe1424506eefe18700573acf4ffc365653dc9fb  results/gate4_events_filtered_GSE300734.tsv
bd926e8467fbd749ded7004b7c43806bd7ac60a3d94b8ea522d797c8cf63c373  results/gate4_progression.tsv
3d7de375ab9659e633e66665348311b38965ffed07605b853f6f6d5183820e2e  results/gate4_zsa_summary.tsv
80d6ad2a195e3c1ef774a59677ddab84eb9daa29f7818b5b6b812d93f05ba08e  results/gate4_decision.md
```

### 2026-09-17 — POST HOC Gate 4 cross-fitted check: plan frozen before running

Plan `gz/plan/POSTHOC_gate4_crossfit_frozen.md` (16 leave-one-E2C-and-one-L2C-out folds in P1; reproduction check
first; reading SURVIVES / ATTENUATED BUT INFLATED / CONSTRUCTION fixed with the −0.20 and 15-of-16 thresholds). P2 and P3
cannot be split and stay "not cross-checked". Script `gz/src/16_posthoc_gate4_crossfit.py`. Pre-run checksums:
```
ae70df8563ea735be7f7fec8b22a69bd6daec708dea4049c835ac2b8ba6c5976  plan/POSTHOC_gate4_crossfit_frozen.md
7d0f55712adfdc3ca136e14d2c77aa1e3ecfd5b8c144cd56a4dc56b8b12afa15  src/16_posthoc_gate4_crossfit.py
```

### 2026-09-17 — POST HOC Gate 4 cross-fitted check: **ATTENUATED BUT INFLATED**

**Run history.** Attempt 1 stopped at the reproduction check with a pandas error: the Gate 4 reference table holds
`P control` once per dataset and the lookup did not filter to GSE280522 (log kept as
`logs/16_posthoc_gate4_crossfit_attempt1_checkbug.log`; the computed check values already equalled Gate 4 and no fold
had run). One line fixed (filter `gse == GSE280522`); script checksum after the fix
`87510c817a29c3a5ea2d47badd0827d5e54fc3eca677d4214ef939da109b3fbe` (pre-run `7d0f5571…fa15`). Attempt 2 ran to the end
(log `logs/16_posthoc_gate4_crossfit.log`, ≈ 25 min). Plan unchanged (`ae70df85…5976`).

**Reproduction check passed**: with no hold-out the code returns P control 1.000000, A485 0.556809, A485+Dux 0.786949
and the identical 511-event set.

**16 folds (mean; fold spread 2.5–97.5%, not a CI — folds share libraries).**

| quantity | mean | fold spread |
|---|---|---|
| fold ZSA events | 571 | 429 – 808 (76% overlap with the Gate 4 set) |
| P control, in-sample | 1.000 | by construction |
| **P control, held out** | **0.790** | 0.380 – 1.072 |
| P A485 | 0.517 | 0.468 – 0.575 |
| P A485+Dux | 0.734 | 0.626 – 0.822 |
| **I A485 (vs held-out control)** | **−0.273** | −0.546 – +0.100 |
| I A485+Dux | −0.056 | −0.308 – +0.246 |
| R = P(A485+Dux) − P(A485) | +0.217 | 0.146 – 0.251 |

In-sample inflation of the control arm = 0.210, i.e. about a fifth of the Gate 4 control value, and roughly half of the
Gate 4 A485 interaction (−0.443 → −0.273). I A485 is negative in 14 of 16 folds. By the frozen reading
(≤ −0.20 **and** ≥ 15/16 negative for SURVIVES) the outcome is **ATTENUATED BUT INFLATED**: the direction may be
reported, the Gate 4 magnitude may not. It missed SURVIVES on the fold count (14 vs 15); that is recorded as the result,
not re-read.

Held-out control progression depends on which libraries are held out (by held-out L2C: SRR30827779 0.563, SRR30827780
0.929, SRR30827781 0.913, SRR30827782 0.757; by held-out E2C: SRR30827763 0.572, others 0.750–0.924). With one library per
stage held out, single-library heterogeneity dominates the fold spread.

**Rescue.** R does not involve the control arm and is unaffected by the selection bias; it is positive in 16 of 16
folds (minimum +0.144) and the A485+Dux arm sits close to the held-out control level (I A485+Dux −0.056; below held-out
control in 10 of 16 folds). The Gate 4 bootstrap interval for R, +0.230 [−0.016, 0.473], still touches zero; the fold
count is not a substitute for that interval.

**Consequences.** Gate 4 verdict line unchanged (ZSA ZGA-DEPENDENT by the rule). Reportable: in P1, splicing
progression is lower under A485 than in held-out controls (cross-fitted −0.27, 14/16 folds) and returns towards control
with Dux. Not reportable: the −0.44 magnitudes, and any magnitude for P2 or P3, which remain not cross-checked.

```
ae70df8563ea735be7f7fec8b22a69bd6daec708dea4049c835ac2b8ba6c5976  plan/POSTHOC_gate4_crossfit_frozen.md
87510c817a29c3a5ea2d47badd0827d5e54fc3eca677d4214ef939da109b3fbe  src/16_posthoc_gate4_crossfit.py
924e2b0e767149a690c1b7e60e30308a95ac2e0ebec85c46796938d4cbe57401  results/posthoc_gate4_crossfit_folds.tsv
ca9d863137cdf7c8aa357e214a0129a8a5756c1ab2f5373a0bbf52680d5743a8  results/posthoc_gate4_crossfit.md
```

## 2026-09-17 — Manuscript draft v1 (gz line), English and Korean; no new analysis

`gz/results/MANUSCRIPT_gz_v1_EN.md` (submission language) and `gz/results/MANUSCRIPT_gz_v1_KR.md` (same content, for
review). Target Aging Cell Research Article. English word counts: abstract 248; main text 5,435 (introduction 614,
results 1,923, discussion 840, methods 1,509, legends 549); 5 figures, 2 tables, 27 references.

- Every number is copied from the dated entries above; no analysis was run for the draft.
- References: all 27 retrieved from PubMed by PMID or DOI on 2026-09-17; entries with five or fewer authors list them in
  full as in PubMed. One title in `aging_cell_reference_literature.tsv` was a paraphrase (Chen & Zhang 2019); the
  PubMed title is used. GSE300734 has no publication and is cited only as a GEO accession.
- Figure 1 legend written to match the existing panels `gz/figures/FigGZ_1–3`; the GSE66582 series is assigned to a
  Supplementary Figure S1 that does not yet exist. Figures 2–5 exist only as review-page charts, not as publication files.
- Claim boundary kept: no causal verbs for observational results; post-hoc analyses labelled in section titles; the
  Gate 4 magnitude is not reported, only the cross-fitted direction; "rejuvenation" appears only in cited titles.
- Open points recorded in the draft: the title rests on a post-hoc result (an alternative title based on the prespecified
  results is given); embryo-related statements about Tyshkovskiy 2026 and Zakar-Polyák 2024 rest on the reference-table
  summaries and need checking against the papers; the hpi reference of GSE162345 needs checking against Tomikawa 2021.
- Stress-gene hypothesis (rescue arm) not tested, by decision on 2026-09-17: the central claim does not depend on it,
  a test would be weak (4 libraries per arm, non-physiological overexpression) and would add a further post-hoc layer.

Korean review page (artifact, private) built from the Korean draft with section notes and review charts.

```
64a53e60a2d9dd9703172fa7e2efafbf28ccf58572aee5866aa90e2a410f6339  results/MANUSCRIPT_gz_v1_EN.md
62dc3f542b5a01e7bde62eca2a98881ad6d117ac4e57ce084025cdb83c690f6b  results/MANUSCRIPT_gz_v1_KR.md
```

### 2026-09-17 — Figure-style survey, table limit, figure plan v1 (no analysis)

Two research agents read the figure legends of 16 related papers (8 clock papers, 7 ZGA/splicing papers read; Sakamoto
2024 and the published Xiao 2025 were not accessible — Xiao was read in its bioRxiv version). Plan written:
`gz/results/FIGURE_PLAN_gz_v1_KR.md` (5 main figures; library-level points before effect sizes; rescue trajectories in the
form of Xiao 2025 Fig 6H; contribution bars and cross-study scatter in the form of Tyshkovskiy 2026; ΔPSI quadrant plots
not recommended for the main text because the Gate 4 selection bias would carry into them).

Aging Cell allows at most 6 figures and 2 tables (2020 author checklist). The splicing event-count table added earlier the
same day as Table 3 was therefore moved to the Supporting Information as Table S2 in both drafts; text and legend
references updated.

**Disclosure.** While searching for an open copy of Sakamoto 2024, one agent sent the user's e-mail address to the
Unpaywall API (api.unpaywall.org) as the contact parameter, without the user's consent. The agent reports that nothing
else was sent. Recorded here and reported to the user.

## 2026-09-17 — Publication Figures 1–5 (display only; no new inference)

Decisions by the user: five main figures; schematics drawn in matplotlib. Following the figure plan, the post-hoc zygotic
timing panel (old FigGZ_2) moved to Supplementary Figure S2.

- `gz/src/17_figdata.py` writes `gz/figures/source_data/`: library-level values re-extracted from stored results, each
  checked before writing (zygotic-set group means equal `posthoc_rescue_programme_scores.tsv` to 1e-9; 23 GSE280522
  libraries; 13 zygotic genes). New displays: heat-map z-scores (13 zygotic + 17 two-cell/DUX-target genes; gene names
  from Entrez symbols, else GTF gene_name, two novel genes keep their Ensembl ID), key-gene log2(CPM + 1), and a PCA of
  PSI on the Gate 4 filtered events (missing values = event mean; display only). In every primary dataset all late
  two-cell libraries lie above all early two-cell libraries on PC1 (P1 variance 20%, P2 37%, P3 21%).
- Noted: key-gene log2FC differs slightly between the clock's Entrez-level matrix (text, e.g. Klf9 4.48) and all-gene
  Ensembl CPM (figure axis, Klf9 4.29); the legend states the basis.
- `gz/src/18_figures_publication.py` draws `gz/figures/pub/Figure1–5.{pdf,png}`: 167 mm wide (axes placed in mm),
  FreeSans, 600 dpi PNG. Figure 1 bootstrap intervals recomputed with the original seeds. Colour meaning fixed and
  validated for colour-vision deficiency: red/blue expression up/down; dark/light grey contributions pushing the clock
  down/up; black control, amber ZGA-blocking arm, teal A485+DUX (amber has contrast 2.6 on white, so points carry
  outlines and direct labels).
- Manuscript v1 (EN, KR): figure references and legends rewritten for the new panels, Methods 4.12 "Figures" added; the
  splicing event-count table is Table S2. English main text 5,839 words including legends (limit 7,500).
- Korean review page republished with the publication figures in place of the review charts.

```
8cf77f68c7286154ec3625727c6e6b8025887027cdaee173637b7013e53081d7  src/17_figdata.py
a3c437160da7f7f80fade49503d8a6f611792eac6c7fe5410f67c02bca956ff2  src/18_figures_publication.py
cdf4a9a686d20ffb7528dbb391cdddbae32e0dbe50a0c0564554b8b9fc087483  figures/pub/Figure1.pdf
7636083245d569bafb0d9cf2655012cdf9878ed02d339c27ad86493eb3dc25af  figures/pub/Figure2.pdf
d65727ec49d9708ad36c3e831f6a06b9375bc8bcd8c6a8e173ce0fdb9041f1fc  figures/pub/Figure3.pdf
a26af2c30cbb18e66afabc4f1bec0aaba019a2f22a78bea1f5f8d2b1cf04121d  figures/pub/Figure4.pdf
35953a49d0a59c7d431e85b2cf7bdb803685578c37bd8bfc6a12c3379c4ad849  figures/pub/Figure5.pdf
c796ad1973ec28d3bba5d1b3a9b96afee9d356dcc0bb2456b4c743ee2f859f80  results/MANUSCRIPT_gz_v1_EN.md
bf2dbd9f487e9fd7dc5e01b4afc087d3154e29b7f41411ba9d360223115c0cf2  results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-18 — Figures redrawn without line plots; two Oct4/Sox2 citations added (no analysis)

**Figures.** At the user's request, connected-line displays were replaced: Figure 1a and 2a are now stage grids with
schematic embryos (tiles for ZGA timing, maternal RNA shading, sampled stages, perturbation presence, collection dots);
Figure 1b–c show stage means as diamonds with 95% bootstrap bars and no connecting lines; Figures 2b, 3d and 4b show every
library as open (E2C) or filled (L2C) circles with mean bars, grouped by arm (D printed above each pair in 2b); Figure 5d
shows the per-fold cross-fitted differences I and R as points (16 each) instead of fold-joining lines. Stage marker
unified across figures (open = early, filled = late two-cell). Forest plots kept. No value changed; legends in both drafts
rewritten to match. Figure 2 height 191 mm, all widths 167 mm.

**Citations.** Chandramohan et al. 2026 (Dev Cell 61:621–637.e5, PMID 41338198) — full text not accessible (Cell Press and
ScienceDirect returned 403; not in Europe PMC); GEO GSE275652 (E4.5 scRNA-seq, 1 control and 2 maternal-zygotic Sox2 KO
samples) read for design. Hou et al. 2025 (eLife 13:RP100735, PMID 40014376; GEO GSE264615/GSE264614), abstract read in
PubMed. Both cited in Discussion for abstract-level statements only: Oct4/Sox2 after the 8-cell stage gate epiblast
developmental capacity and epithelialisation timing; Oct4/Sox2 take part in the morula-to-early-ICM transcriptome
reorganisation. Added as an open question whether the second clock decrease (16-cell → early blastocyst, GSE45719,
−0.241, Gate 2b R1 entry) depends on the Oct4–Sox2 wave. References now 29. English main text 6,023 words incl. legends.

```
beb2de97fa984fa554798d8c8d83b03822fa023e98a4fd8add726b2c14721eee  src/18_figures_publication.py
b4fcbcf4d76ba7ec7f417878b36979468e3c424e88d5ba80c705ba71c73692f5  figures/pub/Figure1.pdf
52ad2609a391e9ad00ce655ce836e9ddede71ccfa273e794cf8cbd7fde548506  figures/pub/Figure2.pdf
8ae0d46b6d1f7a3a3aaa7ecc7b4ec90099afa7f0d2c2703d18c42e68471d10b7  figures/pub/Figure3.pdf
a6f44072eee221ef742aa3d18c0df76c517a2774e9040bfd56c4b457d67b1d05  figures/pub/Figure4.pdf
3bb1df316219be478f68ca435c988b1c1dde7af8deaa8631c0277f0feeb7ba05  figures/pub/Figure5.pdf
08aa9709b59f5069f88c588994ee75f1f6e2798c7f2cc499bf1f03e4c216e62b  results/MANUSCRIPT_gz_v1_EN.md
196bfaac0a5115d222d27e201aabb5ab643f47ceebf38bfe5ae35e488caf7ae1  results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-18 — Four panels moved to Supporting Information; Supplementary Figures S1–S7 drawn (no analysis)

**Panel moves (user decision via question prompt).** Figure 2d (secondary single-stage forest) → Figure S3; Figure 3d
(expression of five key contributors) → Figure S4; Figure 4d (per-gene contribution scatter) → Figure S5; Figure 5a
(PSI principal-component display) → Figure S6. Remaining panels re-laid out: Figure 2c widened; Figure 3b full height,
figure height 150 → 138 mm; Figure 4b, 4c widened; Figure 5 relettered (a per-library progression, b per-fold
cross-fitted I and R, fold summary moved into the row labels, c forest with cross-fitted rows).

**New supplementary figures.** S1 GSE66582 library clock values (V0/V2; points and mean bars; stored table
`gate2b_R2_library_tage.tsv`; E2C → 2C change −0.260 / −0.182 recomputed and equal to Gate 2b R2). S2 post-hoc zygotic
timing redrawn from `gate1alt_intervals.tsv` in the publication style (diamonds with 95% bootstrap bars instead of bars;
pooled Spearman recomputed 0.524, equal to the Gate 1-alt entry). S7 library QC: `src/19_figdata_supp.py` re-applies the
Gate 3 rule (p_pseudoaligned ≥ 30; log10 detected genes ≥ dataset median − 3 MAD) to all 76 libraries and asserts that
the passing set equals `results/gate3_library_tage.tsv` per dataset (passed; excluded SRR30827777, SRR34172919,
SRR26923904, as before). FigGZ_2 is superseded by Figure S2 for the manuscript; the file is kept.

**Manuscripts.** EN and KR: Results references updated (2.3 → S3; 2.4 → 3b + S4; 2.5 → S5; 2.6 PCA → S6, progression → 5a,
cross-fit → 5b, c; QC exclusion → S7); Figure 2–5 legends rewritten to the remaining panels; Methods 4.12 PCA reference →
S6 and a sentence on why library points are not joined; new "Supporting Information figure legends" section (S1–S7).
English main text 6,016 words incl. main legends (supplementary legends excluded). Review page republished (version 4).
All figures 167 mm wide (heights: F1 176, F2 190, F3 138, F4 184, F5 146, S1 62, S2 118, S3 80, S4 76, S5 62, S6 56,
S7 96 mm). PDF checksums change at every render (embedded creation date); the PNG content is deterministic.

```
7db3b2ae3acacb57d2730d392f8d0351efa5f2be06b635c019e0908951f46a0e  src/18_figures_publication.py
61f7b0d88cd1cbf70cf0beff947e23eb1219790e0c6b173b24a22a6b062c2659  src/19_figdata_supp.py
c14748bfee3c5b5f1929e04f2b4c266f48051b3034aa2325b7e86857f0cf7e49  figures/source_data/figS7_library_qc.tsv
767958e928fbad082885bca4846af6e3dc0d354a81eac1097d57b55cc4e6308f  figures/pub/Figure1.pdf
bc209a4781f32f4ae7d78b07f66b27b91278421f3e02adc24d8500a32ced9feb  figures/pub/Figure2.pdf
82b46173a7c247214b3c58cdc0f6838e594d02902cf183c3ec121128d6b931ab  figures/pub/Figure3.pdf
208b8b2145301812905f7d404c56fae00bf805c2ba928e43c3473b61e5684f73  figures/pub/Figure4.pdf
1a3596b3163c473612340787518c4f19f4de9d48c4bd6016d956eab090232f34  figures/pub/Figure5.pdf
755c4769509add8697410edf9997a8f721c7a790e6739dcdf7c43ff15367e576  figures/supp/FigureS1.pdf
3bedb1a99a6c2e3a1354cfc6992b9841e2c2af62f4499a9cbe97433529e0ae76  figures/supp/FigureS2.pdf
35583c0e3a9346cbb8dd7ad16e9bab161060eafa13cb2e42447d093868c083fa  figures/supp/FigureS3.pdf
c0b030e260839c71feacf250d97bbe8e3addcf7e1eae76e0b1f148acf2406a36  figures/supp/FigureS4.pdf
fbfabb5233fd292f194584806f2bf8221528d6f92a45f05c52e1d16313129921  figures/supp/FigureS5.pdf
41fc62cf5f031c00627688b1d7dd9ff6514bc9356cc514e1deabe3144880ddfd  figures/supp/FigureS6.pdf
b9313e1c1af7900411af209ad3c282c0de3a0a7472ae92520cd63cc2392439da  figures/supp/FigureS7.pdf
5a1d3aad153363fb23728f2be3792de8fb58b88cb05444f1043cef11c7071003  results/MANUSCRIPT_gz_v1_EN.md
8efcd5e5f5f46f25439fcedbae8afd945a01ec7cb7eefdaf87434ca72a9bef76  results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-18 — Absolute paths removed from scripts; public code repository prepared (no analysis)

**Why.** The user created the public GitHub repository `vincentjshim-rgb/2C_ZGA` for the code. Scripts carried the
machine-specific prefix `/home/shim/.../analysis`, so they could not run from a clone.

**Change.** In the 20 Python scripts of `gz/src/` that contained the prefix, two lines were inserted after the module
docstring (`import os`; `ANALYSIS = dirname(dirname(dirname(abspath(__file__))))`) and each literal prefix was replaced by
`f'{ANALYSIS}'` (38 occurrences). A reverse substitution reproduced every original file exactly (checked in the edit
script). Shell scripts 07, 10 and 10u now derive `B` from their own location, take the index from `$IDX` (default
`analysis/data/processed/mm_full.idx`) and kallisto from `$KALLISTO` (default `kallisto` on PATH); `bash -n` passed.
Originals kept in the session scratchpad (`src_before_relpath/`). Frozen plans and this log were not edited.

**Check.** Run from `/`: `00_validate_tage_py.py` OVERALL PASS; `19_figdata_supp.py` assertions passed;
`18_figures_publication.py 3` re-rendered. `results/00_tage_py_validation.tsv`, `figures/source_data/figS7_library_qc.tsv`
and `figures/pub/Figure3.png` are byte-identical to the pre-change files.

**Repository.** Git initialised at `analysis/` with a whitelist `.gitignore`: tracked = `gz/src`, `gz/plan`, `gz/results`,
`gz/figures`, `gz/metadata`, `gz/requirements_gz.txt`, `src/figstyle.py`, `src/50_tx2gene_full.py`,
`metadata/tx2gene_full.tsv`, this log and the Gate 3–4/post-hoc run logs (205 files, 43 MB). Not tracked: all data
(~196 GB), `gz/metadata/gene_orthologs.gz` (124 MB, over the GitHub limit), manuscript drafts, summaries and the
figure plan, third-party SUPPA2 and tAge code and the tAge model files. `gz/src/tage_py.py` is held back: the tAge MGB
Open Access License (section 4) requires modified materials to be shared back with the licensor before distribution;
decision left to the user. Line endings kept as deposited (`core.autocrlf false`) so recorded checksums stay valid.
Commit identity: GitHub noreply address (the personal e-mail is not written into commits). README.md added.

```
2d9adc713c70b84cda2c36c4253537e5774884e67292d5dbd431ec34eff30ff8  gz/src/00_validate_tage_py.py
a3e67a1131dbc95db34bc260ee129cfb51337100d21ea5b7d4a64c4f38377db4  gz/src/01_gate1_cross_species.py
dac2e4c49de95fa478f6ce84c1b5e12954fe6e72b541e2abbd75d812b5f28c81  gz/src/02_fig_gate1_trajectories.py
6d9ea309f8294d67cf49f8dd4ea5bfc6f6de2e4e8500004e17b9dfe8948454f0  gz/src/03_gate1alt_zga_timing.py
9fbe8cb75818242325d0868b8ef0d790fb0462b4e21b50f8bbd1110bdaec2bab  gz/src/04_fig_gate1alt_zga_vs_clock.py
e0f5e611868b8f1e15989a571c03c67cc3e2435639e76badee0afcf8c2105bd1  gz/src/05_gate2a_composition.py
c87af651872bd7f881095edb5ea01eb00eea79f03f355e0b6b7a9f138441b0dd  gz/src/06_gate2b_R1_deng.py
3b915c405a23588c0c592313ea73c40a404c714a07b8792b7c4218d1a7c6609b  gz/src/08_gate2b_R2_wu.py
47e1bf79915b8d01ab7b667dabf4ed4ed711afc30b0de489c77c7a7d3c4c72fd  gz/src/09_fig_gate2.py
7c1a7d715a7b5ee0a04254c6afd2f80a2aaf95dbfcf727be54de788911640adc  gz/src/10u_make_worklists.py
983dd4d2d4f3b766a88c8ed6b343c48e849f7f34fd2d7a1d1668010f8ef70d30  gz/src/11_gate3_analysis.py
7f11fb8478c308840fe28b36db4dd39ca82f5dc6726d49a6c114f1f5995188eb  gz/src/12u_quant_concordance.py
9c08b0930ab341c733a302a42ecafcd8653e8b7051c09ef32ca76a88d40b27c7  gz/src/13_posthoc_gate3_contribution.py
2a0a765947079ea8598a99e86428369a87912d81c383379f3e301c0939cad4c2  gz/src/14_posthoc_rescue.py
92309592cbff415be53197da5a9a624e6d2b6d8230657bada6fb8bc8cab8df79  gz/src/15_gate4_analysis.py
b18de89c8c194dac0c1843fbff2dfc7fc0add98f2f7c4f07248247b392ad14fc  gz/src/16_posthoc_gate4_crossfit.py
274fdc1cc71a2eebcb7b265e9b2486a7ef4e2759bf6f353b6eacf861135a2d5a  gz/src/17_figdata.py
a9cba70c76e29175fe426986d76fd18eb5ea606d1647a2250d137bebd9d04307  gz/src/18_figures_publication.py
096cb99271c07b3523256679cebad8c52fc0ce0acc659e35ccc6d5bc03884fec  gz/src/19_figdata_supp.py
cf5dc612af90885599526d454600d381421616fa151c63356b64353dd75296fa  gz/src/07_gate2b_R2_quant.sh
93abe2b312e8a4258074b2423339f2ac5caf035b3ee2e65c1c67f5f295554699  gz/src/10_gate3_download_quant.sh
25ebe0899f05140d34ec925dc426623207d4922e20484ac6c5f49915f0d3b5b2  gz/src/10u_gate3_download_quant.sh
6671b9fc041979ffb1bb07fb3cc671da98af832deb15779fa60e88507c1fb0c2  README.md
f4d6599148fd28b3ae8f38686dc14bd71b71e5e750af421eed66c111230060e2  .gitignore
```

## 2026-09-18 — Supplementary Tables S3–S5 built from stored results (no analysis)

`gz/src/20_supp_tables.py` reshapes stored result files into three Supporting Information tables and checks each against
the stored summaries before writing: S3, the clock value of every scored embryo or library (557 rows: GSE225056 334 incl.
18 rhesus scored for description only, GSE45719 39 pseudobulks × V0/V2, GSE66582 15 × V0/V2, Gate 3 115 rows), with GSM
and SRR where they exist; S4, per-gene contributions for the 1,839 clock genes (P1 control/A485/A485+DUX, P3
control/Brg1 matKO, coefficient, expression change, dynamic-gene class); S5, the 16 cross-fitted folds. Checks passed:
the P1 control drop recomputed from S3 equals `gate3_arm_drops.tsv` (< 1e-12); S4 contributions sum to the same drop
(< 1e-9) and agree with `posthoc_rescue_contributions.tsv`; S5 reproduces 14 of 16 negative interactions and 16 of 16
positive rescue differences. Output `gz/results/supp_tables/` (three TSV plus one xlsx with a README sheet).
Both manuscripts now cite S3 in Section 2.3, S4 in 2.4 and S5 in 2.6, carry the three table legends in the Supporting
Information section, and name them in Data availability.

```
3e855ae98d44aacec2729a9de8764ed8f1dc2bc6b899c4c0b50e8dc846efd177  gz/src/20_supp_tables.py
ee72e5aeb306fe3aff124edfca047ed89730dea962edb36429405c6ddbbb84d3  gz/results/supp_tables/SupplementaryTables_S3-S5.xlsx
d9e9714014ea3447e049a71ecd0b60a72d76c7c3f4488a24146b27c0bf093d9d  gz/results/supp_tables/TableS3_clock_values.tsv
576f962990dd0a7fde98d3358e6e07b4a04cf21d39d660472955fb3c566da32a  gz/results/supp_tables/TableS4_gene_contributions.tsv
8d04fa628c2a94ca66b4105f9c0d9a33ddd5d9acb2e37fdb46c6592e9c5d260e  gz/results/supp_tables/TableS5_crossfit_folds.tsv
3e630990c33a46a98ed513af92ac36024e3d2d7ec25d970d97cd4d48101961f1  gz/results/MANUSCRIPT_gz_v1_EN.md
54c8a5a0db205d4ec64190fe4d222e130e238069b967e3c9d825c51dcc148392  gz/results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-18 — Source-text verification, literature check, title change (no new analysis)

**Source-text verification (agent-assisted, every item re-checked here against the source or the local data).** Corrected
in both drafts: (i) GSE225056 has no blastocysts — stages end at the 16-cell (mouse) or morula (cow, pig, rabbit),
verified in the deposited matrices; (ii) Table 1's second library number was the number downloaded for the contrast, not
the number deposited — series totals added (GSE221985 40, GSE300734 25, GSE235547 88, GSE248499 25, GSE162345 30,
GSE66582 32, GSE45719 317, GSE225056 361); (iii) GSE162345 rewritten from Tomikawa 2021 STAR Methods — C2C12 nuclei fused
into intact, non-enucleated demecolcine-arrested two-cell embryos at 21–24 or 30–33 hpi, collected 45/54 hpi, hpi = hours
post insemination, and the libraries are host-plus-donor mixtures (now stated); (iv) Methods 4.4 no longer says the
features passed through YuGene (the scaleddiff models use the scaled branch; checked in `src/tage_py.py`), states that
values are in normalised-age units (fraction of species maximum lifespan; `clocks_metadata.csv` lifespan_scaled = TRUE with
no species factor applied) and that cow/pig/rabbit lie outside the clock's training species; (v) the splicing claim
"extends them to pharmacological and maternal perturbations" was wrong — Zhang 2024 already used α-amanitin and maternal
Btg4/Pabpn1l loss (verified in PMC11005748); the sentence now claims only p300/CBP inhibition and maternal Brg1 loss.
Also: GSE248499 → Matoba 2024 is correct (GEO PubMed 38729154; the paper's data-availability names the accession), and the
SCNT arm used is the deposit's "SCNT, cont"; GSE300734 now has a matching publication (Shi et al. 2026, Fundamental
Research, doi 10.1016/j.fmre.2026.03.007, Crossref-verified), cited with an explicit note that its text was not accessible
and the accession could not be confirmed inside it; GSE235547 libraries used are the total RNA-seq ones.

**Literature check (two agents, verdicts verified here).** The per-gene decomposition is not new: Tyshkovskiy 2026
contains "Gene contribution analysis" (clock coefficient × effect size) and "Module contribution analysis"
(tAge = C0 + Σ_module Σ_g w_g·expr_g), correlates contribution vectors across models, and reports anti-correlated gene
contributions before and after the E10 minimum in mouse embryogenesis — read directly in PMC13233323. The drafts now cite
those methods and claim only what is added: the signed budget within one comparison, the perturbed and rescued arms, and
the preimplantation window (that paper's embryo series GSE39897 runs egg → newborn with no preimplantation stages; the
words "2-cell", "blastocyst", "preimplantation" do not occur in it). No published work was found that applies an ageing
clock to ZGA-blocked embryos, that dissociates a scalar clock from its restored components, or that cross-fits a
reference-defined developmental score.

**Title.** Changed on the user's decision to "Gene-level decomposition links the two-cell transcriptomic-age decrease in
mouse embryos to minor zygotic genome activation"; the former title became the alternative title.

**References.** 30 → 42. Added (all resolved through Crossref or Europe PMC on 2026-09-18): Ambroise & McLachlan 2002,
Chernozhukov 2018, Deng 2026 NAR, Higgins-Chen 2022, Isaev & Knowles 2025 (bioRxiv), Kriukov 2024, Li 2024 Cell,
Sehgal 2025, Tomusiak 2024, Weston 2019, Xing 2020, Zhang/Tyshkovskiy 2026 (bioRxiv). Every reference is cited in the text
and every citation is listed (checked by script). English main text incl. main figure legends 6,902 of 7,500 words.

```
78bd328e7baba39b7a5edddda4059f4f6ee4c4e3578828eac9b9b2468694c26f  gz/results/MANUSCRIPT_gz_v1_EN.md
1a06ad05959ea52b16fc614f12e4a89af13fbce8f7e799c925e74d9e9df57e65  gz/results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-18 — Methods shortened; cross-species precedence checked (no new analysis)

**Methods.** Section 4 reduced from 1,863 to 1,585 words in both drafts by removing process detail that is in this log and
the repository (ENA/SRA mirroring and SRA-Lite note, cross-machine index rebuild, tAge commit hash and v1.0.0 bug note,
mortality-clock validation value, the YuGene branch sentence, unused arms of GSE248499/GSE235547) and by stating the
reference libraries of every dataset once in Section 4.4. No reported method or number changed.

**Cross-species precedence (agent-assisted, key claims re-verified here).** Europe PMC title/abstract search for rabbit
or Oryctolagus with clock/biological-age terms returns 0 records (re-run here); the whole clock-and-embryo field is 25
records. No published work applies an ageing clock to the preimplantation embryos of more than two species: Kerepesi 2021
(mouse), Kerepesi & Gladyshev 2023 (human, and reports no preimplantation change), Tyshkovskiy 2026 (mouse only, no
preimplantation stages). No study tests coincidence between a clock decrease and the ZGA interval. Three defensive
citations added after verification: Schaetzlein & Rudolph 2005 (mouse and cow telomere lengthening at the
morula-to-blastocyst transition, PMID 15745634), Hao et al. 2026 (bioRxiv 10.64898/2026.08.25.746714, abstract read:
morula as the nadir of age-associated methylation entropy in human), Li et al. 2025 (Biology of Reproduction 113:541-556,
five-species SCNT/IVF ZGA comparison without an age axis). Introduction and Discussion now state that prior clock
analyses stayed within one species, that none tested the ZGA interval, and that rabbit had no molecular-age readout;
the cross-species claim is framed as breadth plus a negative result. References 42 → 45. English main text 6,772 words.

```
9e3dd7c613d926f102a911503963aefe3291449e32d3fe00c527728bf03eb21b  gz/results/MANUSCRIPT_gz_v1_EN.md
02937a99a265dd71a682495629b34cb43c844bf666fecfba58e7af2c6187c6cb  gz/results/MANUSCRIPT_gz_v1_KR.md
```
