
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

## 2026-09-18 — Post hoc: global remodelling vs gene-specific change (plan frozen before the numbers)

**Plan** `gz/plan/POSTHOC_global_remodelling_frozen.md`, written and frozen before any value was computed, with the
reading rule and the decision about what may enter the paper. Script `gz/src/21_posthoc_global_remodelling.py`
(seed 20260920); the script asserts that every recomputed arm drop equals `gate3_arm_drops.tsv` (< 1e-9).

**Results (V0).** Split D = G + S with G = (Σ coefficients) × (mean feature change):
P1 GSE280522 control D = −0.2375, G = +0.0072, S = −0.2447 (|S|/|G| = 34); P3 GSE300734 control D = −0.1714,
G = +0.0122, S = −0.1836 (|S|/|G| = 15). The global component is therefore ~3–7% of the drop and of the opposite sign.
Coefficient-reassignment null (2,000 draws, observed feature changes held fixed): null SD 0.125 (P1) and 0.133 (P3);
the observed drop sits at the 2.6th percentile (z = −1.93) in P1 and the 8.1st percentile (z = −1.40) in P3. The A485
interaction +0.1809 sits at the 96.8th percentile (z = +1.88); the Brg1 interaction +0.1091 at the 91.5th (z = +1.34).
Direction-restricted partial drops: in the P1 control arm the drop comes almost entirely from genes whose feature rose
(−0.217 of −0.238).

**Verdict and a gap in the frozen rule.** The rule required, for "NOT EXPLAINED", both |S| > |G| and the observed drop
below the 2.5th percentile of the null in both datasets. |S| > |G| holds by a wide margin in both, but the percentile
condition fails (2.6% in P1, 8.1% in P3), and the plan's third label ("EXPLAINED" if |G| ≥ |S| in both) does not fit
either. The rule as written therefore does not classify this outcome; recorded here rather than reinterpreted. The
substantive reading is two-sided: the decrease is not a uniform global shift of the features, and it is also not
exceptional relative to a random assignment of the same coefficients to genes given the spread of feature changes in
this window.

**Effect on the paper.** Following the pre-decision for anything other than "NOT EXPLAINED": no new figure, table or
Results paragraph. One sentence was added to the Discussion (both halves of the reading) and one sentence to Methods 4.8
describing the control. English main text 6,772 → 6,904 words.

```
e9372d167ea0065b21559c0113ae9619570a2f83f9114028fc6f1053cb8c00db  gz/plan/POSTHOC_global_remodelling_frozen.md
0568ccc9874cecb6527ccd978562d634f768fbdafde4e43d871b90a4c1ac113c  gz/src/21_posthoc_global_remodelling.py
69a5cfd9e0dafbd36a8b86bb8f664f0992218c4fc1a292ccf603eb227693da40  gz/results/posthoc_global_remodelling.tsv
35481fc52b2c54d69ac49708b407428efff889f037c94075ed006fe1705bcfe2  gz/results/posthoc_global_remodelling.md
e8f42506f0ff4efeead1396cc101920f22323efc7f23a385df0255d4bda4a3ee  gz/results/MANUSCRIPT_gz_v1_EN.md
4651bc14a923026670944f3cdf7268b4f09e31951cb2e938287405ed2a5afade  gz/results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-18 — Submission package assembled (no analysis)

Aging Cell's own format was checked empirically rather than from the guidelines page (Wiley returns 403 to automated
requests): the reference style of a published Aging Cell article (PMC11561706) is APA 7th with full author lists,
19 authors + ellipsis + last author beyond 20, volume(issue), pages and DOI, and its section order is Abstract,
Introduction, Results, Discussion, Methods, Author contributions, Funding information, Conflict of interest statement,
Acknowledgements, Data availability statement, References.

`scratchpad/build_submission.py` writes `gz/results/MANUSCRIPT_gz_v1_EN_submission.md` from the working draft: title
page with the required elements (author fields left as placeholders), references regenerated in APA 7th from Crossref
metadata for all 45 entries (full author lists, issue numbers where Crossref has them; 39 of 45), the data and code
availability text moved out of Methods into its own statement, and the statement sections added as templates. Main text
6,847 words including figure legends; abstract 248 words.

The availability statement now names the public repository and states that the Python port of the clock preprocessing is
not redistributed there (source-package licence) but is available on request; the same sentence was added to the Korean
draft. Package in `gz/results/submission/`: manuscript, cover letter and supporting information as .docx (pandoc),
Figures 1-5 and S1-S7 as vector PDF, Table S1 as TSV, Tables S3-S5 as the xlsx workbook, plus a README listing what the
author still has to complete. `.gitignore` extended so the submission package, the supporting-information markdown, the
cover letter and the Korean checklist stay out of the public repository; only the code, plans, results, figures and this
log are published there.

```
4e107fcb5ff16a5f82d6d4825b0a786c0dd3cd0b3bc500784b3497be0e2fe91a  gz/results/MANUSCRIPT_gz_v1_EN_submission.md
a2f3df1fadd7c5d214877ba6f34e9bfce65b571481156edc25d12d1c20a4e1aa  gz/results/submission/Manuscript_AgingCell.docx
d2bf9d9da8c918eaf97879a9d7c671cc1915b639f1306a6da9c25035baecd7a5  gz/results/submission/Supporting_Information.docx
32a2c077619e5ddc3b7efa5b87a022a8db16d98c980bf0cc369bd1f97d5bbfe5  gz/results/COVER_LETTER_draft.md
```

## 2026-09-19 — Figures 2–5 and S3 redrawn in standard forms (display only; no value changed)

**Why.** The author judged the forest plots (Figures 2c, 5c) and the "net" flow bars (3a, 4c) unclear and atypical, and
asked for figure types standard in the splicing and clock literature, keeping only panels the logic chain needs.

**What changed.** Figure 2: forest panel removed; D and the interaction I with its 95% bootstrap interval are printed on
the per-library panel (2b), R in the footnote. Figure 3: (a) is now a ranked-contribution strip over a running-sum curve
(GSEA-style) for the P1 and P3 control arms — the sum falls to −1.40 (P3 −1.48) and returns to −0.24 (−0.17); (b) paired
bars, control vs A485, for the 20 largest downward contributors, with log2FC printed; (c) unchanged; the axis label
"β × Δ feature" replaced by "contribution to the clock change" with the definition in the legend. Figure 4: (c) is now
the running sum of the three arms in the control gene order, so the curves compare the same genes: at the 666th gene
−1.40 / −0.76 / −1.19 (the 85% recovery of the text), end points −0.24 / −0.06 / −0.06. Figure 5: (a) PSI heat map of the
511 control-defined events × 23 P1 libraries (z-scored per event); (b) cumulative distribution of the per-event oriented
ΔPSI per arm (medians 0.32 / 0.16 / 0.25; the legend states that the control distribution is shifted by construction);
(c) progression strip and (d) cross-fitted folds unchanged; forest panel removed. Supplementary S3: per-library points
instead of a forest. `src/22_figdata_v2.py` writes the new source data and checks each against the stored results
(running sums end at the Gate 3 drops; progression recomputed from the PSI matrix equals Gate 4; 85% recovery
reproduced). Heights (mm): F2 124, F3 136, F4 184, F5 152; all 167 mm wide.

**Manuscripts.** Legends of Figures 2–5 and S3 rewritten in both drafts; Results references updated (2b; 3a once; 4c on
the 85% sentence; 5a,b new sentence on the PSI heat map and ΔPSI distribution in 2.6; progression → 5c; cross-fit → 5d);
Methods 4.12 describes the cumulative curves. English main text incl. legends 7,088 words (legends grew by ~190 words).
Submission package (docx, PDFs, SI) rebuilt; review page rebuilt.

```
545bfa9807a9792d16eef2bbe17ff30dc46f407c7940ed5fda39deb7c4e18579  gz/src/18_figures_publication.py
ad84f86d16446ec3a53763ae0fd3361fa14f6b0a56f2fa7772eabac4183ded5d  gz/src/22_figdata_v2.py
f80c183d80f573bb46c1b557afc8b5ad9e4c9919d0662f665e60acea6eb53a41  gz/figures/source_data/fig3a_4c_cumulative_contributions.tsv
4b2bafc056790cdfd55462c9b61490f7ac1d66926fc7362073434bfe020c0fd3  gz/figures/source_data/fig5a_psi_zsa_events_P1.tsv
80db09477dd8fd0f9cf6b901a31e89f9f6d1e6eccea7c4b2c148a2a7e00a1bb9  gz/figures/source_data/fig5b_dpsi_per_arm_P1.tsv
b1ec41ba08134deb1a2487e3de03707dcda80b4ab637508e1164ca4a61f56672  gz/figures/pub/Figure2.png
4411843ff260fb9e2d141a6e2be59f1b1c06e914490f7ef79c4c15cd945f4df5  gz/figures/pub/Figure3.png
c29c766787cc06d35e4a1d81ecd779f80626e96b56a38957b4fb346c49e5d786  gz/figures/pub/Figure4.png
bab5ef77f4f3d49248d76505369f824be71d48616ac04eb27fb815d0ea112dd7  gz/figures/pub/Figure5.png
476d9e4759eabb8c1a023b57829f99fed76795c52de9a1bbfbd1bd35a8bed782  gz/figures/supp/FigureS3.png
657f7fc68be39d75490418ebf0c350c2254ee7f5a5291e0a97d0ffd0876eab18  gz/results/MANUSCRIPT_gz_v1_EN.md
6cfa8f2fe65e936ca281656e56d9f77b3106a6020592b90bdd8510f413d44715  gz/results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-19 — Figures 3 and 4 merged; verdict wording; clearance-route follow-up; exec bug fixed

**Figures.** The decomposition figure (a–c) and the rescue figure (d–f) were merged into one Figure 3 (167 × 224 mm) at
the author's request; the splicing figure is now Figure 4. No panel content changed. Both drafts renumbered (old 4a/b/c →
3d/e/f; old 5a–d → 4a–d), legends merged, title-page count 4 figures. Submission package and local review page rebuilt.

**Verdict wording.** The word "mixed" is kept where it is the frozen plan's own label (Table 2 and the rule statement in
Methods 4.7) and replaced in prose by what happened: the abstract now says the smallest dataset failed the prespecified
rule; 2.1 states the criterion was met in one of four species; 2.3 states that all three V0 interactions were positive
and that the smallest dataset reversed sign in V2, so the rule requiring agreement between the two gene sets was not met
overall; the Discussion opens with the two rules that were not met in full. Nothing was removed: every prespecified
verdict remains in Table 2 and in the text, and the plans and this log are public.

**Clearance route (post hoc; plan `gz/plan/POSTHOC_clearance_route_frozen.md` written before the number).** Mouse
GSE225056 composition-only simulation decomposed gene by gene (`gz/src/23_posthoc_clearance_route.py`; Δ_sim = −0.1337
equals `gate2a_simulation.tsv` to 1e-9). Of the 1,839 clock genes, 235 belong to the Gate 2a maternal set (defined from
the oocyte across all later stages), whereas the within-window V2 definition in GSE280522 contains eight maternal clock
genes — the two "maternal" sets differ in scope, which is why Sections 2.2 and 2.4 are not in conflict. 55% of the
simulated change ran through the 235 maternal clock genes and 45% through the features of the other genes (shifted by
the normalisation when abundant transcripts are removed). Output `gz/results/posthoc_clearance_route.tsv`. Not yet in
the manuscript; candidate sentence for 2.4.

**Second decrease (for the story-check gap 5).** Already in stored results: GSE45719 16-cell → early blastocyst −0.241
(95% CI −0.296 to −0.188; `gate2b_R1_intervals.tsv`); GSE66582 8-cell → ICM −0.245 (n = 2 and 3, direction only;
`gate2b_R2_library_tage.tsv`). Not yet in the manuscript.

**Bug fixed.** After the 2026-09-18 path change, `src/01_gate1_cross_species.py` derives ANALYSIS from `__file__`; the two
scripts that `exec` its header (`05_gate2a_composition.py`, `03_gate1alt_zga_timing.py`) passed an empty namespace and
would have failed on re-run (found when script 23 reused the same pattern). Both now pass `__file__`. No stored result is
affected (both were last run before the path change, and script 23 reproduces the Gate 2a simulation value exactly).

```
f4009d4e35b64eebb644a15f115eb8b6842f2a45c9be05f906a7f809fea7c9cd  gz/src/18_figures_publication.py
cbf0e21209956822ddb7ae2f48512e7bd42bd2d151e744059226079c694543f6  gz/src/23_posthoc_clearance_route.py
05641c83684235311b9a33ce7f0e1e61ae933455c37b0b98658a2af23c5052c8  gz/src/05_gate2a_composition.py
963ea5d118c9f7b886e6b8f801560a81fe7cfc3d2379f19b29598bd6e3a1aa90  gz/src/03_gate1alt_zga_timing.py
1218b9409bf8ec9b1b4d2b224712d409d50d8e42ad0d7c803730c671beef960d  gz/plan/POSTHOC_clearance_route_frozen.md
110d13c5c9e7825ee39ca7c8c26347bd02e4854b12f1e407e32212a2544d2941  gz/results/posthoc_clearance_route.tsv
0ea2684e4fa3f534f5518e919e97ad36981dcc3a496d399c11247ff58bf05a2e  gz/figures/pub/Figure3.png
bab5ef77f4f3d49248d76505369f824be71d48616ac04eb27fb815d0ea112dd7  gz/figures/pub/Figure4.png
98ba75eba5b9bd62f3190728ea873f7deed0674d1c9dd959006f35cffaf09ac2  gz/results/MANUSCRIPT_gz_v1_EN.md
8b83e6547d2b5a7dcab370a3ec06b45363356954d240a11e2dd36458b912bf16  gz/results/MANUSCRIPT_gz_v1_KR.md
```

## 2026-09-19 — Six story-check gaps closed in both drafts (text only)

(1) Discussion: splicing and clock readouts agree in direction but not magnitude across perturbations (ρ = 0.00; Tardbp
+0.04) — one sentence before the recommendations. (2) Results 2.4: why 2.2 and 2.4 do not conflict — the Gate 2a maternal
set (oocyte-defined) holds 235 clock genes, the within-window V2 set eight; 55%/45% split of the simulated clearance
effect from `posthoc_clearance_route.tsv` (post hoc, plan frozen first); matching sentence in Methods 4.8. (3) Discussion:
the splicing readout met its rule where the clock readout did not. (4) Discussion: the cross-species post-hoc correlation
(ρ = 0.52) beside the A485 sentence. (5) Results 2.1: second decrease, GSE45719 16-cell → early blastocyst −0.24
(95% CI −0.30 to −0.19) and GSE66582 ICM −0.25 below the 8-cell libraries (n = 3 and 2, direction only). (6) Abstract:
closing sentence now names the ZGA link and the scalar's blind spot; trimmed elsewhere to stay at 250 words. English main
text incl. legends now as printed by the build. Submission package and local review page rebuilt.

```
4bfb12619696c8f45cb390d3a3edb2d417816ba1be86f5e07aa26c85045acdb4  gz/results/MANUSCRIPT_gz_v1_EN.md
002119eef4bce0ebb7609c104603f48b7da5a1e9a281bd031fa8955b692b4a37  gz/results/MANUSCRIPT_gz_v1_KR.md
737cebd109d83efaf7fc333452aebcd9156e2618055a9bf8dc3c4037f2303712  gz/results/MANUSCRIPT_gz_v1_EN_submission.md
```

## 2026-09-19 — Discussion: cited context for the largest contributors and for SCNT; six references added

Agent-assisted literature check, each citation re-resolved in Crossref here (51 references now). Added to the Discussion in
both drafts: (i) the five largest downward contributors are all zygotically activated classes in DBTMEE (Park et al.
2013 Genes Dev; Park et al. 2015 NAR), consistent with their loss under a minor-ZGA block given that minor ZGA is
required for major ZGA (Abe et al. 2018 PNAS); in the clock paper's own multi-tissue ageing signature (Supplementary
Table 2) Klf9, Neto2, Gpatch4 and Psmb5 decline with age and Pi4k2a rises — stated as a caution that regularised
coefficients are model properties, not marginal age trends; (ii) the source study's embryonic GSEA (cell-cycle, MYC-target
and mRNA-splicing programmes up, inflammatory and interferon programmes down) as the combination that lowers tAge in adult
tissues; (iii) SCNT context: reprogramming-resistant regions and donor-transcript retention (Matoba et al. 2014 Cell),
Kdm4d narrowing but not closing the gap (Matoba et al. 2024), and the split evidence on clone age (Lanza et al. 2000
Science; Ogonuki et al. 2002 Nat Genet); no clock readout of SCNT embryos exists in the literature searched. Nothing
verifiable links Klf9 to reprogramming or rejuvenation, and the five genes are not in the Hendrickson 2017 DUX-target
list; neither claim is made. Word budget: abstract 250; main text trimmed elsewhere (summary paragraph, Oct4–Sox2
sentence, reassignment sentence, GSE162345 description, two gap sentences) so the submission version stays under 7,500.
Noted for proof stage: the published coefficient table lists 1,840 non-zero genes for this clock; the pickle used here
has 1,839.

```
8c92fc0b6373319fd99aad6ad9193fa7a06700ef60e195658d6b7d2cff5e1fe4  gz/results/MANUSCRIPT_gz_v1_EN.md
f610405e60145f137ff9d254cb4e413508fbd81b8b7ee5bdbc7a5849db65941a  gz/results/MANUSCRIPT_gz_v1_KR.md
e92ad48b8cbaa4527edfe35de608b30d58ca176600fa6f6ee6f209934ceaf759  gz/results/MANUSCRIPT_gz_v1_EN_submission.md
```

## 2026-09-19 — Coefficient count discrepancy resolved

The earlier note ("published table lists 1,840 non-zero genes; pickle has 1,839") was a miscount by the literature agent
that included the intercept row. Checked directly against Supplementary Table 5 of Tyshkovskiy et al. 2026
(`41586_2026_10542_MOESM7_ESM.xlsx`, sheet "(A) Composite clocks", column "Chronological Age, Multi-species,
Multi-tissue, Scaling"): 1,839 non-zero coefficients, the same 1,839 Entrez IDs as the model file used here, maximum
absolute difference 4.9 × 10⁻¹⁷, identical intercept 0.03212. The model used is exactly the published clock; nothing to
check at proof stage.

## 2026-09-20 — Figure canvas clipping found and fixed (all 11 figures checked)

While building the plain-language review page (`gz/results/review/easy_review_KR.html`), text was seen running off the
right edge of Figure 1d. All 11 rendered figures were then checked programmatically for ink within 3 px of the canvas
border (a border-touching pixel means content was cut off, since every panel has a margin). Six figures failed:

| Figure | What was cut | Fix in `gz/src/18_figures_publication.py` |
|---|---|---|
| 1d | footnote "…sim. maternal clearance only" | centred 2-line footnote → left-aligned 4-line footnote at the panel's left edge |
| 2a | "prespecified QC (3 of 76 excluded)" | 2-line note → 3 shorter lines; definition box widened to 132.5–165 mm and its text moved to x = 135 so the R line is no longer flush with the border |
| 3f | x-axis label "Clock genes in control order…" cut at the bottom | canvas height 224 → 228 mm |
| 4b | x-axis label "(oriented by the control change)" | label shortened to "(oriented to control)"; panel moved 141 → 139 mm, width 23 → 24 mm (at 136 mm its y-label collided with the colour-bar label, so 139 was chosen) |
| S1 | footnote "…95% bootstrap inte[rval]" | split into 2 lines, va='center' → 'top' at y = 111 |
| S4 | footnote "…95% boo[tstrap]" | 2 lines → 4 lines |

After re-rendering, the border check is clean for all 11 figures. No data, no analysis and no numbers changed — only text
wrapping, three panel/box coordinates and one canvas height. The affected figures were re-exported to `figures/pub/`,
`figures/supp/` and copied to `gz/results/submission/`. `SUBMISSION_CHECKLIST_KR.md` updated (Figure 3 is 167 × 228 mm).

This class of defect is invisible in the normal workflow because the clipped text is simply absent from the rendered
file; the border-ink check is now the way to catch it and should be re-run after any figure edit.

## 2026-09-20 — External fact-check of every literature-derived statement

All background claims in the Introduction and Discussion were checked against the cited papers (four parallel
literature audits: ageing clocks, ZGA/DUX, splicing, perturbation datasets and SCNT). Twenty corrections were applied
to `MANUSCRIPT_gz_v1_EN.md` and mirrored into the Korean version. Nine references were added, two removed, two fixed;
the list is now 58 entries and fully alphabetical, and every in-text citation resolves.

Wrong as written (now corrected):

1. **DBTMEE classes.** "All five induced contributors are classified as zygotically activated (major ZGA, two-cell
   transient or transitional)" was false. The database's own `cluster_gene_v2.tsv` (downloaded from
   dbtmee.hgc.jp/download/data/tables.tar.gz) gives Klf9 and Neto2 = "Major ZGA", Pi4k2a and Gpatch4 =
   "2-Cell Transient", **Psmb5 = "MGA"**. "Transitional" is not one of the ten DBTMEE classes. Each gene is now named
   with its actual class.
2. **"MYC-target" programme** attributed to Tyshkovskiy et al. 2026 does not occur in that paper, and the same
   paragraph reports the cell-cycle module clock as *elevated* in early development — the opposite of what we wrote.
   Replaced with the paper's actual statement (immune/lipid vs cell-cycle/splicing; inflammation and interferon are the
   top contributors to the decrease).
3. **Pig molecular-age negative claim** was false: telomere length has been measured across the morula-to-blastocyst
   transition in pig (Dang-Nguyen et al., 2012; also Li et al., 2023 with ZGA). The related "two species" claim for
   telomere lengthening was also wrong.
4. **Schaetzlein & Rudolph (2005) is a review**; the primary finding is Schaetzlein et al. (2004) PNAS 101:8034-8038.
5. **Higgins-Chen et al. (2022)** shows that individual CpG *measurements* are noisy (up to nine-year replicate
   deviations), not that clock *weights* are noisy. Cited correctly in both places.
6. **"Each embryonic application stayed within one species"** is false — Kerepesi et al. 2021 is itself a
   mouse-and-human embryo study. Changed to "at most two species per study".
7. **Kerepesi & Gladyshev (2023)** reports no change across human preimplantation and only the decrease to the
   epiblast; "a pattern later reported for human embryos" overstated it.
8. **GSE162345 timing.** Fusion is at single time points, 21 hpi (early two-cell) or 30 hpi (late two-cell), not
   "21-24 or 30-33 hpi" (verified in PMC8609233). Also disclosed: in the source study the two-cell arm is a comparator
   for a four-cell system, and at 30 hpi only 8 genes were activated from the donor genome against >1,000 at 21 hpi.
   Tomikawa et al. (2022), whose protocol covers four-cell embryos only, was removed as support for this design.
9. **GSE248499** was described as "nuclear transfer with histone demethylases"; the source study is a G9a-*inhibitor*
   study (G9a is a methyltransferase). The Kdm4d/Kdm3a arms do exist and are what we use; the description now says so.
10. **Li et al. (2025)** compares nuclear-transfer reprogramming across species, not zygotic genome activation;
    Oomen et al. (2025) — already our Gate 1 dataset — is the correct citation for cross-species ZGA.
11. **Hou et al. (2025)** does not mention epithelialisation and places OCT4/SOX2 function at the early inner cell
    mass; the sentence now attributes epithelialisation timing to Chandramohan et al. (2026) alone.
12. **A485 framing.** A-485 is a selective acetyl-CoA-competitive inhibitor of the p300/CBP HAT domain
    (Lasko et al., 2017; embryo pharmacology in Wang et al., 2022), not a transcription inhibitor, and in Xiao et al.
    (2025) its ZGA effect runs through failure to induce Dux and is bypassed by exogenous Dux. It is no longer
    presented as independent confirmation of a transcription requirement. For Obox3, the source study's evidence is
    gain-of-function and the OBOX requirement was shown for a six-gene knockout (Ji et al., 2023), so that leg was
    softened too. Only α-amanitin inhibits Pol II directly.
13. Smaller fixes: two-wave timing now cited to Aoki et al. (1997) and Abe et al. (2018) rather than Abe (2015,
    one-cell only) and Xiao (2025, mechanism); Dux-loss tolerance attributed to OBOX redundancy (De Iaco et al., 2020;
    Ji et al., 2023; Guo et al., 2024); Zscan4 restricted to the *late* two-cell embryo and a rare ES-cell
    subpopulation; Zhang et al. (2024) marked as a re-analysis of published RNA-seq; Deng et al. (2026) qualified
    (skipping-specific deficit; the "too much" arm is TDP43 overexpression); Jiang et al. (2026) article number added;
    Isaev & Knowles version pinned; GSE300734's unlinked-deposit status stated in the text.

Also corrected in the Korean version, which had lagged the reviewer-workflow revision: the cross-fitted fold ranges
(0.32-1.08, -0.56 to +0.15, 0.14-0.25, re-derived from `posthoc_gate4_crossfit_folds.tsv`) and the removal of the
r = -0.57 correlation.

Checked and found correct, with no change needed: Abe 2015/2018, Hendrickson 2017, Yu 2016, Sha 2020, Zhang 2024 (all
three perturbations, verified verbatim in the full text), Deng 2026 (including the bidirectional claim), Wyatt 2022,
Li 2024, Matoba 2014 (184/222 RRRs reactivated; 1,212 to 475 DEGs), Lanza 2000, Ogonuki 2002, Nie 2023, Sakamoto 2024
dataset composition, the three dataset-provenance citations, the five tool citations, and every DOI checked.
The bioRxiv prefix 10.64898 is openRxiv's legitimate prefix for preprints posted from 1 December 2025.

Still open: `MANUSCRIPT_gz_v1_EN_submission.md` and `SUPPORTING_INFORMATION.md` predate this revision and the
reviewer-workflow revision; the Korean version still carries the pre-renumbering supplementary figure and table
numbers; and the word budget (7,500 for Aging Cell) is exceeded and must be addressed before the package is rebuilt.

## 2026-09-20 — Figure tidying, typographic minus, and a selection volcano (Supplementary Figure S7)

Three kinds of change, no analysis re-run.

**Occlusion and clipping.** Figure 3f's legend sat on top of the three running-sum curves; it was moved above the axes
as a single row. Figure 4c's bottom legend was anchored at the panel centre and its first handle ran off the left edge;
it is now anchored 2 mm from the left. Figure 1b's per-species annotation sat on the shaded ZGA band and Figure 1c's
dataset note sat near a data point; both now carry a translucent white background. Figure 4a's footnote, extended to
state the selection rule, was split over two lines so it no longer reaches the right edge. `24_check_figure_bounds.py`
is clean for all 12 figures.

**Typographic minus.** Python writes negative numbers with an ASCII hyphen while matplotlib's axis formatter uses
U+2212, so annotations and tick labels disagreed throughout (for example "D -0.24" beside an axis reading "−0.25").
`fix_minus()` now rewrites every text artist before saving, replacing a hyphen with a minus only where it precedes a
digit and does not follow an alphanumeric, so "16-cell", "held-out" and "cross-fitted" are untouched.

**Supplementary Figure S7 (new).** For each primary dataset, every event passing the Gate 4 expression and missingness
filters is plotted as the control late-minus-early change in PSI against the gene-corrected empirical p value of that
change, with the prespecified thresholds (|ΔPSI| ≥ 0.10, p < 0.05) drawn and the selected events marked. This adds no
analysis: `25_figdata_volcano.py` copies the two columns from the Gate 4 outputs and asserts that the thresholds
reproduce the stored selection flag exactly (511, 1,376 and 77 events). Its purpose is to show the effect size and p
value behind a selection that the text previously stated only as counts. It is cited in Section 2.6 before the
principal-component display, so under the first-citation rule the volcano is S7 and the PCA display becomes S8.

Considered and not done: p values or FDR on the clock panels (Figures 1–3). Group sizes are two to four libraries, the
frozen plans specify effect sizes, bootstrap intervals and direction across datasets rather than significance tests,
and the central decomposition is an algebraic identity whose terms sum to the reported difference to 1e-9, so a p value
on it would not be meaningful. The Methods continue to state that the within-species permutation test is the only p
value used for inference.

## 2026-09-20 — Text-collision audit of the figures (new check, 26 real collisions fixed)

The border check of `24_check_figure_bounds.py` only sees ink at the canvas edge, so it could not detect a label
sitting on another label or spilling out of its own panel. `26_check_figure_overlaps.py` renders each figure in
memory and measures the bounding box of every text artist, reporting text-over-text, text-over-a-foreign-panel,
panel-over-panel and text-past-the-canvas, with positions in millimetres. Tick labels that matplotlib creates but
never draws (a tick outside the view, or any tick of a drawing canvas whose axis is switched off) are excluded, as is
text carrying an opaque background, which is a deliberate overlay.

First run: 26 real collisions.

| Figure | Collision | Fix |
|---|---|---|
| S1 | the 6-7 interval labels per panel, rotated 60°, overlapped each other by up to 8.2 pt | rotated to 90° (vertical), so the horizontal footprint is the cap height rather than the string length |
| S1 | y tick labels such as "0.0075" reached 1.2 mm into the neighbouring species panel — 22 mm panels on a 29 mm pitch leave only 7 mm | panels narrowed to 20 mm (9 mm gap) and the y locator reduced from 4 ticks to 3 |
| 3a | the stacked panels' y labels, "Contribution per gene" and "Running sum of contributions (= clock change, late − early two-cell)", touched at the left edge | the lower label shortened to "Running sum of contributions", matching panel f; the identity it stated is in the figure legend |
| 2 | the teal interaction label "I +0.17 [0.09, 0.26]" touched the footnote | footnote moved 2 mm down |

Both checks now pass for all 12 figures: no ink at the canvas border, and no text or panel collisions. Run them after
any figure edit; `26` is the stricter of the two and imports the figure module rather than reading the exported files,
so it needs no render to disk.

## 2026-09-20 — Panel-by-panel figure review and redesign

A reviewer-style pass over every panel, of the kind the author uses on other manuscripts. Six criticisms were judged
valid and acted on; two were judged wrong and are recorded here with the reason.

Acted on:

1. **Figure 2 was drawn on three different y scales** (P1 −0.33 to 0.03, P2 −0.10 to 0.16, P3 −0.26 to 0.06), so a
   drop of −0.24 and a drop of −0.10 occupied a similar length on the page. All three panels now share one scale, with
   tick labels on the first only. This is the most consequential change of the pass: the reader can now compare D by
   eye. The legend says the scale is shared.
2. **Figure 1d was the only bar chart left in the paper** after the earlier redesign. It is now points with bootstrap
   intervals, like every other panel.
3. **Figure 1a stated methods, not the question.** The "Readouts" box, which listed the clock, the decomposition and
   SUPPA2, is now a "The question" box: the clock falls inside the two-cell stage where no cell divides; does the fall
   depend on ZGA, and what is the value made of; six perturbation series and an exact decomposition under frozen
   plans. The box was enlarged from 22 to 28 mm so the text sits inside it.
4. **"major ZGA" and the largest-drop annotation were repeated in all four species panels of Figure 1b.** The band is
   now labelled once, as "major ZGA (literature)", and the per-species annotation carries only the interval.
5. **The grey explanatory paragraph beside Figure 3c** ("Contribution of a gene = clock coefficient x change in its
   preprocessed expression...") took up as much room as a panel and repeats the figure legend. Removed.
6. **Figures 3a and 3f are the same kind of cumulative curve** and a reader meets them twice without being told why.
   Each panel now carries a short title saying what it answers: "What the control decrease is made of" and "What the
   perturbed and rescued arms do to those same genes".

Judged wrong, not acted on:

- *"Figure 4a and 4b show the same 511 events twice."* They do not answer the same question: the heat map shows the
  per-library pattern and the cumulative distribution shows the shift with its median. The pair is the standard
  combination in splicing papers.
- *"Figure 1b and 1c are both mouse clock trajectories, so merge them."* The two are scored against different
  references (oocytes of the same species in GSE225056, zygotes in GSE45719) and cover different stages. Merging
  would require rescoring one of them against the other's reference, which would change a frozen Gate 1 / Gate 2b
  result after the fact. Left as two panels.

Left as a flagged judgement call: Supplementary Figure S8 (principal-component display) carries little information and
no test, and Figure 4c makes the same point quantitatively. It costs nothing in the main text and answers the
"is there global structure?" question in one glance, so it is kept for now.

Both figure checks pass for all 12 figures after the redesign.

## 2026-09-20 — Supplementary figures consolidated from eight to six

The author's review asked for fewer supplementary figures and for the splicing side to carry event-type detail.

- **Supplementary Figures S5 and S6 merged into one two-panel figure (S5).** Panel a is the expression of the five
  largest contributors, panel b the per-gene contribution profiles of the perturbed and rescued arms. Both are post-hoc
  detail on the same GSE280522 arms, so they belong together; the text now cites S5a and S5b.
- **The principal-component display was dropped.** It carried no test — the legend said so — and Figure 4c makes the
  same point per library and quantitatively. Its two citations and the Supporting Methods sentence about it were
  removed with it.
- **The splicing-selection volcano moved from S7 to S6**, keeping the numbering in first-citation order:
  S1 measured zygotic timing, S2 GSE66582, S3 quality control, S4 secondary series, S5 contributors and profiles,
  S6 selection volcano.
- **Figure 1a's "The question" box and Figure 2a's D / I / R definition box were removed** at the author's request;
  both definitions are in the figure legends and the Methods, so the boxes were duplication.
- **Figure 4 gained panel e**: for each kind of local event, the percentage of filtered events the control arm
  activates and, among those, the percentage whose PSI rises. This is Table S4 drawn, from the same source file as the
  volcano. Checked before the sentence was written: all seven event classes are represented in GSE280522 and
  GSE221985 but only five of seven in GSE300734 (no retained intron, no mutually exclusive exons), and 43-47% of the
  activated events rise in PSI. A first draft of the sentence said every class was represented in every dataset; that
  was corrected against the data before it entered the manuscript.

Both figure checks pass for the ten figures. `26_check_figure_overlaps.py` was updated to the new figure list.

## 2026-09-20 — Junction-level splicing started (post hoc, plan frozen first)

Limitation 3 of the manuscript — that splicing was quantified from pseudoalignment-based transcript estimates rather
than junction reads — is testable, so it is being tested. `plan/POSTHOC_junction_psi_frozen.md` (sha256
ad09490439c1a0101da9241e32acaf79962644722a29a73e268b54a5ab3dfead) was written and checksummed before any junction was
counted. It fixes the dataset (GSE280522 only, the 23 libraries already scored), the event class (skipped exons, the
only class with an unambiguous junction formula), the depth rule (>= 10 unique junction reads per library-event,
>= 20 events per library), the three readings with their thresholds, and the rule by which the illustrated events are
chosen. No outcome may change a gate verdict or a reported interaction.

Reference: Ensembl GRCm39 release 112 primary assembly and GTF — the same release used for the kallisto index and for
SUPPA2 event generation, so event coordinates are comparable. STAR 2.7.3a index built with `--sjdbOverhang 149`
(26 GB). Alignment is two-pass per sample, unique junctions only.

One convention hazard was handled in advance: SUPPA writes exon boundaries and STAR writes the first and last base of
the intron, a one-base difference. The script does not assume either; it tries both offsets against the annotated
junctions of the first library, takes whichever matches more, and records the choice. That decision does not depend on
any result.

Scale of the analysis, fixed before running: 4,076 skipped-exon events pass the Gate 4 filters in GSE280522, of which
155 are control-activated.

## 2026-09-21 — Junction-level splicing: result, and two implementation errors found on the way

All 23 GSE280522 libraries aligned (STAR 2.7.3a, two-pass, unique junctions; unique mapping 42.7–68.5%,
78k–242k junctions per library). The junction coordinate offset was chosen from the data as the plan specified:
+1 (277 of 400 annotated inclusion junctions matched at +1, none at 0), i.e. SUPPA's exon boundaries versus STAR's
first and last intron base, as expected.

**Reading 1 — agreement of the two quantifications: supports the Gate 4 quantification.** Pearson r = 0.779,
Spearman rho = 0.783 over 59,482 library-event pairs scored by both methods, against a prespecified threshold of
rho >= 0.6. The transcript-TPM PSI used for Gate 4 and a junction count of the same events agree.

**Reading 2 — reproduction: reproduces the Gate 4 direction.** Progression on junction PSI, control = 1 by
construction: control 1.00, A485 0.39, A485 + DUX 0.90, giving I = −0.61 and R = +0.51, against the prespecified
rule I < 0 and R > 0. The direction of both the interaction and the rescue holds under a different aligner, a
different statistic and a different quantification.

**Two errors in my implementation, both found after seeing a result, both corrected; recorded because of that.**

1. *Unbalanced event sets.* The first run scored each library on whatever events passed the depth rule in that
   library. Junction depth differs systematically by arm (1,499 to 3,491 events per library; unique mapping 44–59%
   in A485 late two-cell against 60–69% in control early two-cell), so libraries were being averaged over different
   events. The score is now computed on complete cases only — events scored in all 23 libraries — which is what the
   Gate 4 rule already required of its own quantification (PSI in >= 80% of the libraries of every arm-by-stage
   group).
2. *Division by a near-zero control change.* The score divides by the control arm's change. Of the 78
   control-activated SE events complete in every library, 33 had a junction control change below 0.02 and the
   smallest was exactly 0, so a handful of events with meaningless denominators dominated the mean. Gate 4 defined
   its events as |dPSI| >= 0.10 *in its own quantification*; the faithful translation applies the same threshold
   inside the junction quantification, which is now done.

With the first error alone the readings were I = +0.25, R = −0.38; with the first corrected but not the second,
I = +0.66, R = −0.47; with both corrected, I = −0.61, R = +0.51. The first two are artefacts of the two errors
above, not properties of the data, and are recorded here so the sequence is on the record.

**Limits of this test.** Only 31 events survive every restriction (skipped exons, control-activated, scored in all
23 libraries, junction control change >= 0.10) out of 155 control-activated skipped-exon events, so the magnitude
rests on few events. The junction score is in-sample by the same construction as the transcript score — events and
scale are set on the control arm and the control arm is then scored — so its magnitude is inflated in the same way
the in-sample transcript magnitude was (−0.44 against a cross-fitted −0.27) and is not cross-fitted here. The test
is read as a test of direction and of quantification method, not of magnitude. No gate verdict and no reported
interaction changes.

Outputs: results/posthoc_junction_psi.tsv (59,482 library-event rows with junction counts, junction PSI and TPM
PSI), results/posthoc_junction_psi.md. Events selected by rule for a read-level illustration: Camk2b-region
ENSMUSG00000028693 (dPSI −0.57), ENSMUSG00000071533 (+0.51), ENSMUSG00000026511 (−0.50).

## 2026-09-21 — Two figures added: the clock matrix (Figure 1e) and the junction validation (Supplementary Figure S7)

**Figure 1e — the matrix the clock value is a weighted sum of.** The published trajectory is one number per stage;
this panel draws the matrix behind it. `src/30_figdata_clock_matrix.py` re-runs the Gate 1 preprocessing unchanged on
mouse GSE225056 (75 libraries after the same quality-control rule) and writes the stage mean of the preprocessed
feature of all 1,839 clock genes, with each gene's coefficient and its contribution to the early-to-late two-cell
change. Rows are ordered by that contribution, so the blocks carrying the decrease sit at the top and bottom and the
middle is visibly empty.

Verification: the clock pipeline centres the features before the estimator (StandardScaler with `with_std=False`), so
a coefficient-weighted sum of uncentred features differs from the stored clock value by one constant, identical at
every stage. Stage-to-stage differences are therefore exact, and the script asserts them against
`results/gate1_embryo_tage.tsv`: maximum absolute difference 3.85e-16, constant offset −0.0079. The within-two-cell
change reproduces the stored −0.1907.

Worth recording as an independent observation: in GSE225056 the same interval decomposes into −1.549 over 464
downward genes and +1.359 over 448 upward genes, a net −0.19. The residual structure the paper reports for
GSE280522 (−1.40 and +1.17 giving −0.24) therefore also holds in the cross-species dataset, which was not part of the
decomposition analysis.

**Supplementary Figure S7 — junction-level validation.** Panel a is the density of junction PSI against the
transcript PSI used for Gate 4 over the 58,958 library-event pairs scored by both (r = ρ = 0.78). Panel b shows the
three events the frozen plan selects, as read coverage with junction arcs, in the late two-cell libraries of each arm.
The read-level picture matches the summary statistics: for ENSMUSG00000026511 the two inclusion junctions average
414 and 514 unique reads in control, 74 and 124 under A485, and 318 and 415 with DUX co-expression.

Source data: `figJ_junction_vs_tpm.tsv` (binned density, so the figure ships a table rather than 59,482 points) and
`figJ_sashimi.tsv`, written by `src/29_figdata_junction.py` from the STAR BAMs and SJ.out.tab files.

Figure count is now 4 main and 7 supporting. Both figure checks pass for all 11. Three collisions introduced by the
new panels were found by `26_check_figure_overlaps.py` and fixed: Figure 1d's footnote touching panel e's title, the
panel e value plot's y-label reaching into the contribution sidebar, and the Supplementary Figure S7 panel-a footnote
running into the sashimi tracks.

## 2026-09-21 — Story pass: a broken sentence, an unattributed interpretation, and six digressions removed

A reading of the whole manuscript for breaks in the argument, for background that rested on inference rather than a
source, and for material that leaves the central chain (the clock falls inside the two-cell stage → does that depend on
ZGA → what is the value made of). Fifteen changes; the reference list is now exactly at the journal's limit of 45 as a
consequence, not as the aim.

**Broken or unattributed:**

1. An earlier edit had left "DNA-methylation clocks indicated that epigenetic age declines ... and examined at
   single-cell resolution", whose subject had gone. Repaired.
2. "This coincidence invites the interpretation that the embryonic clock decrease reflects reprogramming" asserted a
   reading of the field without saying whose it was. It is Kerepesi et al.'s ground-zero framing, and is now
   attributed to them.
3. "Module- and system-resolved clocks answer this by retraining (Sehgal et al., 2025)" claimed that Systems Age
   answers the three limitations just listed (composition, feature noise, estimate uncertainty). It does not make that
   claim; the sentence was my inference. Replaced by a plain statement of what this paper does.
4. The Introduction summarised the maternal knockouts as "the same direction and wider intervals in two
   maternal-knockout datasets", which contradicts Section 2.3, where the Tardbp interaction reverses under V2. Now
   states the Brg1 result and the near-zero interaction separately.
5. One sentence in the Introduction carried three unrelated facts about DUX, Dux deletion and Zscan4. Split into three.

**Digressions removed — each left the central chain:**

6. Pan-tissue and causality-enriched methylation clocks (Lu et al., 2023; Ying et al., 2024), cited to show "how
   differently such models can be built and read". True, but not a step in this argument.
7. A methylation-entropy preprint placing a nadir at the human morula (Hao et al., 2026) — a different measurement,
   never used again.
8. Single-cell methylation age (Trapp et al., 2021), already covered by the two clock reports named beside it.
9. Preimplantation splicing dynamics and a programmed splicing failure at ZGA (Xing et al., 2020; Wyatt et al., 2022).
   The splicing background is carried by the two papers that describe zygotic splicing activation itself.
10. Splicing as an axis of adult ageing (Zhang, Tyshkovskiy, et al., 2026). An adult phenomenon and a different
    question from progression through ZGA.
11. Cross-species comparison of nuclear-transfer reprogramming (Li et al., 2025) — a different comparison.
12. The Oct4–Sox2 wave after the 8-cell stage and what it sets (Chandramohan et al., 2026; Hou et al., 2025). A later
    stage and a different question (timing against capacity); the paper never tests it. The open question that rested
    on it now asks whether the later decrease has the same structure, which is a question about our own data.
13. Clone molecular age (Lanza et al., 2000; Ogonuki et al., 2002). Adult clone telomeres and clone lifespan do not
    bear on a two-cell embryo clock value; reasoning from them to our SCNT result was a leap. The passage now states
    what Matoba et al. showed and that no clock readout of nuclear-transfer embryos exists, and leaves the sign open.

**Overstated link softened:**

14. "Across species, the post-hoc correlation ... points the same way" placed a cross-species correlation across 24
    intervals alongside a perturbation interaction as if they were the same kind of evidence. It now says that the
    correlation is consistent with a link but is not a perturbation test.

**Flow:**

15. The second Discussion paragraph was doing five jobs in about 700 words (relation to prior clock work, the
    cross-species negative result, A485 pharmacology, the splicing literature, the two readouts compared). Split at
    its seams into four. A literature-search statement in the Introduction was reworded so it reads as a search
    result rather than a claim about the field. The closing paragraph had "the same decomposition" twice.

After the pass: 45 references (limit 45), abstract 248 words (limit 250), main text with legends 8,607 (limit 7,500).
Every in-text citation resolves to an entry and every entry is cited; the list is alphabetical.

## 2026-09-21 — Trimmed to the journal's limits

Three passes, after the story pass had already removed the digressions.

1. **Numbers moved to the tables and figures that already hold them.** The four-way direction split of the
   decomposition (Table S3), the per-gene fold changes and loss percentages of the eight largest contributors
   (Figure 3b), the filtered-event denominators (Table S4), the per-species library counts (Table 1) and the
   read-length condition (Table 1) are cited rather than repeated.
2. **Design detail moved out of the Results.** The GSE162345 nuclear-transfer design now sits in Methods 4.7, where
   the other arm definitions are, and the Results state only why the series was moved to the secondary set.
3. **Methods compressed again and the remainder moved to the Supporting Methods**: the reference libraries per
   dataset, the Gate 1 library counts, the secondary-arm selection, the reassignment control's two numbers, the
   read-length condition with its per-dataset values, and the SUPPA2 settings.
4. **Two figure legends shortened** (Figures 1 and 3) by removing axis descriptions that the axes already carry.
5. **Four asides in the Discussion tightened**, including the clock's ageing-signature caution, which repeated a
   point made two sentences earlier.

Final counts: 45 references (limit 45), abstract 248 words (limit 250), and the main text — Introduction, Results,
Discussion and Methods — 7,147 words (limit 7,500). With the four main figure legends the total is 8,172.

Which of those two numbers the journal applies is genuinely ambiguous and is recorded here rather than resolved. The
live author guidelines say a research article "should be no more than 7,500 words in length, excluding Author
Checklist, Title Page, Tables, Figures, and References"; figure legends are named as counting only in the parallel
sentence for Reviews. The manuscript is compliant on the reading that "Figures" covers their legends and is 672 words
over on the reading that it does not. Cutting those 672 would mean removing content rather than repetition, so it is
left as a decision for submission.

## 2026-09-23 — A white band in Figure 1e turned out to be a reporting gap, not a rendering fault

The author asked why the middle of the Figure 1e heat map was white. It was not a colour-scale problem and the
values are not small: **927 of the 1,839 weighted clock genes are not detected in GSE225056 at all.** Their rows are
exactly zero, their contribution is exactly zero, and the sort by contribution therefore places all 927 together in
the middle. Drawing absence as white made a technical gap read as a biological statement ("the middle contributes
almost nothing"), which is how the panel and its footnote had been written.

Three consequences, all now in the manuscript.

1. **The panel.** Undetected genes are masked to grey, counted in the panel, and named in the footnote as genes the
   dataset does not detect that the clock pipeline median-imputes. The colour limit is now the 98th percentile of the
   detected values rather than of all values, which the zeros were dragging down. A first attempt printed the count
   inside the grey band; that puts text on the plotting area, so the three row blocks are instead labelled with
   brackets in a canvas to the left of the map, the convention Figure 4a already uses. The heat map moved to x 38
   and the value plot to x 108 to make room. Masking the undetected rows to grey was then dropped as well: it was
   honest but spent half the panel on rows carrying no data and contributing exactly zero to every value, and
   absence of data is a sentence rather than half a heat map. The panel now plots only the 912 genes the dataset
   detects, in two labelled blocks (464 pushing the value down, 448 up), and the 927 that are missing are stated in
   the footnote, in Methods 4.4 and in the Supporting Methods.

2. **Coverage was being reported in a way that overstates it.** The figure reported for Gate 1 — clock feature
   coverage 0.56–0.69 — is the fraction of *all* clock input features present, computed by `tage_py.feature_coverage`
   over the full feature list. Over the 1,839 genes that carry a non-zero coefficient, which are the genes that
   determine the value, coverage is:

   | dataset | weighted genes detected | |
   |---|---|---|
   | GSE225056 mouse, deposited 5′ counts | 912 / 1,839 | 0.50 |
   | GSE280522, kallisto | 1,355 / 1,839 | 0.74 |
   | GSE300734, kallisto | 1,509 / 1,839 | 0.82 |

   Both numbers are correct for what they measure, but only the first was in the paper. Methods 4.4 now gives the
   weighted-gene coverage and the Supporting Methods explain the difference. Absent genes are median-imputed and
   contribute nothing to a difference, so no reported value changes; what changes is the reader's ability to see that
   the cross-species comparison rests on half of the weighted clock. The Discussion now names this as a third reason
   that test was the weakest of the four, alongside the cleavage-division asymmetry and the pig interval.

## 2026-09-23 — Supplementary Figure S7b was not a sashimi plot

The junction panel drew straight-line triangles of constant width with no gene model, which is not the convention a
reader of a splicing paper expects. Redrawn properly: arcs are quadratic Béziers whose width scales with the square
root of the junction count, the three arms of an event share one coverage scale so the depth difference between arms
is real rather than normalised away, each event carries an exon/intron model with the alternative exon marked, and
the junction counts sit at the arc apices. Figure 1 grew to 254 mm and Supplementary Figure S7 to 156 mm; both
figure checks pass for all eleven figures.

## 2026-09-23 — How the perturbation series were chosen was never stated

The author asked whether there were problems with sample collection and whether there was a good reason to use these
particular samples. The answer existed in the frozen plan, the addenda and `metadata/perturb/`, but not in the
manuscript, which a reviewer would have asked for. Fourteen GEO series were screened against the frozen requirement
for a primary dataset — early *and* late two-cell libraries in the same study, in a perturbed arm *and* a matched
control arm:

- **three qualified** (GSE280522, GSE221985, GSE300734) and are the primary set;
- **two** (GSE248499, GSE235547) have a perturbation and its control at the late two-cell stage only, so they give
  single-stage contrasts and are secondary;
- **GSE162345** was listed as primary in the frozen plan but its deposit is not an early-to-late two-cell pair; moved
  to secondary before any score (addendum 2, 2026-09-16);
- **GSE195760** has both stages, but every two-cell arm is a nuclear transfer, so there is no control arm against
  which to form an interaction; descriptive only;
- **seven** (GSE166338, GSE196671, GSE214878, GSE229740, GSE262039, GSE269417, GSE298245) have no early/late
  two-cell split at all, so the drop cannot be computed.

Methods 4.1 now carries this in summary and the Supporting Methods carry the table. Collection problems already
disclosed elsewhere in the paper are unchanged: the GSE162345 design mismatch, the GSE45719 exclusions, the two
distinct SCNT controls in GSE248499, the three libraries failing prespecified quality control, GSE300734 carrying no
linked publication, and two to four libraries per group.

Body word count after these additions: 7,335 of 7,500.

## 2026-09-27 — Human replication of the decomposition: the structure holds, the gene identity does not

Plan `plan/POSTHOC_human_decomposition_frozen.md` (sha256 64315b27…), frozen before any human clock score, with
`addendum1` (sha256 ad08d2d7…) written after the first run failed in preprocessing and before any score existed.

**Why.** The central mouse result is that a scalar clock value is a small residual of large opposing per-gene
contributions. That is a claim about how this clock reads early embryos. If it is a property of the clock and not of
one mouse dataset, it should appear in human embryos at a different interval. The Discussion already named this as
the first open question.

**Data.** GSE36552 (Yan et al. 2013), 90 single cells summed to 20 embryo pseudobulks by the embryo number in the
deposited titles; the 34 embryonic stem cell samples were excluded by the plan. The deposit carries `Uniq_reads_num`,
unique read counts per gene symbol, alongside RPKM; the counts are used and the parse was checked against the total
read count in each file's header. E-MTAB-3929 (Petropoulos et al. 2016) was the first choice but EBI did not respond
from this machine; recorded in the plan.

**Addendum 1, the one correction.** The frozen plan asserted that human needs no ortholog mapping because it is a
training species. That is wrong about the data format: the clock's 10,487 features are *mouse* Entrez identifiers —
10,487 of 10,487 are in `Gene_table_mouse.csv` and 0 of 10,487 in the human, rat or monkey tables. Human symbols are
therefore carried symbol → human Entrez → mouse Entrez → mouse Ensembl using the tables shipped inside the clock
package, so `preprocess` runs unchanged. This is still better placed than the cow, pig and rabbit analyses of Gate 1,
which used one-to-one Ensembl orthologs assembled here: the human mapping is the correspondence the model was built
with. The plan's claim that human needs no mapping is withdrawn and must not be repeated in the manuscript.

**Reading 3 — how much of the clock this dataset carries.** 1,370 of the 1,839 weighted genes (0.74), against 912 in
mouse GSE225056, 1,355 in GSE280522 and 1,509 in GSE300734. Above the plan's floor of 300, so reading 2 is read.

**Reading 1 — is there a decrease?** 8-cell → morula, the interval fixed in the plan: **D = −0.098
[−0.232, +0.036]**, from 3 and 2 embryo pseudobulks. The direction matches the human decrease reported by
Zakar-Polyák et al. (2024) for this interval and it is the largest single drop in the trajectory here, but with two
morula embryos the interval includes zero. Stage means relative to oocytes: oocyte +0.026, zygote +0.031, 2-cell
+0.010, 4-cell +0.146, 8-cell +0.105, morula +0.007, late blastocyst −0.039.

**Reading 2 — is it a residual?** 696 genes contribute −1.243 and 674 contribute +1.145, summing to the reported
−0.098. **Residual fraction 0.08**, against 0.17 in GSE280522 and 0.12 in mouse GSE225056, and a plan threshold of
0.50. Prespecified reading: **same structure as mouse**. The cancellation is in fact more complete in human than in
either mouse dataset — the one-sided sums are twelve times the net.

**Follow-up, computed after the planned items and labelled as such.** It is *not* the same genes. Over the 1,098
clock genes with a non-zero contribution in both the human interval and the mouse GSE280522 control window, the
per-gene contributions are essentially uncorrelated (Pearson r = 0.072, Spearman rho = 0.112) and only 10 of the 50
largest downward contributors are shared. Two mouse studies measuring the same interval in the same species agreed at
r = 0.74. So the near-cancellation reproduces across species and intervals while the gene identity behind it does
not — which is what the manuscript already says about reading the largest contributors as a property of this clock in
a given window rather than as a ranking of genes, now shown rather than asserted.

**Limits.** Three 8-cell and two morula embryos; the interval on D includes zero; the interval spans a cleavage, so
it is not the division-free window the mouse result uses; single cells from 2013 Smart-seq summed to pseudobulks.
Nothing here changes a gate verdict or any mouse number.

Outputs: `results/posthoc_human_decomposition.{md,tsv}`, `results/human_embryo_tage.tsv`, `src/31_human_decomposition.py`.

---

## 2026-09-29 — The manuscript now leads with the near-cancellation, not with ZGA dependence

**Why.** Asked whether the manuscript was at the level the target journal expects, I said no, and the reason was where
the weight sat. The title and abstract led with ZGA dependence, and that claim rests on effectively one dataset:
GSE280522 gives I = +0.18 with an interval excluding zero, GSE300734 gives +0.11 with an interval that includes zero,
and GSE221985 gives +0.04 and flips sign under V2 — which is why the prespecified rule was not met. A reviewer reads
the title first, tests the weakest load-bearing claim, and the paper falls there. The decomposition finding does not
have that problem: it is an algebraic identity, it agrees between two independent studies at r = 0.74, and as of
2026-09-27 it holds in a second species at a different interval. So the claim that carries the most evidence is now
the claim in the title. Nothing was removed and no number changed; what changed is which result the paper is about.

**What changed.**

- Title: "Gene-level decomposition links the two-cell transcriptomic-age decrease in mouse embryos to zygotic genome
  activation" → "A transcriptomic clock reads early mouse and human embryos as a near-cancellation of opposing gene
  sets". Running title and alternative title follow.
- Abstract rewritten (250/250). It now opens on what the number is made of and closes on the two-species result.
- Introduction: the fifth question ("does the same structure appear in human embryos?") added to the list of what the
  study asks, and to the summary of what it finds.
- New Results 2.7, from the 2026-09-27 analysis, labelled post hoc with its plan frozen before any human score.
- Discussion now opens on the residual in both species, and a new paragraph separates what generalises (the
  cancellation) from what does not (the genes), with r = 0.07 against r = 0.74 as the contrast.
- Methods 4.13 describes the human pseudobulks and the three-step identifier mapping; the availability section moved
  to 4.14.
- Figure 5 added (167 × 72 mm): the human trajectory, the running sum over the prespecified interval, and the
  per-gene contributions against mouse. Five main figures now, still under the limit of six.
- The closing paragraph of the Discussion had listed "whether the same decomposition holds in human embryos" as an
  open question. It is no longer open; it is replaced by whether the genes differ because the species differ or
  because the intervals do.

**Trimming.** The reframing added about 620 words. Four passes took the body from 7,836 to 7,492 and the abstract
from 262 to 250, with references already at 45/45. The cuts were made where the same thing was said twice (the
Discussion paragraph that reopened with Section 2.6's own closing sentence), where a detail is stated elsewhere (the
three excluded run accessions, which the Figure S3 legend names), and where the cross-species gate and the
pharmacology of each perturbation are now setting rather than claim. In the submission draft the body is 7,357 words
because the availability text moves into its own statement section.

**Word-count ambiguity, recorded rather than resolved.** The journal excludes "Author Checklist, Title Page, Tables,
Figures, and References" from the 7,500. Whether "Figures" covers figure legends is not stated. The title page now
reports the two numbers separately — 7,357 body, plus 1,153 of legends and 966 of tables — so the editor can apply
their own rule. If legends must be counted, about 1,010 further words have to go.

**Rebuilt from the edited manuscript, so nothing can drift:** `MANUSCRIPT_gz_v1_EN_submission.md`,
`SUPPORTING_INFORMATION.md`, the three .docx, the figure bundle (12 figures as PDF, 600-dpi PNG and 600-dpi LZW
TIFF), `submission/figures/README.md`, `README_package.txt`, `SUBMISSION_CHECKLIST_KR.md`, and the two Korean review
pages that are built from the figures (`review/public_KR.html`, `review/easy_review_KR.html`). Both figure checks
pass: all 12 figures clear of the canvas border, no text or panel collisions.

**Two new scripts,** replacing steps that had been done by hand and had therefore gone stale (the figure README said
four main figures; `README_package.txt` still carried the old title):

- `src/33_export_submission_figures.py` — exports every rendered figure into the bundle and rewrites its README from
  the PDF page sizes and the manuscript's own legend titles. Its first version silently skipped Figure 3, whose
  legend title wraps across two lines, so the regex now crosses newlines and the script asserts that the figure
  numbers it found are 1..N and S1..SN with no gaps.
- `src/34_build_submission_package.py` — builds the three .docx with pandoc and rewrites `README_package.txt` from
  the manuscript title, the title-page counts and the contents of the bundle.

`build_easy_review.py` now re-encodes its embedded figures from `figures/pub` and `figures/supp` on every build, for
the same reason.

**Cover letter** rewritten to match: it now leads with the two-species near-cancellation, and no longer says "minor
zygotic genome activation", which the manuscript itself stopped claiming on 2026-09-19 (A485 cannot separate the two
waves).

**Not done.** The Korean manuscript (`MANUSCRIPT_gz_v1_KR.md`) is still at its 2026-09-20 state and now trails the
English one by the story pass, the fact-check corrections, the junction analysis, the dataset-selection section, the
coverage numbers, the supplementary renumbering, 58 → 45 references and this reframing. `review/manuscript_review_KR.html`
renders that file and was left unrebuilt rather than reprinting a stale source in a fresher-looking page.

---

## 2026-10-06 — The Korean manuscript becomes a summary, and the Korean pages become reproducible

**Why.** `gz/results/MANUSCRIPT_gz_v1_KR.md` was a full parallel translation that had stopped at its 2026-09-20
state. By today it trailed the English manuscript by the story pass, the fact-check corrections, the junction
analysis, the dataset-selection section, the coverage numbers, the supplementary renumbering, 58 → 45 references and
the 2026-09-29 reframing. Keeping two full manuscripts in step costs a rewrite on every edit, and the English one is
what gets submitted — so a stale translation is a standing risk of quoting a number the paper no longer reports.
Asked to choose, the author chose a summary.

**What exists now.** `gz/results/MANUSCRIPT_SUMMARY_KR.md` — the paper in Korean, about four pages: the one-paragraph
claim, what the introduction establishes, each Results section with its numbers, what the Discussion does and does not
claim, the gate table, what each of the five figures is for, the five limitations, the datasets, the submission
counts, and a file map. Every figure is stated as taken from the English submission draft, with the rule that the
English file wins if the two disagree.

`MANUSCRIPT_gz_v1_KR.md` keeps its 589 lines, with a header saying it was superseded on this date, what it is missing
and that it is kept as provenance of the translation. It was not deleted.

**`review/manuscript_review_KR.html` rebuilt** from the summary rather than from the retired file, with the five main
figures placed after the sections that discuss them. It had been sitting at its 2026-09-19 build, rendering a source
that was already stale.

**The three Korean page builders moved into the repository** (`gz/src/35`–`38`). They had been living in the session
scratchpad, which does not survive the session, so the pages the author actually reads were not reproducible. Paths
are now derived from `__file__` rather than hard-coded, and their intermediates (panel crops, re-encoded figures) go
to `gz/figures/review_assets/`, added to `.gitignore`. One crop had the wrong source height for Supplementary
Figure S7 (150 mm against the actual 156 mm), which is fixed by the move.

Pages: `public_KR.html` (general readers, nine questions), `easy_review_KR.html` (figure by figure, why each exists),
`manuscript_review_KR.html` (the summary with figures).

---

## 2026-10-06 — Inside the ground-zero framework: the stage-anchored reading of ZGA is dropped

**Prompt.** The author met V. N. Gladyshev, who said that ZGA may be a regulatory switch but that its relation to the age
minimum may differ between species, and pointed to the Xenopus preprint of his group (Zhang, Tarkhov, …, Peshkin &
Gladyshev, bioRxiv 10.1101/2022.08.02.502559, v1 posted 2022-08-04). I read the preprint in full (22 pages, figures
included; text saved to the session scratchpad). It has appeared since only as a conference abstract (Innovation in
Aging 7, Suppl. 1, igad104.2485, 2023) and is cited as a preprint.

**What the preprint reports.** A mammalian methylation array yields 1,068 usable CpGs in Xenopus laevis; a bagged
Elastic Net clock trained on 40 adult skin samples (185 CpGs, MAE 1.82 years) applied to embryos 5–31.5 hpf falls
rapidly around 10–12 hpf, at the onset of gastrulation (~9.5 hpf), where it reaches its minimum and then rises slowly;
methylation entropy is lowest at 10–12 hpf and global methylation plateaus at ~10 hpf; mean transcript abundance
(GSE73430) is lowest at NF9, just before gastrulation, and rises sharply during it; and of 5,546 embryos followed
individually, 104 (~2%) failed, 66 of them within 3 h of the onset of gastrulation. Ground zero is defined there by the
convergence of those four readouts, not by one clock. For mouse the preprint restates the E4.5–E9.5 window of
Kerepesi et al. (2021), whose own text gives E4.5–E10.5, most probably E6.5/E7.5.

**What it changes here.** Nothing numerical. Two things in how the paper reads itself were wrong or weak:

1. The Introduction said "the coincidence has been read as reprogramming resetting molecular age, the minimum being
   described as a ground zero" — conflating our two-cell window with the minimum, which the ground-zero papers place at
   gastrulation (mouse, frog), after implantation (human) and near E10 (tAge). The paper now says the two-cell decrease
   is a feature of the descending limb, at the first wholesale remodelling of the transcriptome, and is not the minimum.
2. Gate 1 asked whether the largest decrease coincided with a literature major-ZGA stage. Under "a switch whose position
   differs by species" that criterion asks different questions in different species (two-cell in mouse, 2→4-cell in
   pig, 8→16-cell in cow and rabbit). The verdict labels are untouched (1 of 4, mixed; 1-alt partial, post hoc), but
   the Results and Discussion now read the cross-species result as: the stage did not hold, the relation did — over the
   24 intervals the measured zygotic gain correlates with the decrease (ρ = 0.52), the pattern expected if the clock
   responds to the remodelling rather than to a landmark.

Added to the Discussion, after the human paragraph: a scalar that is a 17% or 8% residual of opposing flows moves
whenever either flow is touched (A485, DUX), so where its minimum falls through a remodelling event depends on the
balance of the flows as much as on the state the clock was trained to read; locating a transcriptomic ground zero
needs the decomposition, as the Xenopus study needed four readouts. The closing open questions now include whether the
descent into the post-implantation minimum has the same structure — testable with the same identity.

**References.** Added Gladyshev (2021, Trends Mol Med) and Zhang et al. (2022, bioRxiv); dropped Schaetzlein et al.
(2004), Dang-Nguyen et al. (2012) and Falco et al. (2007) together with the telomere and Zscan4 sentences they
supported. 44 of 45. The reference-check note at the end of the manuscript records this.

**Word budget.** The reframing added ~330 words; three passes took them back out — mainly where a point was made twice
(the ρ = 0.52 sentence that stood in both D2 and D3; the identity's provenance in both 2.4 and D6; the GSE162345
move in 2.3, 4.1 and 4.2), descriptive clauses nothing used (the ICM libraries of GSE66582; the eight-against-a-thousand
genes of GSE162345, which stays in this log) and a source-study module comparison. Working file 7,482 / 7,500;
submission draft 7,347 (availability text in its own section); abstract 250 / 250.

**Not done, proposed instead.** The framework suggests one analysis the identity can do directly: decompose the
descent into the ground-zero minimum itself — human blastocyst → post-implantation epiblast (GSE109555, Zhou et al.
2019, 65 peri-implantation embryos; GSE136447, Xiang et al. 2020, 3D-cultured embryos to the primitive-streak anlage)
and mouse E3.5 → E6.5 (GSE100597, Mohammed et al. 2017) — and ask whether it is also a near-cancellation. Deposited
counts would suffice, as for GSE36552 and GSE45719. It needs a frozen plan and the author's decision first; nothing was
downloaded.

**Derived files rebuilt** from the manuscript: submission text (src/39), supporting information (src/40), the three
.docx and README (src/34), and the three Korean pages (src/35–38). The submission-text and SI builders had been
living in the session scratchpad with the Crossref cache; they now sit in `gz/src/` with the cache at
`gz/metadata/refs_crossref.json`.

**For external feedback.** `gz/results/for_feedback_VNG/` holds a one-page note in English (what was done, where it
sits relative to ground zero, four questions), the two .docx files, the 12 figures as PDF and four data tables. The
folder is excluded from the public repository because it contains the manuscript. Nothing has been sent.

---

## 2026-10-08 — The feedback package is cut down to the story and the solid numbers

**Why.** The author wants to send V. N. Gladyshev the argument first and the full data only if asked, and wants the
points on which his advice is sought set out separately.

**What is in `gz/results/for_feedback_VNG/` now.** `send/`: a one-page note (`NOTE_for_VNG.md` — the argument in
six steps, the reading inside the ground-zero framework, three lines on what is not yet solid, and the points for his
view grouped as A interpretation / B method / C species and data / D positioning) and one figure (`story_figure.pdf`,
`src/41_story_figure_for_feedback.py`: the mouse control running sum, the three arms over the same genes, the human
running sum, and human against mouse per gene — all from the published source-data tables, nothing recomputed), plus
Figures 2, 3 and 5 as optional attachments. `hold/`: the manuscript and SI .docx, the twelve figures, the four data
tables and the earlier long note, to be sent on request. A Korean guide (`보내기_안내_KR.md`) says which numbers went
into the note and why, which were left out, and lists the A–D points in Korean. The folder stays out of the
repository. Nothing has been sent.

**Numbers admitted to the note** — each either an identity, a prespecified verdict, or replicated: the within-two-cell
decrease in three datasets; −1.40 / +1.17 / 17%; r = 0.74 between studies; I = +0.18 [0.11, 0.25]; 85% / 113% with
the scalar unmoved; the splicing direction surviving cross-fitting; the human 8% residual and r = 0.07. The three
weaknesses that could change his advice are stated (cross-species stage test 1 of 4; the smallest dataset's sign flip
in one variant; human n = 5 with an interval including zero). Gene names, the secondary series, fold values, the
dataset-selection account and the methods detail were left out.

**Later the same day.** The note is an opinion request only — no coauthorship, acknowledgement or reviewer matters are
raised, at the author's instruction. Questions that exist because our own data are thin are marked ▲ in the note
(A1 no post-implantation window of our own; B3 half-covered atlas; C1 human n = 5; C2 perturbation groups of 2–4 with
the verdict on an n = 2 dataset; C3 the 24-interval ρ = 0.52; C4 no E3.5 → E8.5 series), and the Korean guide tabulates
what data would resolve each.
