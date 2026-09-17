# 2C_ZGA — transcriptomic age across the mouse two-cell stage

Analysis code, prespecified plans and audit trail for the manuscript in preparation
*"The two-cell transcriptomic-age decrease in mouse embryos is a residual of opposing gene contributions"*.

The study applies a published multi-species transcriptomic clock (tAge; Tyshkovskiy et al., 2026) to public
preimplantation RNA-seq data and asks whether the clock decrease within the two-cell stage depends on zygotic genome
activation (ZGA). No new data were generated and no clock was trained.

> **Status.** Manuscript in preparation; results are not peer reviewed. The manuscript draft is not included.

## What is here

| Path | Content |
|---|---|
| `gz/plan/` | Analysis plans frozen before scoring (Gates 1–4), dated addenda, post-hoc plans written before measurement |
| `gz/src/` | All analysis and figure scripts, numbered in execution order |
| `gz/results/` | Gate decisions (`*_decision.md`), post-hoc reports and result tables |
| `gz/figures/pub/`, `gz/figures/supp/` | Figures 1–5 and Supplementary Figures S1–S7 (PDF and 600 dpi PNG, 167 mm wide) |
| `gz/figures/source_data/` | Every value drawn in the figures |
| `gz/metadata/` | GEO/ENA sample tables; `perturb/gate3_accessions_geo_sra.tsv` maps every run to GSM, SRX and BioSample (Table S1) |
| `metadata/tx2gene_full.tsv` | Ensembl GRCm39 release-112 transcript-to-gene map (cDNA + ncRNA, 145,918 transcripts) |
| `logs/AUDIT_LOG.md` | Dated record of every decision, deviation, failure and checksum, including the project's earlier direction |
| `logs/*.log` | Console logs of the Gate 3–4 and post-hoc runs |

Labels used throughout: **pre** = analysis frozen before any score was seen; **post hoc** = planned after a gate
verdict (each post-hoc plan was still written before its own measurement and is labelled in the script docstring).

## Scripts

| Script | Step | Plan |
|---|---|---|
| `00_validate_tage_py.py` | Python tAge port vs the package's golden values | — |
| `01_gate1_cross_species.py`, `02_fig_…` | Gate 1: timing of the clock decrease vs major ZGA in four species (GSE225056) | `GATE1_PLAN_frozen.md` + addendum 1 |
| `03_gate1alt_zga_timing.py`, `04_fig_…` | Gate 1-alt: measured zygotic timing (post hoc) | `GATE1ALT_PLAN_frozen.md` |
| `05_gate2a_composition.py` | Gate 2a: composition artefact test | `GATE2A_PLAN_frozen.md` |
| `06_gate2b_R1_deng.py` | Gate 2b R1: GSE45719 | `GATE2B_PLAN_frozen.md` |
| `07_gate2b_R2_quant.sh`, `08_gate2b_R2_wu.py` | Gate 2b R2: GSE66582 download, kallisto, clock | `GATE2B_PLAN_frozen.md` |
| `09_fig_gate2.py` | Gate 2 review figure | — |
| `10u_make_worklists.py`, `10u_gate3_download_quant.sh` | Gate 3 downloads (ENA, md5-verified) and kallisto (76 runs); `10_…` is the first version | `GATE3_PLAN_frozen.md` |
| `11_gate3_analysis.py` | Gate 3: interaction of the two-cell drop with ZGA perturbation | `GATE3_PLAN_frozen.md` + addenda 1–2 |
| `12u_quant_concordance.py` | Technical concordance of repeated quantifications | — |
| `13_posthoc_gate3_contribution.py` | Exact per-gene decomposition of the drop (post hoc) | `POSTHOC_gate3_contribution_frozen.md` |
| `14_posthoc_rescue.py` | Rescue arm: zygotic programme and contributions (post hoc) | `POSTHOC_rescue_frozen.md` |
| `15_gate4_analysis.py` | Gate 4: SUPPA2 splicing progression | `GATE4_PLAN_frozen.md`, `GATE4_IMPLEMENTATION_NOTES.md` |
| `16_posthoc_gate4_crossfit.py` | Cross-fitted check of the Gate 4 score (post hoc) | `POSTHOC_gate4_crossfit_frozen.md` |
| `17_figdata.py`, `19_figdata_supp.py` | Figure source data; each value is checked against the stored result before writing | — |
| `18_figures_publication.py` | Figures 1–5 and S1–S7 from `source_data/` only | — |

Scripts locate the project from their own path (`analysis/` = three levels above the script), so the repository can
be cloned anywhere. Run them from any directory, e.g. `python gz/src/11_gate3_analysis.py`.

## Requirements

* Python 3.11 with the pinned packages in `gz/requirements_gz.txt` (scikit-learn 1.5.2 matters: the clock pickles were
  validated with it).
* kallisto 0.51.1 on `PATH` (or set `KALLISTO=/path/to/kallisto`).
* SUPPA2 (Trincado et al., 2018; https://github.com/comprna/SUPPA) cloned to `gz/tools_SUPPA`, with one import changed
  for current statsmodels in `lib/diff_tools.py`:
  `statsmodels.sandbox.stats.multicomp.multipletests` → `statsmodels.stats.multitest.multipletests`.
* tAge (https://github.com/Gladyshev-Lab/tAge, v1.1.0, commit 0dba58f) cloned to `gz/tools_tAge`, and the Elastic Net
  clock files (Zenodo record 18763485) placed in `gz/tools_tAge_models/`. Both are distributed by their authors under
  the MGB Open Access License and are not redistributed here.
* **`gz/src/tage_py.py` (Python port of the tAge preprocessing) is not yet included.** The tAge licence asks that
  modified versions be shared back with the licensor before distribution; the port will be added after that step.

## Data

Nothing under `data/` is tracked (about 200 GB). To rebuild:

1. **Reference.** Ensembl release 112 `Mus_musculus.GRCm39.cdna.all.fa.gz` and `ncrna.fa.gz`, concatenated (cDNA then
   ncRNA) into `data/raw/mm_cdna_ncrna.fa.gz`; `kallisto index -t 8 -i data/processed/mm_full.idx` on that file
   (750,479 contigs, 145,918 targets; `logs/kallisto_index_full_ubuntu.log`). `src/50_tx2gene_full.py` builds the
   transcript-to-gene map shipped in `metadata/` (its final comparison block refers to a superseded map that is not
   distributed). The GTF `Mus_musculus.GRCm39.112.gtf` goes to `gz/data/annotation/`.
2. **Deposited matrices.** GSE225056 5′ count files and GSE45719 per-cell expression files from GEO into `gz/data/`.
3. **Raw reads.** The 76 perturbation runs in `gz/metadata/perturb/gate3_selected_runs.tsv` (GSE280522, GSE221985,
   GSE300734, GSE162345, GSE248499, GSE235547) and the 15 GSE66582 runs, fetched by the download scripts above. Read
   counts of all 76 runs match the NCBI SRA records (`gate3_accessions_geo_sra.tsv`).
4. **Orthologs.** NCBI `gene_orthologs.gz` (https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz) into `gz/metadata/`
   (not tracked: over the GitHub file-size limit).

Checksums of inputs, frozen plans and outputs are recorded in `logs/AUDIT_LOG.md`.

## Licence

No licence has been chosen yet; until one is added, all rights are reserved by the author.
