"""Build the .docx files of the submission package and rewrite its README.

Inputs (all generated earlier, none edited here):
  results/MANUSCRIPT_gz_v1_EN_submission.md   main text with title page and statements
  results/SUPPORTING_INFORMATION.md           supporting methods, figures and tables
  results/COVER_LETTER_draft.md               cover letter

Outputs under results/submission/: the three .docx files and README_package.txt, whose title, figure counts and
file list are read from the manuscript and the bundle rather than typed, so they cannot go stale.

Run with /usr/bin/python3 after 33_export_submission_figures.py.
"""
import os
import re
import subprocess

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
RES = f'{B}/results'
SUB = f'{RES}/submission'
PANDOC = '/home/shim/anaconda3/bin/pandoc'

ms = open(f'{RES}/MANUSCRIPT_gz_v1_EN_submission.md', encoding='utf-8').read()
title = re.search(r'^\*\*Title:\*\* (.+)$', ms, re.M).group(1).strip()
counts = ' '.join(re.search(r'^\*\*Counts:\*\* (.+?)\n\n', ms, re.M | re.S).group(1).split())

# ---- the three .docx, built from the markdown with pandoc (cwd = results, so the image paths resolve)
for src, dst in [('MANUSCRIPT_gz_v1_EN_submission.md', 'Manuscript_AgingCell.docx'),
                 ('SUPPORTING_INFORMATION.md', 'Supporting_Information.docx'),
                 ('COVER_LETTER_draft.md', 'Cover_letter.docx')]:
    subprocess.run([PANDOC, src, '-o', f'submission/{dst}'], cwd=RES, check=True)
    print(f'{dst}  {os.path.getsize(f"{SUB}/{dst}"):,} bytes')

# ---- README_package.txt, listing what is actually in the bundle
fig = open(f'{SUB}/figures/README.md', encoding='utf-8').read()
nmain = len(re.findall(r'^\| Figure\d+ \|', fig, re.M))
nsupp = len(re.findall(r'^\| FigureS\d+ \|', fig, re.M))
nsupp_tab = len([f for f in os.listdir(SUB) if re.match(r'TableS\d+_', f)])

README = f"""Submission package — Aging Cell, Original Research Article
{title}

{counts}

Manuscript_AgingCell.docx       main text: title page, abstract, introduction, results, discussion, methods,
                                statements, references, Tables 1-2, figure legends, supporting-information legends
Cover_letter.docx               cover letter (author details still to be filled in)
Supporting_Information.docx     Supporting Methods, Supporting Figures S1-S{nsupp} with legends, Supporting Tables

figures/                        every figure as its own file ({nmain} main, {nsupp} supporting)
  figures/pdf/                  vector PDF, the preferred upload format
  figures/png/                  600 dpi PNG
  figures/tiff/                 600 dpi LZW TIFF, for systems that require TIFF
  figures/README.md             one line per figure with its size and title

TableS1_run_to_sample.tsv       Table S1: every run mapped to GEO sample, SRA experiment, BioSample, layout, reads
TableS2_clock_values.tsv        Table S2: clock value of every scored embryo and library
TableS3_gene_contributions.tsv  Table S3: per-gene contributions to the within-two-cell drop
TableS4_splicing_events.tsv     Table S4: splicing events in the control arms
TableS5_crossfit_folds.tsv      Table S5: cross-fitted splicing progression per fold
SupplementaryTables_S2-S5.xlsx  Tables S2-S{nsupp_tab} as one workbook with a README sheet

Still to complete before upload: authors, affiliations, ORCID, corresponding author, CRediT roles, funding,
conflict-of-interest and acknowledgement statements, optional suggested reviewers, and the archived-code DOI.
Source data for every figure: analysis/gz/figures/source_data/ in https://github.com/vincentjshim-rgb/2C_ZGA
"""
open(f'{SUB}/README_package.txt', 'w', encoding='utf-8').write(README)
print(f'\nREADME_package.txt rewritten\n  {title}\n  {counts}')
