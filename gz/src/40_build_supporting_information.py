#!/usr/bin/env python3
"""Regenerate SUPPORTING_INFORMATION.md from the manuscript, so numbering can never drift."""
import re

import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
RES = f'{B}/results'
t = open(f'{RES}/MANUSCRIPT_gz_v1_EN.md', encoding='utf-8').read()
methods = open(f'{RES}/SUPPORTING_METHODS.md', encoding='utf-8').read().strip()

title = re.search(r'^# (.+)$', t, re.M).group(1).strip()

# supporting figure legends, in manuscript order
blk = t.split('## Supporting Information figure legends', 1)[1].split('\n---\n', 1)[0]
legends = re.split(r'\n\n(?=\*\*Figure S\d:)', blk.strip())
legends = [' '.join(l.split()) for l in legends if l.strip().startswith('**Figure S')]

# supporting tables section
tables = t.split('## Supporting Information tables', 1)[1].split('\n---\n', 1)[0].strip()

out = [f'# Supporting Information', '', f'**{title}**', '', '## Contents', '',
       'Supporting Methods; Supporting Figures S1–S7 with legends; Supporting Tables S1–S5. Table S1 (run-to-sample',
       'accession table) and Tables S2, S3 and S5 are supplied as separate tab-separated files and as one workbook',
       '(`SupplementaryTables_S2-S5.xlsx`); Table S4 is reproduced below.', '',
       '---', '', methods, '', '---', '', '## Supporting Figures', '']

for i, leg in enumerate(legends, start=1):
    out += [f'### Figure S{i}', '', f'![Figure S{i}](submission/figures/png/FigureS{i}.png)', '', leg, '']

out += ['---', '', tables.replace('## Supporting Information tables', '').strip(), '']

open(f'{RES}/SUPPORTING_INFORMATION.md', 'w', encoding='utf-8').write('\n'.join(out))
wc = lambda s: len(re.findall(r"[A-Za-z0-9α-ωΔ][\w'’.\-/]*", s))
print('supporting figures found:', len(legends))
for i, l in enumerate(legends, 1):
    print(f'  S{i}: {l[:78]}')
print('SI words:', wc('\n'.join(out)))
