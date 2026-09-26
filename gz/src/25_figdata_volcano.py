"""Source data for the splicing-selection volcano (Supplementary Figure S8).

Every event that passed the Gate 4 expression and missingness filters in a primary dataset, with the two quantities
the frozen plan used to select control splicing-activation events: the control late-minus-early change in PSI and the
gene-corrected empirical p value from SUPPA2 diffSplice. No new analysis: the columns are copied from the Gate 4
outputs and the script asserts that the selection flag is exactly |dPSI| >= 0.10 and p < 0.05.
"""
import os

import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
B = f'{ANALYSIS}/gz'
GSE = {'GSE280522': 'P1', 'GSE221985': 'P2', 'GSE300734': 'P3'}

out = []
for gse, tag in GSE.items():
    d = pd.read_csv(f'{B}/results/gate4_events_filtered_{gse}.tsv', sep='\t')
    d.columns = ['event', 'event_type', 'dPSI_control', 'p_diffsplice', 'zsa']
    rule = (d.dPSI_control.abs() >= 0.10) & (d.p_diffsplice < 0.05)
    assert (rule == d.zsa.astype(bool)).all(), gse      # the plotted thresholds are the selection rule
    d.insert(0, 'dataset', tag)
    d.insert(1, 'gse', gse)
    out.append(d)

R = pd.concat(out, ignore_index=True)
R.to_csv(f'{B}/figures/source_data/figS8_volcano_events.tsv', sep='\t', index=False)
print('wrote figS8_volcano_events.tsv:', len(R), 'events')
print(R.groupby('dataset').agg(filtered=('event', 'size'), selected=('zsa', 'sum'),
                               min_p=('p_diffsplice', 'min'), max_p=('p_diffsplice', 'max')).to_string())
