"""POST HOC (plan: plan/POSTHOC_junction_psi_frozen.md). Junction-level PSI for the skipped-exon events of
GSE280522, from STAR splice junctions, and the three readings the plan fixed before any junction was counted:

  1. agreement between junction PSI and the transcript-TPM PSI used for Gate 4
  2. whether the Gate 4 interaction and rescue difference reproduce on junction PSI
  3. the three events selected by rule for the read-level illustration

Nothing here changes a gate verdict. seed 20260921.
"""
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
GSE = 'GSE280522'
ALIGN = f'{B}/align/{GSE}'
MIN_READS = 10          # plan: an event is scored in a library only at >= 10 unique junction reads
MIN_EVENTS = 20         # plan: a library needs >= 20 scorable events
NBOOT, SEED = 2000, 20260921

# ---------------------------------------------------------------- events and their junctions
E = pd.read_csv(f'{B}/results/gate4_events_filtered_{GSE}.tsv', sep='\t')
E.columns = ['event', 'event_type', 'dPSI_control', 'p_diffsplice', 'zsa']
SE = E[E.event_type == 'SE'].copy()
coord = SE.event.str.extract(r';SE:([^:]+):(\d+)-(\d+):(\d+)-(\d+):([+-])')
coord.columns = ['chrom', 'e1', 's2', 'e2', 's3', 'strand']
SE = pd.concat([SE.reset_index(drop=True), coord.reset_index(drop=True)], axis=1).dropna()
for c in ('e1', 's2', 'e2', 's3'):
    SE[c] = SE[c].astype(int)
print(f'{len(SE)} skipped-exon events in the filtered set ({int(SE.zsa.astype(bool).sum())} control-activated)')

# ---------------------------------------------------------------- STAR junctions per library
Mc = pd.read_csv(f'{B}/figures/source_data/fig4a_heatmap_columns.tsv', sep='\t')
Mc = Mc[['run', 'arm', 'stage']].drop_duplicates()
runs = [r for r in Mc.run if os.path.exists(f'{ALIGN}/{r}/SJ.out.tab')]
assert runs, 'no SJ.out.tab found — run src/27_junction_align.py first'
print(f'{len(runs)} aligned libraries')

SJ = {}
for r in runs:
    d = pd.read_csv(f'{ALIGN}/{r}/SJ.out.tab', sep='\t', header=None,
                    names=['chrom', 'start', 'end', 'strand', 'motif', 'annot', 'n_unique', 'n_multi', 'overhang'],
                    dtype={'chrom': str})
    SJ[r] = dict(zip(zip(d.chrom, d.start, d.end), d.n_unique))

# SUPPA writes exon boundaries; STAR writes the first and last base of the intron. Decide the offset
# from the data rather than assuming it: whichever convention finds more of the annotated junctions.
ref = SJ[runs[0]]
hits = {}
for off in (0, 1):
    h = sum((row.chrom, row.e1 + off, row.s2 - off) in ref for row in SE.head(400).itertuples())
    hits[off] = h
OFF = max(hits, key=hits.get)
print(f'junction coordinate offset chosen from the data: +{OFF} ({hits})')


def counts(row, sj):
    i1 = sj.get((row.chrom, row.e1 + OFF, row.s2 - OFF), 0)
    i2 = sj.get((row.chrom, row.e2 + OFF, row.s3 - OFF), 0)
    s = sj.get((row.chrom, row.e1 + OFF, row.s3 - OFF), 0)
    return i1, i2, s


rows = []
for r in runs:
    sj = SJ[r]
    for row in SE.itertuples():
        i1, i2, s = counts(row, sj)
        tot = i1 + i2 + s
        if tot < MIN_READS:
            continue
        inc = (i1 + i2) / 2
        rows.append((r, row.event, i1, i2, s, inc / (inc + s)))
J = pd.DataFrame(rows, columns=['run', 'event', 'inc1', 'inc2', 'skip', 'psi_junction'])
per_lib = J.groupby('run').size()
keep = per_lib[per_lib >= MIN_EVENTS].index
J = J[J.run.isin(keep)]
print(f'{len(J):,} library-event pairs scored in {J.run.nunique()} libraries '
      f'({J.event.nunique():,} events, median {int(per_lib.median())} per library)')

# ---------------------------------------------------------------- transcript-TPM PSI for the same pairs
P = pd.read_csv(f'{B}/data/gate4/{GSE}.psi', sep='\t', index_col=0)
P.columns = [c.split('/')[-1].replace('.psi', '') for c in P.columns]
T = P.reindex(J.event.unique())
tpm = T.stack(dropna=True).rename('psi_tpm').reset_index()
tpm.columns = ['event', 'run', 'psi_tpm']
J = J.merge(tpm, on=['event', 'run'], how='left')
J = J.merge(Mc, on='run', how='left')
J.to_csv(f'{B}/results/posthoc_junction_psi.tsv', sep='\t', index=False)

ok = J.dropna(subset=['psi_tpm'])
r_p = pearsonr(ok.psi_junction, ok.psi_tpm)
r_s = spearmanr(ok.psi_junction, ok.psi_tpm)
reading1 = ('supports the Gate 4 quantification' if r_s.statistic >= 0.6 else
            'partially supports the Gate 4 quantification' if r_s.statistic >= 0.4 else
            'does not support the Gate 4 quantification')

# ---------------------------------------------------------------- progression on junction PSI
rng = np.random.default_rng(SEED)
ctrl = J[J.arm == 'control']
piv = J.pivot_table(index='event', columns='run', values='psi_junction')
ce = [r for r in Mc.run[(Mc.arm == 'control') & (Mc.stage == 'E2C')] if r in piv.columns]
cl = [r for r in Mc.run[(Mc.arm == 'control') & (Mc.stage == 'L2C')] if r in piv.columns]
zsa_se = set(SE.loc[SE.zsa.astype(bool), 'event'])
# The score is a mean over events, so every library must be averaged over the SAME events: junction depth
# differs by arm (see the report), and scoring each library on whatever it happens to cover would compare
# different event sets. This mirrors the Gate 4 rule, which required PSI in >= 80% of the libraries of every
# arm-by-stage group. Complete cases only.
complete = piv.index[piv.notna().all(axis=1)]
cand = [e for e in complete if e in zsa_se]
base_all = piv.loc[cand, ce].mean(axis=1)
d_all = piv.loc[cand, cl].mean(axis=1) - base_all
# The score divides by the control change, so the control change has to be a real one. Gate 4 defined its
# events as |dPSI| >= 0.10 in its own quantification; the faithful translation applies the same threshold
# inside this quantification. Without it the mean is dominated by events whose junction control change is
# ~0 (33 of 78 below 0.02, minimum exactly 0), which is a division artefact, not a measurement.
DPSI_MIN = 0.10
ev = [e for e in cand if abs(d_all[e]) >= DPSI_MIN]
base = base_all[ev]
d = d_all[ev]
sgn = np.sign(d)
scale = d.abs()


def progression(run):
    v = piv.loc[ev, run]
    return float((sgn * (v - base) / scale).mean(skipna=True))


prog = {r: progression(r) for r in piv.columns}


def arm_mean(arm, stage):
    rs = [r for r in Mc.run[(Mc.arm == arm) & (Mc.stage == stage)] if r in prog]
    return np.mean([prog[r] for r in rs]), rs


res = {}
for arm in ['control', 'A485', 'A485+Dux']:
    e_m, e_r = arm_mean(arm, 'E2C')
    l_m, l_r = arm_mean(arm, 'L2C')
    res[arm] = (l_m - e_m, e_r, l_r)
c0 = res['control'][0]
scaled = {a: (v[0] / c0 if c0 else np.nan) for a, v in res.items()}
I = scaled['A485'] - scaled['control']
R = scaled['A485+Dux'] - scaled['A485']
reading2 = 'reproduces the Gate 4 direction' if (I < 0 and R > 0) else 'does not reproduce the Gate 4 direction'

# ---------------------------------------------------------------- events for the read-level illustration
depth_ok = J.groupby('event').run.nunique() == J.run.nunique()
cand = SE[SE.zsa.astype(bool) & SE.event.isin(depth_ok[depth_ok].index)]
top3 = cand.reindex(cand.dPSI_control.abs().sort_values(ascending=False).index).head(3)

with open(f'{B}/results/posthoc_junction_psi.md', 'w') as fh:
    fh.write('# Post hoc: junction-level PSI in GSE280522 (plan: plan/POSTHOC_junction_psi_frozen.md)\n\n')
    fh.write(f'Skipped-exon events in the Gate 4 filtered set: {len(SE)} '
             f'({int(SE.zsa.astype(bool).sum())} control-activated).\n')
    fh.write(f'Libraries aligned and scored: {J.run.nunique()} of 23. '
             f'Library-event pairs at >= {MIN_READS} unique junction reads: {len(J):,}.\n')
    fh.write(f'Junction coordinate offset chosen from the data: +{OFF}.\n\n')
    fh.write('## 1. Agreement of the two quantifications\n\n')
    fh.write(f'Pearson r = {r_p.statistic:.3f}, Spearman rho = {r_s.statistic:.3f} '
             f'over {len(ok):,} library-event pairs scored by both.\n')
    fh.write(f'Prespecified reading: **{reading1}**.\n\n')
    fh.write('## 2. Does the Gate 4 result reproduce on junction PSI?\n\n')
    fh.write(f'Progression (control = 1 by construction): '
             f'control {scaled["control"]:.2f}, A485 {scaled["A485"]:.2f}, A485+DUX {scaled["A485+Dux"]:.2f} '
             f'over the {len(ev)} control-activated SE events that are scored in every library and whose junction '
             f'control change is at least {DPSI_MIN:.2f} (of {len(cand)} scored in every library).\n')
    dep = J.groupby('run').event.nunique()
    fh.write(f'Junction depth differs by arm ({dep.min():,} to {dep.max():,} events scored per library), so the '
             f'score is computed on complete cases only; scoring each library on whatever it covers would average '
             f'different event sets and is reported here as the reason for that restriction.\n')
    fh.write(f'I = {I:+.2f}, R = {R:+.2f}. Prespecified reading: **{reading2}**.\n')
    fh.write('Compare with the cross-fitted transcript-level values (I = -0.27, R = +0.22); the in-sample '
             'transcript-level magnitude is not the comparator.\n\n')
    fh.write('## 3. Events selected by rule for the read-level illustration\n\n')
    for row in top3.itertuples():
        fh.write(f'- `{row.event}` dPSI_control = {row.dPSI_control:+.2f}\n')

print(open(f'{B}/results/posthoc_junction_psi.md').read())
