"""Gate 2b, dataset R2 (plan: plan/GATE2B_PLAN_frozen.md). GSE66582 (Wu 2016) bulk embryos, kallisto re-quant.
Direction-only readout (n = 2 per stage). Reference = MII oocytes."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp

B = f'{ANALYSIS}/gz'
MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
RUNS = {'SRR2927026': 'MII', 'SRR2927027': 'MII', 'SRR1840514': 'Zygote', 'SRR1840515': 'Zygote',
        'SRR2927028': 'Early-2C', 'SRR2927029': 'Early-2C', 'SRR1840516': '2C', 'SRR1840517': '2C',
        'SRR1840518': '4C', 'SRR1840519': '4C', 'SRR1840520': '8C', 'SRR1840521': '8C',
        'SRR1840522': 'ICM', 'SRR1840523': 'ICM', 'SRR1840526': 'ICM'}
ORDER = ['MII', 'Zygote', 'Early-2C', '2C', '4C', '8C']

t2g = pd.read_csv(f'{ANALYSIS}/metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
t2g = dict(zip(t2g.tx, t2g.gene_id))
cols = {}
for r, st in RUNS.items():
    a = pd.read_csv(f'{B}/data/GSE66582/quant/{r}/abundance.tsv', sep='\t', index_col=0)
    g = a.index.to_series().map(t2g)
    cols[r] = a.est_counts.groupby(g.values).sum()
C = pd.DataFrame(cols).fillna(0)
stage = pd.Series(RUNS)[C.columns]
det = np.log10((C >= 1).sum(axis=0))
print('genes %d; libraries %d; detected genes (log10):' % C.shape, det.round(2).to_dict())


def score(Cx):
    pp = tp.preprocess(Cx, reference_samples=list(stage.index[stage == 'MII']))
    return tp.predict(MODEL, pp['scaled_diff'], species=None), tp.feature_coverage(MODEL, pp['scaled_diff'])


def readout(v, label):
    m = {s: v[stage == s].mean() for s in ORDER + ['ICM']}
    rows = [dict(analysis=label, interval=f'{a}->{c}', delta=m[c] - m[a],
                 rep_values=';'.join('%.3f' % x for x in v[stage == c])) for a, c in zip(ORDER[:-1], ORDER[1:])]
    I = pd.DataFrame(rows)
    target = m['2C'] - m['Early-2C']
    consistent = bool(target < 0 and I.delta.idxmin() == I.index[I.interval == 'Early-2C->2C'][0])
    return I, m, consistent


v0, cov = score(C)
I0, m0, cons0 = readout(v0, 'primary')
cpm = C / C.sum(axis=0) * 1e6
c1, c2 = list(stage.index[stage == 'MII']), list(stage.index[stage == 'Zygote'])
later = ['Early-2C', '2C', '4C', '8C']
med1 = cpm[c1].median(axis=1)
mat = ((cpm[c1] >= 10).mean(axis=1) >= 0.8) & ((cpm[c2] >= 10).mean(axis=1) >= 0.8) & \
      pd.concat([cpm[list(stage.index[stage == s])].median(axis=1) <= med1 / 8 for s in later], axis=1).any(axis=1)
zyg = ((cpm[c1] < 1).mean(axis=1) >= 0.8) & ((cpm[c2] < 1).mean(axis=1) >= 0.8) & \
      pd.concat([(cpm[list(stage.index[stage == s])] >= 10).mean(axis=1) >= 0.5 for s in later], axis=1).any(axis=1)
v2, _ = score(C.loc[~(mat | zyg)])
I2, m2, cons2 = readout(v2, 'V2_drop_all_dynamic')

I = pd.concat([I0, I2])
I.to_csv(f'{B}/results/gate2b_R2_intervals.tsv', sep='\t', index=False)
pd.DataFrame({'stage': stage, 'tAge_primary': v0, 'tAge_V2': v2}).to_csv(f'{B}/results/gate2b_R2_library_tage.tsv', sep='\t')
verdict = 'CONSISTENT' if cons0 else 'INCONSISTENT'
pd.set_option('display.width', 200)
print('clock coverage %.3f; dynamic genes removed %d' % (cov, int((mat | zyg).sum())))
print(I.round(4).to_string(index=False))
print('stage means primary', {k: round(x, 3) for k, x in m0.items()})
print('R2 VERDICT (direction only):', verdict, '| V2:', 'CONSISTENT' if cons2 else 'INCONSISTENT')
with open(f'{B}/results/gate2b_R2_decision.md', 'w', encoding='utf-8') as fh:
    fh.write(f'# Gate 2b R2 (GSE66582) decision, direction only (n = 2 per stage)\n\nVERDICT: **{verdict}** (V2: {"CONSISTENT" if cons2 else "INCONSISTENT"})\n\n')
    fh.write(I.round(4).to_string(index=False) + '\n')
