"""Gate 2b, dataset R1 (plan: plan/GATE2B_PLAN_frozen.md). GSE45719 (Deng 2014) single-cell Smart-seq,
embryo pseudobulk. Does transcriptomic age drop within the 2-cell stage? seed 20260917."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys, glob, os, re
import numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp

B = f'{ANALYSIS}/gz'
MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
rng = np.random.default_rng(20260917)
PREFIX = [('zy', 'Zygote'), ('early2cell', 'Early-2C'), ('mid2cell', 'Mid-2C'), ('late2cell', 'Late-2C'),
          ('4cell', '4C'), ('8cell', '8C'), ('16cell', '16C'), ('earlyblast', 'EarlyBl'), ('midblast', 'MidBl'),
          ('lateblast', 'LateBl')]
ORDER = [s for _, s in PREFIX]
PRE_BLAST = ORDER[:7]
EXCL = re.compile(r'split|pooled|smartseq2|C57twocell|fibroblast|liver|BXC', re.I)

cells = {}
for f in sorted(glob.glob(f'{B}/data/GSE45719/*_expression.txt.gz')):
    name = re.sub(r'^GSM\d+_', '', os.path.basename(f)).replace('_expression.txt.gz', '')
    if EXCL.search(name):
        continue
    stage = next((s for p, s in PREFIX if name.startswith(p)), None)
    if stage is None:
        continue
    if stage == '8C' and not re.match(r'^8cell_(1|2|5|8)-', name):
        continue
    if stage == '16C' and not re.match(r'^16cell_(1|4|5|6)-', name):
        continue
    embryo = name if stage == 'Zygote' else name.rsplit('-', 1)[0]
    x = pd.read_csv(f, sep='\t', usecols=[0, 3])
    x.columns = ['symbol', 'reads']
    cells[name] = (stage, embryo, x.groupby('symbol').reads.sum())

meta = pd.DataFrame([(n, s, e) for n, (s, e, _) in cells.items()], columns=['cell', 'stage', 'embryo'])
X = pd.DataFrame({n: v for n, (_, _, v) in cells.items()}).fillna(0)
gm = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')
u = gm.drop_duplicates('Gene.Symbol', keep=False)
sym2ens = dict(zip(u['Gene.Symbol'], u.Ensembl))
ens = X.index.to_series().map(sym2ens)
X = X.loc[ens.notna()].groupby(ens[ens.notna()].values).sum()
P = X.T.groupby(meta.set_index('cell').embryo.reindex(X.columns).values).sum().T
emb = meta.groupby('embryo').agg(stage=('stage', 'first'), n_cells=('cell', 'size'))
emb = emb.loc[P.columns]
det = np.log10((P >= 1).sum(axis=0))
mad = 1.4826 * np.median(np.abs(det - det.median()))
keep = det >= det.median() - 3 * mad
P, emb = P.loc[:, keep], emb[keep]
stage = emb.stage
print('cells used %d; embryos %d (QC excluded %d); genes %d' % (len(meta), keep.sum(), (~keep).sum(), P.shape[0]))
print(emb.groupby('stage').agg(embryos=('n_cells', 'size'), cells=('n_cells', 'sum')).reindex(ORDER).to_string())


def score(C):
    pp = tp.preprocess(C, reference_samples=list(stage.index[stage == 'Zygote']))
    return tp.predict(MODEL, pp['scaled_diff'], species=None), tp.feature_coverage(MODEL, pp['scaled_diff'])


def analyse(v, label):
    groups = {s: v[stage == s].values for s in ORDER}
    means = {s: g.mean() for s, g in groups.items()}
    boot = {s: np.array([rng.choice(g, len(g)).mean() for _ in range(2000)]) for s, g in groups.items()}
    rows = []
    for a, c in zip(ORDER[:-1], ORDER[1:]):
        bd = boot[c] - boot[a]
        rows.append(dict(analysis=label, interval=f'{a}->{c}', delta=means[c] - means[a],
                         ci_lo=np.percentile(bd, 2.5), ci_hi=np.percentile(bd, 97.5)))
    I = pd.DataFrame(rows)
    pre = I[I.interval.isin([f'{a}->{c}' for a, c in zip(PRE_BLAST[:-1], PRE_BLAST[1:])])].index
    bmat = np.column_stack([boot[c] - boot[a] for a, c in zip(PRE_BLAST[:-1], PRE_BLAST[1:])])
    amin = bmat.argmin(axis=1)
    I['p_boot_most_negative_preblast'] = np.nan
    for j, k in enumerate(pre):
        I.loc[k, 'p_boot_most_negative_preblast'] = float(np.mean(amin == j))
    s1b = boot['Late-2C'] - boot['Early-2C']
    s1 = dict(analysis=label, S1=means['Late-2C'] - means['Early-2C'], S1_lo=np.percentile(s1b, 2.5), S1_hi=np.percentile(s1b, 97.5))
    return I, s1, means


v0, cov = score(P)
I0, s10, m0 = analyse(v0, 'primary')

# secondary: Gate 2a V2 analogue (earliest two stages play oocyte / zygote roles)
cpm = P / P.sum(axis=0) * 1e6
c1, c2 = list(stage.index[stage == 'Zygote']), list(stage.index[stage == 'Early-2C'])
later = [s for s in ORDER if s not in ('Zygote', 'Early-2C')]
med1 = cpm[c1].median(axis=1)
mat = ((cpm[c1] >= 10).mean(axis=1) >= 0.8) & ((cpm[c2] >= 10).mean(axis=1) >= 0.8) & \
      pd.concat([cpm[list(stage.index[stage == s])].median(axis=1) <= med1 / 8 for s in later], axis=1).any(axis=1)
zyg = ((cpm[c1] < 1).mean(axis=1) >= 0.8) & ((cpm[c2] < 1).mean(axis=1) >= 0.8) & \
      pd.concat([(cpm[list(stage.index[stage == s])] >= 10).mean(axis=1) >= 0.5 for s in later], axis=1).any(axis=1)
v2, _ = score(P.loc[~(mat | zyg)])
I2, s12, _ = analyse(v2, 'V2_drop_all_dynamic')

I = pd.concat([I0, I2])
I.to_csv(f'{B}/results/gate2b_R1_intervals.tsv', sep='\t', index=False)
pd.DataFrame({'stage': stage, 'n_cells': emb.n_cells, 'tAge_primary': v0, 'tAge_V2': v2}).to_csv(
    f'{B}/results/gate2b_R1_embryo_tage.tsv', sep='\t')

pre0 = I0.dropna(subset=['p_boot_most_negative_preblast'])
kmin = pre0.delta.idxmin()
p2c = pre0.loc[pre0.interval.isin(['Early-2C->Mid-2C', 'Mid-2C->Late-2C']), 'p_boot_most_negative_preblast'].sum()
sig = s10['S1_hi'] < 0
if s10['S1'] > 0 or not sig:
    verdict = 'NOT REPLICATED'
elif pre0.loc[kmin, 'interval'] in ('Early-2C->Mid-2C', 'Mid-2C->Late-2C') and p2c >= 0.6:
    verdict = 'REPLICATED'
else:
    verdict = 'PARTIAL'

pd.set_option('display.width', 250)
print('clock feature coverage %.3f; dynamic genes removed in V2: %d' % (cov, int((mat | zyg).sum())))
print(I.round(4).to_string(index=False))
print('S1 primary', {k: round(v, 4) for k, v in s10.items() if k != 'analysis'}, '| V2', {k: round(v, 4) for k, v in s12.items() if k != 'analysis'})
print('most negative pre-blastocyst interval:', pre0.loc[kmin, 'interval'], '| P_boot(within 2-cell) %.3f' % p2c)
print('R1 VERDICT:', verdict)
with open(f'{B}/results/gate2b_R1_decision.md', 'w', encoding='utf-8') as fh:
    fh.write('# Gate 2b R1 (GSE45719) decision\n\n')
    fh.write(f'VERDICT: **{verdict}**\n\n- S1 (late − early 2-cell), primary: {s10}\n- S1, V2 (dynamic genes removed): {s12}\n')
    fh.write(f'- most negative pre-blastocyst interval: {pre0.loc[kmin, "interval"]}; P_boot within 2-cell {p2c:.3f}\n')
