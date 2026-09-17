"""POST HOC (plan: plan/POSTHOC_rescue_frozen.md). Did the A485+Dux arm of GSE280522 restore zygotic transcription,
and did it restore the clock contributions identified in src/13_posthoc_gate3_contribution.py?
Descriptive; bootstrap CIs for arm differences only. Loading and QC mirror src/11_gate3_analysis.py."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import re, sys
import numpy as np, pandas as pd
B = f'{ANALYSIS}/gz'
sys.path.insert(0, f'{B}/src')
import tage_py as tp

rng = np.random.default_rng(20260918)
MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
SEL = SEL[SEL.gse == 'GSE280522']
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
t2g = dict(zip(t2g.tx, t2g.gene_id))
GT = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv')
GT = GT[GT.Entrez.notna()].drop_duplicates('Ensembl')
ENS2SYM = dict(zip(GT.Ensembl, GT['Gene.Symbol']))
SYM = dict(zip(GT.Entrez.astype('int64').astype(str), GT['Gene.Symbol']))
MARKER = re.compile(r'^(Zscan4|Tcstv|Tdpoz|Dux)', re.I)


def annotate(t):
    m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep', t)
    return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)]


def load():
    cols, meta = {}, []
    for _, r in SEL.iterrows():
        q = f'{B}/data/perturb/GSE280522/quant/{r.run_accession}'
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        p_al = float(re.search(r'"p_pseudoaligned": ([\d.]+)', open(f'{q}/run_info.json').read()).group(1))
        cols[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(t2g).values).sum()
        stage, arm = annotate(r.title_clean)
        meta.append(dict(run=r.run_accession, stage=stage, arm=arm, p_pseudoaligned=p_al))
    C = pd.DataFrame(cols).fillna(0)
    M = pd.DataFrame(meta).set_index('run')
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    M['qc_pass'] = (M.p_pseudoaligned >= 30) & (det >= det.median() - 3 * mad)
    M = M[M.qc_pass]
    return C[M.index], M


def boot_diff(a, b):
    d = np.array([rng.choice(a, len(a)).mean() - rng.choice(b, len(b)).mean() for _ in range(2000)])
    return np.percentile(d, 2.5), np.percentile(d, 97.5)


C, M = load()
cpm = C / C.sum(axis=0) * 1e6
lg = np.log2(cpm + 1)
e = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
l = M.index[(M.arm == 'control') & (M.stage == 'L2C')]
zyg = ((cpm[e] < 1).mean(axis=1) >= 0.8) & ((cpm[l] >= 10).mean(axis=1) >= 0.5)
mat = ((cpm[e] >= 10).mean(axis=1) >= 0.8) & (cpm[l].median(axis=1) <= cpm[e].median(axis=1) / 8)
print(f'zygotic set n = {int(zyg.sum())}, maternal set n = {int(mat.sum())} (control arm, Ensembl space)')

rows = []
for name, gset in [('zygotic', zyg), ('maternal', mat)]:
    sc = lg.loc[gset.values].mean(axis=0)
    per = {arm: sc[M.index[(M.arm == arm) & (M.stage == 'L2C')]].values for arm in ['control', 'A485', 'A485+Dux']}
    e_ctrl = sc[e].values
    for arm, v in per.items():
        rows.append(dict(set=name, arm=f'{arm} L2C', n=len(v), mean=v.mean()))
    rows.append(dict(set=name, arm='control E2C', n=len(e_ctrl), mean=e_ctrl.mean()))
    for a, b in [('A485', 'control'), ('A485+Dux', 'control'), ('A485+Dux', 'A485')]:
        lo, hi = boot_diff(per[a], per[b])
        rows.append(dict(set=name, arm=f'{a} - {b} (L2C)', n=np.nan, mean=per[a].mean() - per[b].mean(),
                         ci_lo=lo, ci_hi=hi))
Z = pd.DataFrame(rows)

sym = pd.Series({g: ENS2SYM.get(g, '') for g in lg.index})
mk = lg.loc[[bool(MARKER.match(str(s))) for s in sym.values]]
mk = mk.loc[mk.max(axis=1) > 0]
mk.index = [f'{g} {ENS2SYM.get(g, "")}' for g in mk.index]
MK = pd.DataFrame({arm: mk[M.index[(M.arm == arm) & (M.stage == 'L2C')]].mean(axis=1)
                   for arm in ['control', 'A485', 'A485+Dux']})
MK['control_E2C'] = mk[e].mean(axis=1)
MK = MK.sort_values('control', ascending=False)

pp = tp.preprocess(C, reference_samples=list(e))
X = pp['scaled_diff']
model, feats = tp.load_clock(MODEL)
sel = model.named_steps['featureselection'].get_support()
coef = pd.Series(model.named_steps['estimator'].coef_, index=np.array(feats)[sel])
coef = coef[coef != 0]
con = {}
for arm in ['control', 'A485', 'A485+Dux']:
    ee = M.index[(M.arm == arm) & (M.stage == 'E2C')]
    ll = M.index[(M.arm == arm) & (M.stage == 'L2C')]
    con[arm] = coef * (X.reindex(coef.index)[ll].mean(axis=1) - X.reindex(coef.index)[ee].mean(axis=1)).fillna(0)
K = pd.DataFrame(con)
K['symbol'] = [SYM.get(i, '') for i in K.index]
top20 = K.nsmallest(20, 'control')

Z.to_csv(f'{B}/results/posthoc_rescue_programme_scores.tsv', sep='\t', index=False)
MK.to_csv(f'{B}/results/posthoc_rescue_markers.tsv', sep='\t')
K.to_csv(f'{B}/results/posthoc_rescue_contributions.tsv', sep='\t')

pd.set_option('display.width', 250)
print('\n=== programme scores, mean log2(CPM+1) (POST HOC, descriptive) ===')
print(Z.round(4).to_string(index=False))
print('\n=== 2-cell / Dux-target markers, mean log2(CPM+1) at L2C ===')
print(MK.round(3).head(25).to_string())
print('\n=== clock contributions of the 20 largest control down-contributors ===')
print(top20[['symbol', 'control', 'A485', 'A485+Dux']].round(4).to_string())
print('\nsum over those 20 genes:', top20[['control', 'A485', 'A485+Dux']].sum().round(4).to_dict())
print('total D per arm:', K[['control', 'A485', 'A485+Dux']].sum().round(4).to_dict())
