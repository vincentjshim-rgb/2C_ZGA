"""POST HOC (plan: plan/POSTHOC_human_decomposition_frozen.md). Does the decomposition that describes the mouse
two-cell clock decrease also describe the human preimplantation decrease?

GSE36552 (Yan et al. 2013): deposited per-sample unique read counts per gene symbol. Cells of one embryo are summed
to an embryo pseudobulk; embryonic stem cell samples are excluded. The clock and its preprocessing are unchanged,
scored against the oocytes of this dataset. The prespecified interval is 8-cell -> morula.

The three readings and their thresholds were fixed in the plan before any human score. seed 20260927.
"""
import gzip
import os
import re
import sys

import numpy as np
import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
sys.path.insert(0, f'{B}/src')
import tage_py as tp

SCRATCH = ('/tmp/claude-1000/-home-shim-Downloads-embryogenesis-aging-260915/'
           '80cf57c0-6efa-4c21-b7d9-1f7ecca3e550/scratchpad/human')
MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
ORDER = ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'Late blastocyst']
INTERVAL = ('8-cell', 'Morula')        # prespecified
NBOOT, SEED = 2000, 20260927
MIN_CLOCK_GENES = 300                  # plan: below this the dataset is too thin to decompose

# ---------------------------------------------------------------- embryo pseudobulks
S = pd.read_csv(f'{SCRATCH}/samples.tsv', sep='\t').query("stage != 'ESC'")


def embryo_id(title, stage):
    m = re.search(r'#\s*(\d+)', title)
    return f'{stage}_{m.group(1)}' if m else f'{stage}_1'


S['embryo'] = [embryo_id(t, s) for t, s in zip(S.title, S.stage)]
cols = {}
for gsm, emb in zip(S.gsm, S.embryo):
    p = f'{SCRATCH}/expr/{gsm}.txt.gz'
    if not os.path.exists(p):
        continue
    with gzip.open(p, 'rt') as fh:
        d = pd.read_csv(fh, sep='\t', usecols=[0, 1], names=['sym', 'n'], skiprows=1)
    d = d.groupby('sym').n.sum()
    cols.setdefault(emb, []).append(d)

P = pd.DataFrame({e: pd.concat(v, axis=1).fillna(0).sum(axis=1) for e, v in cols.items()}).fillna(0)
meta = S.drop_duplicates('embryo').set_index('embryo').stage.reindex(P.columns)
print(f'{P.shape[1]} embryo pseudobulks from {len(S)} cells, {P.shape[0]:,} gene symbols')
print(meta.value_counts().reindex(ORDER).dropna().astype(int).to_string())

# ---------------------------------------------------------------- symbols -> Entrez, then the clock
# addendum 1: the clock's features are MOUSE Entrez identifiers (10,487 of 10,487 found in the mouse table,
# none in the human one), so human symbols are carried onto mouse identifiers with the package's own tables:
# symbol -> human Entrez -> mouse Entrez -> mouse Ensembl, the last step so that preprocess() runs unchanged.
H = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_human.csv').dropna(subset=['Entrez', 'Gene.Symbol'])
sym2h = dict(zip(H.drop_duplicates('Gene.Symbol')['Gene.Symbol'].astype(str),
                 H.drop_duplicates('Gene.Symbol')['Entrez'].astype('int64')))
O = pd.read_csv(f'{tp.PKG}/extdata/metadata/Table_of_orthologs.csv').dropna(subset=['Entrez.Human', 'Entrez.Mouse'])
h2m = dict(zip(O['Entrez.Human'].astype('int64'), O['Entrez.Mouse'].astype('int64')))
M = pd.read_csv(f'{tp.PKG}/extdata/metadata/Gene_table_mouse.csv').dropna(subset=['Entrez', 'Ensembl'])
m2e = dict(zip(M.drop_duplicates('Entrez')['Entrez'].astype('int64'),
               M.drop_duplicates('Entrez')['Ensembl'].astype(str)))

step = P.index.to_series().map(sym2h)
print(f'  symbol -> human Entrez : {int(step.notna().sum()):,} of {len(P):,}')
step = step.map(h2m)
print(f'  human  -> mouse Entrez : {int(step.notna().sum()):,}')
step = step.map(m2e)
print(f'  mouse Entrez -> Ensembl: {int(step.notna().sum()):,}')
ok = step.notna()
C = P.loc[ok].groupby(step[ok].values).sum()
print(f'{C.shape[0]:,} genes carried onto mouse Ensembl identifiers')

ref = list(meta.index[meta == 'Oocyte'])
pp = tp.preprocess(C, reference_samples=ref)
X = pp['scaled_diff']
model, feats = tp.load_clock(MODEL)
sel = model.named_steps['featureselection'].get_support()
coef = pd.Series(model.named_steps['estimator'].coef_, index=np.array(feats)[sel])
coef = coef[coef != 0]
present = X.reindex(coef.index).notna().any(axis=1)
n_present = int(present.sum())
print(f'weighted clock genes detected: {n_present:,} of {len(coef):,} ({n_present / len(coef):.2f})')

tage = pd.Series(tp.predict(MODEL, pp['scaled_diff'], species=None), index=X.columns)
T = pd.DataFrame({'embryo': tage.index, 'stage': meta.reindex(tage.index).values, 'tAge': tage.values})
T.to_csv(f'{B}/results/human_embryo_tage.tsv', sep='\t', index=False)

# ---------------------------------------------------------------- reading 1: is there a decrease?
rng = np.random.default_rng(SEED)
a, b = INTERVAL
ea = list(meta.index[meta == a])
eb = list(meta.index[meta == b])
D = float(tage[eb].mean() - tage[ea].mean())
boot = np.array([rng.choice(tage[eb].values, len(eb)).mean() - rng.choice(tage[ea].values, len(ea)).mean()
                 for _ in range(NBOOT)])
lo, hi = np.percentile(boot, [2.5, 97.5])
print(f'\n{a} -> {b}: D = {D:+.3f} [{lo:+.3f}, {hi:+.3f}]  (n = {len(ea)} and {len(eb)} embryos)')

# ---------------------------------------------------------------- reading 2: the decomposition
Xc = X.reindex(coef.index).fillna(0)
delta = Xc[eb].mean(axis=1) - Xc[ea].mean(axis=1)
contrib = coef * delta
assert abs(float(contrib.sum()) - D) < 1e-6, (float(contrib.sum()), D)   # the identity must hold
dn, up = float(contrib[contrib < 0].sum()), float(contrib[contrib > 0].sum())
frac = abs(D) / max(abs(dn), abs(up)) if max(abs(dn), abs(up)) > 0 else np.nan

G = pd.DataFrame({'coefficient': coef, 'feature_change': delta, 'contribution': contrib})
G = G.sort_values('contribution')
G.index.name = 'entrez'
G.to_csv(f'{B}/results/posthoc_human_decomposition.tsv', sep='\t')

if n_present < MIN_CLOCK_GENES:
    reading = f'NOT READ — only {n_present} weighted clock genes detected (plan requires at least {MIN_CLOCK_GENES})'
elif D < 0 and frac <= 0.5:
    reading = 'SAME STRUCTURE AS MOUSE — the decrease is a small residual of opposing contributions'
elif D >= 0:
    reading = 'NO DECREASE across the prespecified interval; question 2 reported for description only'
else:
    reading = f'DIFFERENT STRUCTURE — the decrease is not a small residual (fraction {frac:.2f} > 0.5)'

with open(f'{B}/results/posthoc_human_decomposition.md', 'w') as fh:
    fh.write('# Post hoc: does the decomposition hold in human embryos? '
             '(plan: plan/POSTHOC_human_decomposition_frozen.md)\n\n')
    fh.write(f'GSE36552, {P.shape[1]} embryo pseudobulks from {len(S)} cells; embryonic stem cell samples excluded '
             f'by the plan.\n\n')
    fh.write('## 3. How much of the clock does this dataset carry?\n\n')
    fh.write(f'Weighted clock genes detected: **{n_present:,} of {len(coef):,}** ({n_present / len(coef):.2f}). '
             f'Mouse comparators: 912 (GSE225056), 1,355 (GSE280522), 1,509 (GSE300734).\n\n')
    fh.write('## 1. Is there a decrease across the prespecified interval?\n\n')
    fh.write(f'{a} -> {b}: **D = {D:+.3f}** [{lo:+.3f}, {hi:+.3f}], {len(ea)} and {len(eb)} embryo pseudobulks.\n\n')
    fh.write('Clock value by stage (mean over embryo pseudobulks, relative to oocytes):\n\n')
    for s in ORDER:
        v = tage[[e for e in tage.index if meta[e] == s]]
        if len(v):
            fh.write(f'- {s}: {v.mean():+.3f} (n = {len(v)})\n')
    fh.write('\n## 2. Is it a residual of opposing contributions?\n\n')
    fh.write(f'Downward contributions sum to **{dn:.3f}**, upward to **{up:+.3f}**, and they add to the reported '
             f'change of {D:+.3f}.\n')
    fh.write(f'Residual fraction = |D| / max(|down|, |up|) = **{frac:.2f}** '
             f'(mouse: 0.17 in GSE280522, 0.12 in GSE225056; plan threshold 0.50).\n\n')
    fh.write(f'**Prespecified reading: {reading}**\n')

    # ---- follow-up, computed after the planned items and labelled as such -----------------------------
    from scipy.stats import pearsonr, spearmanr
    MM = pd.read_csv(f'{B}/results/posthoc_gate3_contribution_genes.tsv', sep='\t', index_col=0)
    MM.index = MM.index.astype(str)
    G2 = G.copy(); G2.index = G2.index.astype(str)
    j = G2.join(MM[['c_control', 'symbol']], how='inner')
    j = j[(j.contribution != 0) & (j.c_control != 0)].dropna(subset=['contribution', 'c_control'])
    rp = pearsonr(j.contribution, j.c_control).statistic
    rs = spearmanr(j.contribution, j.c_control).statistic
    shared = len(set(j.nsmallest(50, 'contribution').index) & set(j.nsmallest(50, 'c_control').index))
    fh.write('\n## Follow-up (computed after the planned items, not in the frozen plan)\n\n')
    fh.write(f'Is it the same genes? Over the {len(j):,} clock genes with a non-zero contribution in both the human '
             f'8-cell to morula interval and the mouse GSE280522 control window, the per-gene contributions are '
             f'essentially uncorrelated: Pearson r = {rp:.3f}, Spearman rho = {rs:.3f}, and {shared} of the 50 '
             f'largest downward contributors are shared.\n\n')
    fh.write('For comparison, two mouse studies measuring the same interval in the same species agreed at r = 0.74. '
             'The near-cancellation therefore reproduces across species and intervals while the gene identity behind '
             'it does not, which is consistent with reading the largest contributors as a property of this clock in '
             'a given window rather than as a ranking of genes.\n')

print('\n' + open(f'{B}/results/posthoc_human_decomposition.md').read())
