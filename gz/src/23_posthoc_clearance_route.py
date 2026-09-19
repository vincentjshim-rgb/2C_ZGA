"""POST HOC (plan: plan/POSTHOC_clearance_route_frozen.md). Split the Gate 2a composition-only simulation (mouse
GSE225056) into the part carried by maternal clock genes and the part carried by all other clock genes. Descriptive."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp

B = f'{ANALYSIS}/gz'
g2 = open(f'{B}/src/05_gate2a_composition.py', encoding='utf-8').read().split('rows, sim_rows, set_rows')[0]
ns = {'__file__': f'{B}/src/05_gate2a_composition.py'}
exec(g2, ns)                                                    # load, qc, ORDER, gene_sets, DATA_ZGA, MODEL, gm
load, qc, ORDER, gene_sets, DATA_ZGA, MODEL, gm = (ns[k] for k in ['load', 'qc', 'ORDER', 'gene_sets', 'DATA_ZGA', 'MODEL', 'gm'])

sp = 'Mouse'
C, stage, _ = load(sp)
keep, _ = qc(C)
C, stage = C.loc[:, keep], stage[keep]
mat, zyg, cpm = gene_sets(C, stage, ORDER[sp])
a, c = DATA_ZGA[sp].split('->')
share = cpm.loc[list(mat)].sum(axis=0) / 1e6
r = share[stage == c].mean() / share[stage == a].mean()
ref = list(stage.index[stage == 'Oocyte'])
oo = C[ref]
sim = oo.copy()
sim.loc[sim.index.isin(mat)] = sim.loc[sim.index.isin(mat)] * r
sim.columns = ['sim_' + x for x in sim.columns]
pp = tp.preprocess(pd.concat([oo, sim], axis=1), reference_samples=ref)
X = pp['scaled_diff']
model, feats = tp.load_clock(MODEL)
sel = model.named_steps['featureselection'].get_support()
coef = pd.Series(model.named_steps['estimator'].coef_, index=np.array(feats)[sel])
coef = coef[coef != 0]
dx = (X.reindex(coef.index)[sim.columns].mean(axis=1) - X.reindex(coef.index)[ref].mean(axis=1)).fillna(0)
contrib = coef * dx
delta_sim = float(contrib.sum())
stored = pd.read_csv(f'{B}/results/gate2a_simulation.tsv', sep='\t').set_index('species').loc[sp, 'delta_sim']
assert abs(delta_sim - stored) < 1e-9, (delta_sim, stored)

# maternal set in Entrez space
ens2ent = dict(zip(gm.Ensembl, gm.Entrez.astype('int64').astype(str)))
mat_entrez = {ens2ent[g] for g in mat if g in ens2ent}
is_mat = pd.Series([i in mat_entrez for i in contrib.index], index=contrib.index)
via_mat, via_other = float(contrib[is_mat].sum()), float(contrib[~is_mat].sum())
rows = [('delta_sim', delta_sim), ('clearance_ratio', float(r)), ('n_clock_genes', float(len(coef))),
        ('n_clock_genes_maternal', float(is_mat.sum())), ('n_clock_genes_with_nonzero_change', float((dx != 0).sum())),
        ('sum_via_maternal_clock_genes', via_mat), ('sum_via_other_clock_genes', via_other),
        ('fraction_via_other', via_other / delta_sim)]
R = pd.DataFrame(rows, columns=['quantity', 'value'])
R.to_csv(f'{B}/results/posthoc_clearance_route.tsv', sep='\t', index=False)
print(R.round(4).to_string(index=False))
