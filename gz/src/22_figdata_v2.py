"""Source data for the redrawn Figures 3-5 (display transforms of stored results; no new analysis).
- Cumulative contribution curves: clock genes ranked by contribution, running sum (Figures 3a, 4c).
- PSI of the control-defined ZSA events in P1 across the 23 libraries (Figure 5a) and per-arm oriented dPSI (5b).
Each re-extraction is checked against the stored result before it is written. Output: figures/source_data/*.tsv"""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import numpy as np, pandas as pd

B = f'{ANALYSIS}/gz'
OUT = f'{B}/figures/source_data'
R = f'{B}/results'
rd = lambda f, **kw: pd.read_csv(f'{R}/{f}.tsv', sep='\t', **kw)
A = rd('gate3_arm_drops')
drop = lambda gse, arm: float(A[(A.gse == gse) & (A.variant == 'V0') & (A.arm == arm)]['drop'].iloc[0])

# ---- cumulative contribution curves ---------------------------------------------------------------------
RC = rd('posthoc_rescue_contributions', index_col=0)            # P1: control, A485, A485+Dux per clock gene
P3 = rd('posthoc_gate3_contribution_genes_P3', index_col=0)
rows = []
for tag, gse, series in [('P1_control', 'GSE280522', RC['control']), ('P1_A485', 'GSE280522', RC['A485']),
                         ('P1_A485+Dux', 'GSE280522', RC['A485+Dux']), ('P3_control', 'GSE300734', P3['c_control']),
                         ('P3_Brg1_matKO', 'GSE300734', P3['c_perturbed'])]:
    arm = tag.split('_', 1)[1]
    v = series.sort_values()                                     # most negative first
    cs = v.cumsum()
    assert abs(cs.iloc[-1] - drop(gse, arm)) < 1e-9, (tag, cs.iloc[-1], drop(gse, arm))
    assert abs(v[v < 0].sum() - cs.min()) < 1e-12
    for k, (g, c) in enumerate(cs.items(), start=1):
        rows.append(dict(series=tag, rank=k, gene=g, contribution=v[g], cumulative=c))
# Figure 4c: the perturbed arms accumulated in the CONTROL gene order, so the curves compare the same genes
order = RC['control'].sort_values().index
n_neg = int((RC['control'] < 0).sum())
for arm in ['A485', 'A485+Dux']:
    v = RC[arm].reindex(order)
    cs = v.cumsum()
    assert abs(cs.iloc[-1] - drop('GSE280522', arm)) < 1e-9
    for k, (g, c) in enumerate(cs.items(), start=1):
        rows.append(dict(series=f'P1_{arm}_ctrlorder', rank=k, gene=g, contribution=v[g], cumulative=c))
CC = pd.DataFrame(rows)
at = CC[CC['rank'] == n_neg].set_index('series').cumulative
print(f'control-defined downward genes n = {n_neg}; running sum at that rank:', at.round(3).to_dict())
assert abs(at['P1_A485+Dux_ctrlorder'] / at['P1_control'] - 0.85) < 0.01       # the 85% recovery quoted in the text
CC.to_csv(f'{OUT}/fig3a_4c_cumulative_contributions.tsv', sep='\t', index=False)
print(CC.groupby('series').cumulative.agg(['min', 'last']).round(3).to_string())

# ---- PSI of control ZSA events, P1 ----------------------------------------------------------------------
def read_expr(path):
    with open(path) as fh:
        c = fh.readline().rstrip('\n').split('\t')
    df = pd.read_csv(path, sep='\t', skiprows=1, header=None, index_col=0)
    df.columns = c
    return df


G = rd('gate4_progression')
g = G[G.gse == 'GSE280522'].set_index('run')
cols = pd.read_csv(f'{OUT}/fig4a_heatmap_columns.tsv', sep='\t')          # the 23 Gate 3 QC-passing P1 libraries
assert set(cols.run) == set(g.index), 'Gate 3 and Gate 4 library sets differ'
ev = rd('gate4_events_filtered_GSE280522', index_col=0)
Z = ev[ev.zsa].copy()
P = read_expr(f'{B}/data/gate4/GSE280522.psi').loc[Z.index, list(g.index)]
assert P.isna().sum().sum() == 0
# check: progression score recomputed from this matrix equals the stored Gate 4 values
ce = g.index[(g.arm == 'control') & (g.stage == 'E2C')]
cl = g.index[(g.arm == 'control') & (g.stage == 'L2C')]
base = P[ce].mean(axis=1)
d = P[cl].mean(axis=1) - base
assert np.allclose(d, Z.dPSI_control)
score = (P.sub(base, axis=0).mul(np.sign(d), axis=0).div(d.abs(), axis=0)).mean(axis=0)
assert np.allclose(score.reindex(g.index), g.progression, atol=1e-9)
Z = Z.sort_values('dPSI_control', ascending=False)
M = P.loc[Z.index]
M.insert(0, 'event_type', Z.event_type)
M.insert(1, 'dPSI_control', Z.dPSI_control)
M.to_csv(f'{OUT}/fig5a_psi_zsa_events_P1.tsv', sep='\t')

# per-arm dPSI (late - early), oriented by the sign of the control change
rows = []
for arm in ['control', 'A485', 'A485+Dux']:
    e = g.index[(g.arm == arm) & (g.stage == 'E2C')]
    l = g.index[(g.arm == arm) & (g.stage == 'L2C')]
    dd = (P[l].mean(axis=1) - P[e].mean(axis=1)) * np.sign(Z.dPSI_control.reindex(P.index))
    rows.append(pd.DataFrame({'event': P.index, 'event_type': Z.event_type.reindex(P.index).values, 'arm': arm,
                              'dPSI_oriented': dd.values}))
D = pd.concat(rows, ignore_index=True)
D.to_csv(f'{OUT}/fig5b_dpsi_per_arm_P1.tsv', sep='\t', index=False)
print(D.groupby('arm').dPSI_oriented.describe()[['mean', '50%']].round(3).to_string())
print('written to', OUT)
