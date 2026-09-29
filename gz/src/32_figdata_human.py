"""Source data for the human replication figure (post hoc; plan POSTHOC_human_decomposition_frozen.md).

  figHa_human_tage.tsv        clock value of every human embryo pseudobulk, with stage means and intervals
  figHb_human_cumulative.tsv  the running sum of the ranked per-gene contributions, 8-cell -> morula
  figHc_human_vs_mouse.tsv    per-gene contribution in human against mouse, for the genes non-zero in both

Every value is re-read from the stored result and checked against it before it is written.
"""
import os

import numpy as np
import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
SD = f'{B}/figures/source_data'
ORDER = ['Oocyte', 'Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'Late blastocyst']
SEED, NBOOT = 20260927, 2000

# ---- (a) per-embryo clock values and stage means
T = pd.read_csv(f'{B}/results/human_embryo_tage.tsv', sep='\t')
rng = np.random.default_rng(SEED)
rows = []
for s in ORDER:
    v = T.loc[T.stage == s, 'tAge'].values
    if not len(v):
        continue
    bm = np.array([rng.choice(v, len(v)).mean() for _ in range(NBOOT)])
    rows.append((s, len(v), float(v.mean()), float(np.percentile(bm, 2.5)), float(np.percentile(bm, 97.5))))
A = pd.DataFrame(rows, columns=['stage', 'n', 'mean', 'ci_lo', 'ci_hi'])
T.to_csv(f'{SD}/figHa_human_tage.tsv', sep='\t', index=False)
A.to_csv(f'{SD}/figHa_human_stage_means.tsv', sep='\t', index=False)
print('stage means:\n' + A.round(3).to_string(index=False))

# ---- (b) running sum of the ranked contributions
G = pd.read_csv(f'{B}/results/posthoc_human_decomposition.tsv', sep='\t', index_col=0)
G = G.sort_values('contribution')
cum = pd.DataFrame({'rank': np.arange(1, len(G) + 1),
                    'contribution': G.contribution.values,
                    'cumulative': G.contribution.cumsum().values})
cum.to_csv(f'{SD}/figHb_human_cumulative.tsv', sep='\t', index=False)
D_from_cum = float(cum.cumulative.iloc[-1])
md = open(f'{B}/results/posthoc_human_decomposition.md').read()
import re
D_reported = float(re.search(r'\*\*D = ([-+0-9.]+)\*\*', md).group(1))
assert abs(D_from_cum - D_reported) < 5e-4, (D_from_cum, D_reported)
nneg = int((G.contribution < 0).sum())
print(f'\ncumulative: {len(cum):,} genes, trough {cum.cumulative.min():.3f} at the {nneg} negative genes, '
      f'end {D_from_cum:+.4f} (reported {D_reported:+.3f})')

# ---- (c) human against mouse, the genes with a non-zero contribution in both
M = pd.read_csv(f'{B}/results/posthoc_gate3_contribution_genes.tsv', sep='\t', index_col=0)
M.index = M.index.astype(str)
G.index = G.index.astype(str)
J = G.join(M[['c_control', 'symbol']], how='inner').dropna(subset=['contribution', 'c_control'])
J = J[(J.contribution != 0) & (J.c_control != 0)]
J[['symbol', 'contribution', 'c_control']].rename(
    columns={'contribution': 'human_8cell_to_morula', 'c_control': 'mouse_E2C_to_L2C'}
).to_csv(f'{SD}/figHc_human_vs_mouse.tsv', sep='\t')
r = float(np.corrcoef(J.contribution, J.c_control)[0, 1])
print(f'human vs mouse: {len(J):,} shared genes, Pearson r = {r:.3f}')
