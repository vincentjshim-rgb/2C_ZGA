"""Source data for the clock-matrix figure: what the clock value is a weighted sum of, across mouse stages.

The published trajectory is one number per stage. This writes the matrix behind it — the preprocessed feature of
every clock gene at every mouse stage of GSE225056, with the gene's clock coefficient — so the figure can show the
rows that carry the within-two-cell decrease instead of only their sum. Descriptive; no verdict depends on it.

The script re-runs the preprocessing of src/01_gate1_cross_species.py unchanged and asserts that the coefficient-
weighted column sums reproduce the stage means already stored in results/gate1_embryo_tage.tsv.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
sys.path.insert(0, f'{B}/src')
import tage_py as tp

MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
ORDER = ['Oocyte', 'Zygote', 'Early-2-cell', 'Late-2-cell', '4-cell', '8-cell', '16-cell']

C = pd.read_csv(f'{B}/data/GSE225056_Mouse_count_5prime.txt.gz', sep='\t', index_col=0)
C.index = C.index.str.strip('"')
C.columns = C.columns.str.strip('"')
C = C[C.index.str.startswith('ENS')]
C.index = C.index.str.replace(r'\.\d+$', '', regex=True)
C = C.groupby(level=0).sum()
stage = pd.Series([re.sub(r'_\d+$', '', c) for c in C.columns], index=C.columns)

# the same quality-control rule as Gate 1
det = np.log10((C >= 1).sum(axis=0))
mad = 1.4826 * np.median(np.abs(det - det.median()))
keep = det >= det.median() - 3 * mad
C, stage = C.loc[:, keep], stage[keep]
print(f'{C.shape[1]} libraries pass quality control, stages: {sorted(set(stage))}')

pp = tp.preprocess(C, reference_samples=list(stage.index[stage == 'Oocyte']))
X = pp['scaled_diff']
model, feats = tp.load_clock(MODEL)
sel = model.named_steps['featureselection'].get_support()
coef = pd.Series(model.named_steps['estimator'].coef_, index=np.array(feats)[sel])
coef = coef[coef != 0]
X = X.reindex(coef.index).fillna(0)
print(f'{len(coef):,} clock genes with a non-zero coefficient')

# stage means of the preprocessed feature, and the clock value they add up to
M = pd.DataFrame({s: X.loc[:, stage.index[stage == s]].mean(axis=1) for s in ORDER if (stage == s).any()})
got = (M.T * coef).T.sum()
ref = pd.read_csv(f'{B}/results/gate1_embryo_tage.tsv', sep='\t')
ref = (ref[ref.species == 'Mouse']
       .groupby('stage')['EN_Chronoage_Multispecies_Multitissue_scaleddiff'].mean().reindex(M.columns))
# The clock pipeline centres the features before the estimator (StandardScaler with with_std=False), so a
# coefficient-weighted sum of the uncentred features differs from the stored clock value by one constant, the
# same at every stage. Differences between stages are therefore exact, and that is what this figure shows.
off = float((got - ref).mean())
err = float(((got - ref) - off).abs().max())
print(f'stage-to-stage differences vs the stored clock values: max |difference| = {err:.2e} '
      f'(constant offset {off:+.4f} from the centring and the intercept)')
assert err < 1e-9, 'the matrix does not reproduce the stored stage-to-stage changes'
M = M - off / len(coef) * 0            # keep the matrix as preprocessed features; the offset is reported only
got = got - off

# order rows by the contribution to the within-two-cell change, the interval the paper decomposes
contrib = coef * (M['Late-2-cell'] - M['Early-2-cell'])
out = M.copy()
out.insert(0, 'coefficient', coef)
out.insert(1, 'contribution_E2C_to_L2C', contrib)
out = out.loc[contrib.sort_values().index]
out.index.name = 'gene'
out.to_csv(f'{B}/figures/source_data/figK_clock_matrix_mouse.tsv', sep='\t')
print(f'wrote figK_clock_matrix_mouse.tsv: {out.shape[0]:,} genes x {len(M.columns)} stages')
print('clock value per stage:', {s: round(float(v), 3) for s, v in got.items()})
print('within-two-cell change:', round(float(got['Late-2-cell'] - got['Early-2-cell']), 4),
      '| downward genes', int((contrib < 0).sum()), 'sum', round(float(contrib[contrib < 0].sum()), 3),
      '| upward', int((contrib > 0).sum()), 'sum', round(float(contrib[contrib > 0].sum()), 3))
