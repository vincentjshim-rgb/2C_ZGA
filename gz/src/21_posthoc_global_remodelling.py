"""POST HOC (plan: plan/POSTHOC_global_remodelling_frozen.md). Is the within-2-cell clock decrease explained by a
genome-wide shift of the preprocessed features, or by gene-specific changes? Exact split D = G + S with
G = (sum of coefficients) x (mean feature change), plus a coefficient-reassignment null and direction-restricted
partial drops. Descriptive; no gate verdict changes. Loading mirrors src/13_posthoc_gate3_contribution.py. seed 20260920."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys
import numpy as np, pandas as pd
B = f'{ANALYSIS}/gz'
sys.path.insert(0, f'{B}/src')
import importlib.util
spec = importlib.util.spec_from_file_location('c13', f'{B}/src/13_posthoc_gate3_contribution.py')

import re
import tage_py as tp

MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
t2g = dict(zip(t2g.tx, t2g.gene_id))
rng = np.random.default_rng(20260920)
NPERM = 2000


def annotate(gse, t):
    if gse == 'GSE280522':
        m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep', t)
        return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)]
    m = re.search(r'(WT|KO)_(E2C|L2C)_RNA', t)                      # GSE300734
    return m.group(2), {'WT': 'control', 'KO': 'Brg1_matKO'}[m.group(1)]


def load_dataset(gse):
    S = SEL[SEL.gse == gse]
    cols, meta = {}, []
    for _, r in S.iterrows():
        q = f'{B}/data/perturb/{gse}/quant/{r.run_accession}'
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        p_al = float(re.search(r'"p_pseudoaligned": ([\d.]+)', open(f'{q}/run_info.json').read()).group(1))
        cols[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(t2g).values).sum()
        stage, arm = annotate(gse, r.title_clean)
        meta.append(dict(run=r.run_accession, stage=stage, arm=arm, p_pseudoaligned=p_al))
    C = pd.DataFrame(cols).fillna(0)
    M = pd.DataFrame(meta).set_index('run')
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    M['qc_pass'] = (M.p_pseudoaligned >= 30) & (det >= det.median() - 3 * mad)
    M = M[M.qc_pass]
    return C[M.index], M


def feature_changes(gse):
    """Per-arm change of the preprocessed clock features between late and early two-cell libraries."""
    C, M = load_dataset(gse)
    pp = tp.preprocess(C, reference_samples=list(M.index[(M.arm == 'control') & (M.stage == 'E2C')]))
    X = pp['scaled_diff']
    model, feats = tp.load_clock(MODEL)
    sel = model.named_steps['featureselection'].get_support()
    coef = pd.Series(model.named_steps['estimator'].coef_, index=np.array(feats)[sel])
    coef = coef[coef != 0]
    d = {}
    for arm in M.arm.unique():
        e = M.index[(M.arm == arm) & (M.stage == 'E2C')]
        l = M.index[(M.arm == arm) & (M.stage == 'L2C')]
        d[arm] = (X.reindex(coef.index)[l].mean(axis=1) - X.reindex(coef.index)[e].mean(axis=1)).fillna(0)
    # expression direction across the control window, for the Figure 3a categories
    cpm = tp.map_mouse_ensembl(tp.filter_genes(C))
    cpm = cpm / cpm.sum(axis=0) * 1e6
    ce = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
    cl = M.index[(M.arm == 'control') & (M.stage == 'L2C')]
    lfc = np.log2((cpm[cl].mean(axis=1) + 1) / (cpm[ce].mean(axis=1) + 1)).reindex(coef.index)
    return coef, pd.DataFrame(d), lfc


def split(coef, delta, lfc, dataset, arm):
    """Exact global/specific split, reassignment null and direction-restricted partial sums."""
    D = float((coef * delta).sum())
    G = float(coef.sum() * delta.mean())
    S = D - G
    perm = np.array([float((coef.values[rng.permutation(len(coef))] * delta.values).sum()) for _ in range(NPERM)])
    rows = [('D_observed', D), ('G_global_component', G), ('S_specific_component', S),
            ('abs_S_over_abs_G', abs(S) / abs(G) if G != 0 else np.nan),
            ('null_mean', float(perm.mean())), ('null_p2.5', float(np.percentile(perm, 2.5))),
            ('null_p97.5', float(np.percentile(perm, 97.5))),
            ('fraction_of_null_below_observed', float((perm < D).mean())),
            ('sum_of_coefficients', float(coef.sum())), ('mean_feature_change', float(delta.mean())),
            ('D_features_up_only', float((coef * delta.where(delta > 0, 0)).sum())),
            ('D_features_down_only', float((coef * delta.where(delta < 0, 0)).sum())),
            ('D_expression_rises', float((coef * delta.where(lfc > 0, 0)).sum())),
            ('D_expression_falls', float((coef * delta.where(lfc <= 0, 0)).sum())),
            ('n_genes', float(len(coef))), ('n_genes_feature_change_nonzero', float((delta != 0).sum()))]
    return pd.DataFrame([dict(dataset=dataset, arm=arm, quantity=q, value=v) for q, v in rows])


A = pd.read_csv(f'{B}/results/gate3_arm_drops.tsv', sep='\t')
out, verdict_bits = [], {}
for gse, tag, perturbed in [('GSE280522', 'P1_GSE280522', 'A485'), ('GSE300734', 'P3_GSE300734', 'Brg1_matKO')]:
    coef, DL, lfc = feature_changes(gse)
    for arm in DL.columns:
        want = A[(A.gse == gse) & (A.variant == 'V0') & (A.arm == arm)]['drop']
        r = split(coef, DL[arm], lfc, tag, arm)
        got = float(r[r.quantity == 'D_observed'].value.iloc[0])
        assert abs(got - float(want.iloc[0])) < 1e-9, (gse, arm, got, want.iloc[0])   # matches Gate 3
        out.append(r)
        if arm == 'control':
            verdict_bits[tag] = r
    # interaction: change of the feature change between arms
    ri = split(coef, DL[perturbed] - DL['control'], lfc, tag, f'interaction_{perturbed}_minus_control')
    out.append(ri)
    verdict_bits[tag + '_I'] = ri
R = pd.concat(out, ignore_index=True)
R.to_csv(f'{B}/results/posthoc_global_remodelling.tsv', sep='\t', index=False)

get = lambda t, q: float(t[t.quantity == q].value.iloc[0])
ok = {}
for tag in ['P1_GSE280522', 'P3_GSE300734']:
    t = verdict_bits[tag]
    ok[tag] = (abs(get(t, 'S_specific_component')) > abs(get(t, 'G_global_component'))
               and get(t, 'D_observed') < get(t, 'null_p2.5'))
verdict = 'NOT EXPLAINED' if all(ok.values()) else ('PARTLY EXPLAINED' if any(ok.values()) else 'EXPLAINED')
tI = verdict_bits['P1_GSE280522_I']
verdict_I = ('NOT EXPLAINED' if (abs(get(tI, 'S_specific_component')) > abs(get(tI, 'G_global_component'))
                                 and get(tI, 'D_observed') > get(tI, 'null_p97.5')) else 'NOT SUPPORTED')

with open(f'{B}/results/posthoc_global_remodelling.md', 'w') as fh:
    fh.write('# Post hoc: global transcriptome remodelling vs gene-specific change '
             '(plan: plan/POSTHOC_global_remodelling_frozen.md)\n\n')
    fh.write(f'VERDICT (control arms, prespecified rule): **{verdict}** by global remodelling\n\n')
    fh.write(f'P1 interaction (A485 - control): **{verdict_I}** by global remodelling\n\n')
    fh.write('D = G + S; G = (sum of coefficients) x (mean feature change) is the drop expected if every clock gene '
             'moved by the average amount; S is the gene-specific remainder. The null re-assigns the coefficient '
             'vector to genes 2,000 times (seed 20260920), holding the observed feature changes fixed.\n\n')
    piv = R.pivot_table(index=['dataset', 'arm'], columns='quantity', values='value')
    cols = ['D_observed', 'G_global_component', 'S_specific_component', 'abs_S_over_abs_G', 'null_mean',
            'null_p2.5', 'null_p97.5', 'fraction_of_null_below_observed', 'D_features_down_only',
            'D_features_up_only', 'D_expression_rises', 'D_expression_falls']
    fh.write(piv[cols].round(4).to_string() + '\n')
print('VERDICT:', verdict, '| interaction:', verdict_I)
pd.set_option('display.width', 250)
print(R.pivot_table(index=['dataset', 'arm'], columns='quantity', values='value')[
    ['D_observed', 'G_global_component', 'S_specific_component', 'abs_S_over_abs_G', 'null_mean', 'null_p2.5',
     'null_p97.5', 'fraction_of_null_below_observed', 'D_features_down_only', 'D_features_up_only']].round(4).to_string())
