"""Supplementary Tables S3-S5 (no new analysis). Re-shapes stored results into one workbook plus TSVs.
S3: clock value of every scored embryo or library (Gates 1, 2b, 3), with GEO/SRA identifiers where they exist.
S4: per-gene contributions to the within-two-cell drop (post-hoc decomposition and rescue analysis).
S5: per-fold values of the cross-fitted Gate 4 check (post hoc).
Each table is checked against the stored summary values before it is written.
Output: results/supp_tables/TableS3-S5.tsv and SupplementaryTables_S3-S5.xlsx"""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import re
import numpy as np, pandas as pd

B = f'{ANALYSIS}/gz'
R = f'{B}/results'
OUT = f'{R}/supp_tables'
os.makedirs(OUT, exist_ok=True)
rd = lambda f, **kw: pd.read_csv(f'{R}/{f}.tsv', sep='\t', **kw)
ACC = pd.read_csv(f'{B}/metadata/perturb/gate3_accessions_geo_sra.tsv', sep='\t').set_index('srr')

# ---- S3 ------------------------------------------------------------------------------------------------
rows = []
g1 = rd('gate1_embryo_tage', index_col=0)
s225 = pd.read_csv(f'{B}/metadata/GSE225056_samples.tsv', sep='\t').set_index('title').gsm
for emb, r in g1.iterrows():
    rows.append(dict(dataset='GSE225056', analysis='Gate 1' + (' (descriptive only)' if r.species == 'Rhesus' else ''),
                     species=r.species, sample_id=emb, gsm=s225.get(f'{r.species}_{emb}', ''), srr='',
                     stage=r.stage, arm='', variant='V0', reference='oocytes of the same species',
                     clock_value=r.EN_Chronoage_Multispecies_Multitissue_scaleddiff))
g2 = rd('gate2b_R1_embryo_tage', index_col=0)
for emb, r in g2.iterrows():
    for v, col in [('V0', 'tAge_primary'), ('V2', 'tAge_V2')]:
        rows.append(dict(dataset='GSE45719', analysis='Gate 2b R1', species='Mouse', sample_id=f'{emb} ({r.n_cells} cells)',
                         gsm='', srr='', stage=r.stage, arm='', variant=v, reference='zygotes', clock_value=r[col]))
ena = pd.read_csv(f'{B}/metadata/SRP055882_ena.tsv', sep='\t').set_index('run_accession')
g3 = rd('gate2b_R2_library_tage', index_col=0)
for run, r in g3.iterrows():
    gsm = re.search(r'GSM\d+', ena.loc[run, 'experiment_title']).group(0) if run in ena.index else ''
    for v, col in [('V0', 'tAge_primary'), ('V2', 'tAge_V2')]:
        rows.append(dict(dataset='GSE66582', analysis='Gate 2b R2', species='Mouse', sample_id=run, gsm=gsm, srr=run,
                         stage=r.stage, arm='', variant=v, reference='MII oocytes', clock_value=r[col]))
g4 = rd('gate3_library_tage')
REF = {'GSE280522': 'control E2C', 'GSE221985': 'control E2C', 'GSE300734': 'control E2C'}
for _, r in g4.iterrows():
    rows.append(dict(dataset=r.gse, analysis='Gate 3 ' + ('primary' if r.gse in REF else 'secondary'), species='Mouse',
                     sample_id=ACC.loc[r.run, 'gsm_title'], gsm=ACC.loc[r.run, 'gsm'], srr=r.run, stage=r.stage, arm=r.arm,
                     variant=r.variant, reference=REF.get(r.gse, 'control libraries'), clock_value=r.tAge))
S3 = pd.DataFrame(rows)
assert len(S3) == 334 + 2 * 39 + 2 * 15 + 115, len(S3)
assert (S3.gsm != '').sum() == len(S3) - 2 * 39, (S3.gsm == '').sum()        # only GSE45719 pseudobulks lack a GSM
# check: Gate 3 P1 control drop (V0) equals the stored arm drop
d = S3[(S3.dataset == 'GSE280522') & (S3.variant == 'V0') & (S3.arm == 'control')].groupby('stage').clock_value.mean()
A = rd('gate3_arm_drops')
want = A[(A.gse == 'GSE280522') & (A.variant == 'V0') & (A.arm == 'control')]['drop'].iloc[0]
assert abs((d['L2C'] - d['E2C']) - want) < 1e-12, (d, want)

# ---- S4 ------------------------------------------------------------------------------------------------
P1 = rd('posthoc_gate3_contribution_genes', index_col=0)
P3 = rd('posthoc_gate3_contribution_genes_P3', index_col=0)
RC = rd('posthoc_rescue_contributions', index_col=0)
S4 = pd.DataFrame({'entrez_id': P1.index, 'symbol': P1.symbol.values, 'clock_coefficient': P1.coef.values,
                   'P1_log2FC_control_L2C_vs_E2C': P1.log2FC_control_L2C_vs_E2C.values,
                   'P1_dynamic_gene_class': P1.label.values,
                   'P1_contribution_control': P1.c_control.values, 'P1_contribution_A485': P1.c_perturbed.values,
                   'P1_contribution_A485_DUX': RC.loc[P1.index, 'A485+Dux'].values,
                   'P3_contribution_control': P3.loc[P1.index, 'c_control'].values,
                   'P3_contribution_Brg1_matKO': P3.loc[P1.index, 'c_perturbed'].values,
                   'P3_log2FC_control_L2C_vs_E2C': P3.loc[P1.index, 'log2FC_control_L2C_vs_E2C'].values})
assert np.allclose(RC.loc[P1.index, 'control'].values, S4.P1_contribution_control.values)
assert np.allclose(RC.loc[P1.index, 'A485'].values, S4.P1_contribution_A485.values)
assert abs(S4.P1_contribution_control.sum() - want) < 1e-9
S4 = S4.sort_values('P1_contribution_control')

# ---- S5 ------------------------------------------------------------------------------------------------
F = rd('posthoc_gate4_crossfit_folds')
S5 = F.rename(columns={'heldout_E2C': 'held_out_control_E2C', 'heldout_L2C': 'held_out_control_L2C',
                       'n_events': 'n_selected_events', 'overlap_with_gate4': 'fraction_shared_with_gate4_events',
                       'P_ctrl_in': 'progression_control_in_sample', 'P_ctrl_out': 'progression_control_held_out',
                       'P_A485': 'progression_A485', 'P_Dux': 'progression_A485_DUX',
                       'I_A485': 'I_A485_minus_held_out_control', 'I_Dux': 'I_A485_DUX_minus_held_out_control',
                       'R': 'R_A485_DUX_minus_A485'})
assert len(S5) == 16 and int((S5.I_A485_minus_held_out_control < 0).sum()) == 14 and int((S5.R_A485_DUX_minus_A485 > 0).sum()) == 16

# ---- write -------------------------------------------------------------------------------------------------
README = pd.DataFrame({'sheet': ['TableS3', 'TableS4', 'TableS5'], 'content': [
    'Clock value (tAge multi-species chronological Elastic Net, scaled differences) of every scored embryo or library. '
    'Values are relative to the reference named in the row (per-gene median of the reference samples). V0 = all genes; '
    'V2 = maternal and zygotic dynamic genes removed (defined within each dataset). GSE45719 rows are embryo pseudobulks '
    'of single cells. Libraries excluded by the prespecified quality-control rule (3 of 76) were not scored. '
    'Rhesus embryos (GSE225056) were scored for description only.',
    'Per-gene contribution (clock coefficient x change in the preprocessed feature, late minus early two-cell) for the '
    '1,839 clock genes present in both primary datasets; contributions sum to the arm drop. P1 = GSE280522, '
    'P3 = GSE300734. Dynamic-gene class: maternal / zygotic / neither, as defined for V2 on the P1 control arm. Post hoc.',
    'Cross-fitted splicing progression in GSE280522 (post hoc). Each fold holds out one control E2C and one control L2C '
    'library, selects events and sets the scale on the remaining controls, and scores all libraries out of sample. '
    'I = arm progression minus held-out control progression; R = A485+DUX minus A485.']})
with pd.ExcelWriter(f'{OUT}/SupplementaryTables_S3-S5.xlsx') as xw:
    README.to_excel(xw, sheet_name='README', index=False)
    S3.to_excel(xw, sheet_name='TableS3', index=False)
    S4.to_excel(xw, sheet_name='TableS4', index=False)
    S5.to_excel(xw, sheet_name='TableS5', index=False)
for name, df in [('TableS3_clock_values', S3), ('TableS4_gene_contributions', S4), ('TableS5_crossfit_folds', S5)]:
    df.to_csv(f'{OUT}/{name}.tsv', sep='\t', index=False)
print('S3', S3.shape, S3.groupby(['dataset', 'variant']).size().to_dict())
print('S4', S4.shape, 'S5', S5.shape)
print('written to', OUT)
