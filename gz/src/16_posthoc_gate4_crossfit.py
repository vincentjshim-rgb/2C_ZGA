"""POST HOC (plan: plan/POSTHOC_gate4_crossfit_frozen.md). Cross-fitted check of the Gate 4 progression score in
P1 GSE280522: 16 folds, each holding out one control E2C and one control L2C library; events selected and scaled on the
remaining 3 + 3 control libraries. Loading, filter and scoring mirror src/15_gate4_analysis.py."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import os, re, sys, itertools, subprocess
import numpy as np, pandas as pd

B = f'{ANALYSIS}/gz'
PY = sys.executable
SUPPA = f'{B}/tools_SUPPA/suppa.py'
IOE = f'{B}/data/annotation/suppa/mm112_all_strict.ioe'
W = f'{B}/data/gate4/crossfit'
os.makedirs(W, exist_ok=True)
GSE = 'GSE280522'
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
TX2G = dict(zip(t2g.tx.str.replace(r'\.\d+$', '', regex=True), t2g.gene_id))
T2G_V = dict(zip(t2g.tx, t2g.gene_id))


def annotate(t):
    m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep', t)
    return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)]


def load():
    S = SEL[SEL.gse == GSE]
    counts, tpm, meta = {}, {}, []
    for _, r in S.iterrows():
        q = f'{B}/data/perturb/{GSE}/quant/{r.run_accession}'
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        p_al = float(re.search(r'"p_pseudoaligned": ([\d.]+)', open(f'{q}/run_info.json').read()).group(1))
        counts[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(T2G_V).values).sum()
        tpm[r.run_accession] = a.tpm.set_axis(a.index.str.replace(r'\.\d+$', '', regex=True))
        stage, arm = annotate(r.title_clean)
        meta.append(dict(run=r.run_accession, stage=stage, arm=arm, p_pseudoaligned=p_al))
    C = pd.DataFrame(counts).fillna(0)
    M = pd.DataFrame(meta).set_index('run')
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    M = M[(M.p_pseudoaligned >= 30) & (det >= det.median() - 3 * mad)]
    return M, pd.DataFrame(tpm)[M.index]


def write_expr(df, path):
    with open(path, 'w') as fh:
        fh.write('\t'.join(df.columns) + '\n')
        df.to_csv(fh, sep='\t', header=False, float_format='%.6g')


def read_expr(path):
    with open(path) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
    df = pd.read_csv(path, sep='\t', skiprows=1, header=None, index_col=0)
    df.columns = cols
    return df


M, TPM = load()
PSI = read_expr(f'{B}/data/gate4/{GSE}.psi')[M.index]
gene_tpm = TPM.groupby(TPM.index.map(lambda t: TX2G.get(t, t))).sum()
GT = gene_tpm.reindex(PSI.index.to_series().str.split(';').str[0].values)
GT.index = PSI.index
lib = lambda arm, st: list(M.index[(M.arm == arm) & (M.stage == st)])
CE, CL = lib('control', 'E2C'), lib('control', 'L2C')
PERT = {a: (lib(a, 'E2C'), lib(a, 'L2C')) for a in ['A485', 'A485+Dux']}
G4_EVENTS = set(pd.read_csv(f'{B}/results/gate4_events_filtered_{GSE}.tsv', sep='\t', index_col=0).query('zsa').index)


def fold(sel_e, sel_l, tag):
    expr_ok = ((GT[sel_e] >= 5).mean(axis=1) >= 0.8) & ((GT[sel_l] >= 5).mean(axis=1) >= 0.8)
    nan_ok = (PSI[sel_e].notna().mean(axis=1) >= 0.8) & (PSI[sel_l].notna().mean(axis=1) >= 0.8)
    for e, l in PERT.values():
        nan_ok &= (PSI[e].notna().mean(axis=1) >= 0.8) & (PSI[l].notna().mean(axis=1) >= 0.8)
    F = PSI[expr_ok & nan_ok]
    write_expr(PSI[sel_e], f'{W}/{tag}_E.psi'); write_expr(TPM[sel_e], f'{W}/{tag}_E_tpm.tsv')
    write_expr(PSI[sel_l], f'{W}/{tag}_L.psi'); write_expr(TPM[sel_l], f'{W}/{tag}_L_tpm.tsv')
    with open(f'{W}/{tag}_diffsplice.log', 'w') as fh:
        r = subprocess.run([PY, SUPPA, 'diffSplice', '-m', 'empirical', '-gc', '-i', IOE,
                            '-p', f'{W}/{tag}_E.psi', f'{W}/{tag}_L.psi',
                            '-e', f'{W}/{tag}_E_tpm.tsv', f'{W}/{tag}_L_tpm.tsv', '-o', f'{W}/{tag}'],
                           stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f'diffSplice failed for {tag}')
    D = read_expr(f'{W}/{tag}.dpsi')
    p = D[[c for c in D.columns if c.endswith('p-val')][0]].reindex(F.index)
    d = F[sel_l].mean(axis=1) - F[sel_e].mean(axis=1)
    z = (d.abs() >= 0.10) & (p < 0.05)
    Z = F[z]
    score = Z.sub(Z[sel_e].mean(axis=1), axis=0).mul(np.sign(d[z]), axis=0).div(d[z].abs(), axis=0).mean(axis=0, skipna=True)
    return set(Z.index), score


def arm_P(score, e, l):
    return score[l].mean() - score[e].mean()


# reproduction check: no hold-out must give the Gate 4 values
ev0, s0 = fold(CE, CL, 'all')
chk = dict(control=arm_P(s0, CE, CL), A485=arm_P(s0, *PERT['A485']), Dux=arm_P(s0, *PERT['A485+Dux']))
ref = dict(control=1.0, A485=0.5568, Dux=0.7869)
print('reproduction check:', {k: round(v, 6) for k, v in chk.items()}, '| events', len(ev0), '(Gate 4:', len(G4_EVENTS), ')')
G4 = pd.read_csv(f'{B}/results/gate4_arm_progression.tsv', sep='\t')
G4 = G4[G4.gse == GSE].set_index('quantity').value
exact = dict(control=G4['P control'], A485=G4['P A485'], Dux=G4['P A485+Dux'])
for k in chk:
    if abs(chk[k] - exact[k]) > 1e-6:
        sys.exit(f'REPRODUCTION FAILED for {k}: {chk[k]} vs Gate 4 {exact[k]}')
if ev0 != G4_EVENTS:
    sys.exit('REPRODUCTION FAILED: event set differs from Gate 4')
print('reproduction check passed (1e-6, identical event set)')

rows = []
for i, (he, hl) in enumerate(itertools.product(CE, CL)):
    se, sl = [x for x in CE if x != he], [x for x in CL if x != hl]
    ev, sc = fold(se, sl, f'f{i:02d}')
    p_out = sc[hl] - sc[he]
    p_in = arm_P(sc, se, sl)
    pa, pd_ = arm_P(sc, *PERT['A485']), arm_P(sc, *PERT['A485+Dux'])
    rows.append(dict(fold=i, heldout_E2C=he, heldout_L2C=hl, n_events=len(ev),
                     overlap_with_gate4=len(ev & G4_EVENTS) / len(ev) if ev else np.nan,
                     P_ctrl_in=p_in, P_ctrl_out=p_out, P_A485=pa, P_Dux=pd_,
                     I_A485=pa - p_out, I_Dux=pd_ - p_out, R=pd_ - pa))
    print(f'fold {i:02d}: events {len(ev):4d}  P_in {p_in:.3f}  P_out {p_out:.3f}  P_A485 {pa:.3f}  I_A485 {pa - p_out:+.3f}')
FD = pd.DataFrame(rows)
FD.to_csv(f'{B}/results/posthoc_gate4_crossfit_folds.tsv', sep='\t', index=False)

summ = {}
for q in ['n_events', 'overlap_with_gate4', 'P_ctrl_in', 'P_ctrl_out', 'P_A485', 'P_Dux', 'I_A485', 'I_Dux', 'R']:
    v = FD[q].values
    summ[q] = (v.mean(), np.percentile(v, 2.5), np.percentile(v, 97.5))
I_cf = summ['I_A485'][0]
n_neg = int((FD.I_A485 < 0).sum())
if I_cf <= -0.20 and n_neg >= 15:
    reading = 'SURVIVES'
elif I_cf < 0:
    reading = 'ATTENUATED BUT INFLATED'
else:
    reading = 'CONSTRUCTION'
inflation = 1 - summ['P_ctrl_out'][0]

print('\n=== summary over 16 folds (mean; fold spread 2.5-97.5%, not a CI) ===')
for q, (m, lo, hi) in summ.items():
    print(f'{q:20s} {m:8.4f}   [{lo:8.4f}, {hi:8.4f}]')
print(f'I_A485 negative in {n_neg}/16 folds; inflation of in-sample control = {inflation:.4f}')
print('READING:', reading)

with open(f'{B}/results/posthoc_gate4_crossfit.md', 'w', encoding='utf-8') as fh:
    fh.write('# POST HOC — cross-fitted check of the Gate 4 progression score (P1 GSE280522)\n\n')
    fh.write('Plan `plan/POSTHOC_gate4_crossfit_frozen.md`. Reproduction check (no hold-out) matched Gate 4 to 1e-6 with an '
             'identical event set.\n\n')
    fh.write(f'READING: **{reading}** — mean I_A485 {I_cf:.4f}, negative in {n_neg}/16 folds; '
             f'in-sample inflation of the control arm {inflation:.4f}.\n\n')
    fh.write('| quantity | mean over folds | fold spread (2.5–97.5%) |\n|---|---|---|\n')
    for q, (m, lo, hi) in summ.items():
        fh.write(f'| {q} | {m:.4f} | {lo:.4f} – {hi:.4f} |\n')
    fh.write('\nFold spread is not a confidence interval: folds share libraries. P2 and P3 are not cross-checked.\n')
