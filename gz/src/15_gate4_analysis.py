"""Gate 4 (plan: plan/GATE4_PLAN_frozen.md + plan/GATE4_IMPLEMENTATION_NOTES.md).
Zygotic splicing activation in P1 GSE280522, P2 GSE221985, P3 GSE300734: SUPPA2 PSI from kallisto transcript TPM,
control ZSA events, splicing progression score, interaction and rescue; clock coupling descriptive. seed 20260919."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import os, re, sys, glob, subprocess
import numpy as np, pandas as pd
from scipy.stats import spearmanr

B = f'{ANALYSIS}/gz'
PY = sys.executable
SUPPA = f'{B}/tools_SUPPA/suppa.py'
W = f'{B}/data/gate4'
os.makedirs(W, exist_ok=True)
rng = np.random.default_rng(20260919)
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
TX2G = dict(zip(t2g.tx.str.replace(r'\.\d+$', '', regex=True), t2g.gene_id))
T2G_V = dict(zip(t2g.tx, t2g.gene_id))
DATASETS = {'GSE280522': ('A485', 'A485+Dux'), 'GSE221985': ('Tardbp_matKO', None), 'GSE300734': ('Brg1_matKO', None)}


def annotate(gse, t):                                   # identical to 11_gate3_analysis.py for these three series
    if gse == 'GSE280522':
        m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep', t)
        return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)]
    if gse == 'GSE221985':
        m = re.search(r'RNA_(E2C|L2C)_(Ctrl|KO)_rep', t)
        return m.group(1), {'Ctrl': 'control', 'KO': 'Tardbp_matKO'}[m.group(2)]
    m = re.search(r'(WT|KO)_(E2C|L2C)_RNA', t)
    return m.group(2), {'WT': 'control', 'KO': 'Brg1_matKO'}[m.group(1)]


def load(gse):
    """QC exactly as Gate 3; returns library table and transcript TPM (version-free IDs)."""
    S = SEL[SEL.gse == gse]
    counts, tpm, meta = {}, {}, []
    for _, r in S.iterrows():
        q = f'{B}/data/perturb/{gse}/quant/{r.run_accession}'
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        p_al = float(re.search(r'"p_pseudoaligned": ([\d.]+)', open(f'{q}/run_info.json').read()).group(1))
        counts[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(T2G_V).values).sum()
        tpm[r.run_accession] = a.tpm.set_axis(a.index.str.replace(r'\.\d+$', '', regex=True))
        stage, arm = annotate(gse, r.title_clean)
        meta.append(dict(run=r.run_accession, stage=stage, arm=arm, p_pseudoaligned=p_al))
    C = pd.DataFrame(counts).fillna(0)
    M = pd.DataFrame(meta).set_index('run')
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    M['qc_pass'] = (M.p_pseudoaligned >= 30) & (det >= det.median() - 3 * mad)
    M = M[M.qc_pass]
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


def run(cmd, log):
    with open(log, 'w') as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f'{cmd[2]} failed, see {log}')


# merged ioe (one header)
IOE = f'{B}/data/annotation/suppa/mm112_all_strict.ioe'
parts = sorted(glob.glob(f'{B}/data/annotation/suppa/mm112_*_strict.ioe'))
parts = [p for p in parts if not p.endswith('all_strict.ioe')]
with open(IOE, 'w') as out:
    for i, p in enumerate(parts):
        with open(p) as fh:
            head = fh.readline()
            if i == 0:
                out.write(head)
            out.writelines(fh)

summary_rows, prog_rows, arm_rows, zsa_n = [], [], [], {}
for gse, (blk, rescue_arm) in DATASETS.items():
    M, TPM = load(gse)
    write_expr(TPM, f'{W}/{gse}_tpm.tsv')
    run([PY, SUPPA, 'psiPerEvent', '-i', IOE, '-e', f'{W}/{gse}_tpm.tsv', '-o', f'{W}/{gse}'], f'{W}/{gse}_psi.log')
    PSI = read_expr(f'{W}/{gse}.psi')[M.index]

    # plan filter
    gene_tpm = TPM.groupby(TPM.index.map(lambda t: TX2G.get(t, t))).sum()
    ev_gene = PSI.index.to_series().str.split(';').str[0]
    ce = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
    cl = M.index[(M.arm == 'control') & (M.stage == 'L2C')]
    gt = gene_tpm.reindex(ev_gene.values)
    gt.index = PSI.index
    expr_ok = ((gt[ce] >= 5).mean(axis=1) >= 0.8) & ((gt[cl] >= 5).mean(axis=1) >= 0.8)
    nan_ok = pd.Series(True, index=PSI.index)
    for (arm, stage), g in M.groupby(['arm', 'stage']):
        nan_ok &= PSI[g.index].notna().mean(axis=1) >= 0.8
    F = PSI[expr_ok & nan_ok]

    # control-only diffSplice on all events, then intersect
    for nm, cols in [('ctrlE2C', ce), ('ctrlL2C', cl)]:
        write_expr(PSI[cols], f'{W}/{gse}_{nm}.psi')
        write_expr(TPM[cols], f'{W}/{gse}_{nm}_tpm.tsv')
    run([PY, SUPPA, 'diffSplice', '-m', 'empirical', '-gc', '-i', IOE,
         '-p', f'{W}/{gse}_ctrlE2C.psi', f'{W}/{gse}_ctrlL2C.psi',
         '-e', f'{W}/{gse}_ctrlE2C_tpm.tsv', f'{W}/{gse}_ctrlL2C_tpm.tsv',
         '-o', f'{W}/{gse}_ctrl'], f'{W}/{gse}_diffsplice.log')
    D = read_expr(f'{W}/{gse}_ctrl.dpsi')
    pcol = [c for c in D.columns if c.endswith('p-val')][0]
    dpsi = F[cl].mean(axis=1) - F[ce].mean(axis=1)
    pval = D[pcol].reindex(F.index)
    zsa = (dpsi.abs() >= 0.10) & (pval < 0.05)
    etype = F.index.to_series().str.split(';').str[1].str.split(':').str[0]
    EV = pd.DataFrame({'event_type': etype, 'dPSI_control': dpsi, 'p_diffsplice': pval, 'zsa': zsa})
    EV.to_csv(f'{B}/results/gate4_events_filtered_{gse}.tsv', sep='\t')
    zsa_n[gse] = int(zsa.sum())
    for t in sorted(etype.unique()):
        m = etype == t
        summary_rows.append(dict(gse=gse, event_type=t, n_filtered=int(m.sum()), n_zsa=int((zsa & m).sum()),
                                 frac_zsa=float((zsa & m).sum() / m.sum()) if m.sum() else np.nan))
    summary_rows.append(dict(gse=gse, event_type='ALL', n_filtered=len(F), n_zsa=zsa_n[gse],
                             frac_zsa=zsa_n[gse] / len(F) if len(F) else np.nan,
                             n_events_total=len(PSI), n_libraries=len(M)))

    # progression score
    Z = F[zsa]
    s = np.sign(dpsi[zsa])
    base = Z[ce].mean(axis=1)
    score = (Z.sub(base, axis=0).mul(s, axis=0).div(dpsi[zsa].abs(), axis=0)).mean(axis=0, skipna=True)
    for run_id, v in score.items():
        prog_rows.append(dict(gse=gse, run=run_id, arm=M.loc[run_id, 'arm'], stage=M.loc[run_id, 'stage'], progression=v))

    def boot(arm):
        e = score[M.index[(M.arm == arm) & (M.stage == 'E2C')]].values
        l = score[M.index[(M.arm == arm) & (M.stage == 'L2C')]].values
        bs = np.array([rng.choice(l, len(l)).mean() - rng.choice(e, len(e)).mean() for _ in range(2000)])
        return l.mean() - e.mean(), bs, len(e), len(l)

    arms = {a: boot(a) for a in ['control', blk] + ([rescue_arm] if rescue_arm else [])}
    for a, (p, bs, ne, nl) in arms.items():
        arm_rows.append(dict(gse=gse, quantity=f'P {a}', value=p, ci_lo=np.percentile(bs, 2.5),
                             ci_hi=np.percentile(bs, 97.5), n_E2C=ne, n_L2C=nl))
    for a in [x for x in arms if x != 'control']:
        bi = arms[a][1] - arms['control'][1]
        arm_rows.append(dict(gse=gse, quantity=f'INTERACTION {a} - control', value=arms[a][0] - arms['control'][0],
                             ci_lo=np.percentile(bi, 2.5), ci_hi=np.percentile(bi, 97.5)))
    if rescue_arm:
        br = arms[rescue_arm][1] - arms[blk][1]
        arm_rows.append(dict(gse=gse, quantity=f'RESCUE {rescue_arm} - {blk}', value=arms[rescue_arm][0] - arms[blk][0],
                             ci_lo=np.percentile(br, 2.5), ci_hi=np.percentile(br, 97.5)))

SUM = pd.DataFrame(summary_rows)
PRO = pd.DataFrame(prog_rows)
ARM = pd.DataFrame(arm_rows)
SUM.to_csv(f'{B}/results/gate4_zsa_summary.tsv', sep='\t', index=False)
PRO.to_csv(f'{B}/results/gate4_progression.tsv', sep='\t', index=False)
ARM.to_csv(f'{B}/results/gate4_arm_progression.tsv', sep='\t', index=False)

# decision (plan bullets, literal)
inter = {g: float(ARM[(ARM.gse == g) & (ARM.quantity == f'INTERACTION {DATASETS[g][0]} - control')].value.iloc[0])
         for g in DATASETS}
rescue = float(ARM[ARM.quantity.str.startswith('RESCUE')].value.iloc[0])
present = sum(n >= 100 for n in zsa_n.values()) >= 2
neg = sum(v < 0 for v in inter.values())
nonneg = sum(v >= 0 for v in inter.values())
if present and neg >= 2 and rescue > 0:
    verdict = 'ZSA ZGA-DEPENDENT'
elif nonneg >= 2:
    verdict = 'ZSA NOT ZGA-DEPENDENT'
else:
    verdict = 'MIXED'

# clock coupling (descriptive)
L3 = pd.read_csv(f'{B}/results/gate3_library_tage.tsv', sep='\t')
L3 = L3[(L3.variant == 'V0') & L3.gse.isin(DATASETS)]
gc = L3.groupby(['gse', 'arm', 'stage']).tAge.mean()
gp = PRO.groupby(['gse', 'arm', 'stage']).progression.mean()
J = pd.concat([gp, gc], axis=1, join='inner')
coup = {g: spearmanr(J.loc[g].progression, J.loc[g].tAge).correlation for g in DATASETS}
coup['pooled'] = spearmanr(J.progression, J.tAge).correlation
A3 = pd.read_csv(f'{B}/results/gate3_arm_drops.tsv', sep='\t')
A3 = A3[(A3.variant == 'V0') & A3.arm.str.startswith('INTERACTION')]
pairs = []
for _, r in ARM[ARM.quantity.str.startswith('INTERACTION')].iterrows():
    c = A3[(A3.gse == r.gse) & (A3.arm == r.quantity)]
    if len(c):
        pairs.append((r.gse, r.quantity, r.value, float(c['drop'].iloc[0])))
PI = pd.DataFrame(pairs, columns=['gse', 'interaction', 'progression_I', 'clock_I_V0'])
rho_i = spearmanr(PI.progression_I, PI.clock_I_V0).correlation if len(PI) >= 3 else np.nan

pd.set_option('display.width', 250)
print(SUM.round(4).to_string(index=False))
print(ARM.round(4).to_string(index=False))
print('ZSA events per dataset:', zsa_n, '| PRESENT:', present)
print('interactions:', {k: round(v, 4) for k, v in inter.items()}, '| rescue:', round(rescue, 4))
print('VERDICT:', verdict)
print('coupling (Spearman, group means, progression vs tAge V0):', {k: round(v, 3) for k, v in coup.items()})
print(PI.round(4).to_string(index=False))
print('Spearman of interactions (progression vs clock):', round(rho_i, 3))
with open(f'{B}/results/gate4_decision.md', 'w', encoding='utf-8') as fh:
    fh.write('# Gate 4 decision\n\n')
    fh.write(f'VERDICT: **{verdict}**' + ('' if present else ' — note: ZSA PRESENT not met; see counts') + '\n\n')
    fh.write(f'ZSA events per dataset: {zsa_n}; ZSA PRESENT: {present}\n\n')
    fh.write(f'Interactions (progression, blocking arm): {inter}; P1 rescue: {rescue}\n\n')
    fh.write('## ZSA summary\n\n' + SUM.round(4).to_string(index=False) + '\n\n')
    fh.write('## Arm progression\n\n' + ARM.round(4).to_string(index=False) + '\n\n')
    fh.write('## Clock coupling (descriptive, no inference)\n\n')
    fh.write(f'Spearman, arm x stage group means, progression vs Gate 3 V0 tAge: {coup}\n\n')
    fh.write(PI.round(4).to_string(index=False) + f'\n\nSpearman of interactions: {rho_i}\n')
