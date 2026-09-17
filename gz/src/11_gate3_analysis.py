"""Gate 3 (plan: plan/GATE3_PLAN_frozen.md + addenda 1, 2). Does the within-2-cell tAge decrease depend on ZGA?
Primary: P1 GSE280522, P2 GSE221985, P3 GSE300734. Secondary single-stage contrasts: GSE248499, GSE235547 (late 2-cell)
and GSE162345 (addendum 2: window-matched alpha-amanitin vs control at 45 and 54 hpi; not a stage pair).
seed 20260918."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys, re, glob, os
import numpy as np, pandas as pd
sys.path.insert(0, f'{ANALYSIS}/gz/src')
import tage_py as tp

B = f'{ANALYSIS}/gz'
MODEL = f'{B}/tools_tAge_models/EN_Chronoage_Multispecies_Multitissue_scaleddiff.pkl'
rng = np.random.default_rng(20260918)
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
t2g = pd.read_csv(f'{ANALYSIS}/metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
t2g = dict(zip(t2g.tx, t2g.gene_id))


def annotate(gse, title):
    t = title
    if gse == 'GSE280522':
        m = re.search(r'RNA_(E2C|L2C)_(DMSO|A485\+Dux|A485)_rep', t)
        return m.group(1), {'DMSO': 'control', 'A485': 'A485', 'A485+Dux': 'A485+Dux'}[m.group(2)]
    if gse == 'GSE221985':
        m = re.search(r'RNA_(E2C|L2C)_(Ctrl|KO)_rep', t)
        return m.group(1), {'Ctrl': 'control', 'KO': 'Tardbp_matKO'}[m.group(2)]
    if gse == 'GSE300734':
        m = re.search(r'(WT|KO)_(E2C|L2C)_RNA', t)
        return m.group(2), {'WT': 'control', 'KO': 'Brg1_matKO'}[m.group(1)]
    if gse == 'GSE162345':
        # addendum 2: "early"/"late" = recipient stage at transfer (21/30 hpi), collected 45/54 hpi. Not a stage pair;
        # the arm carries its collection window so that contrasts stay within a window.
        m = re.search(r'(aAm_)?(early|late)2cell_NT\d', t)
        win = {'early': '45hpi', 'late': '54hpi'}[m.group(2)]
        return f'NT{win}', ('alpha-amanitin_' if m.group(1) else 'control_') + win
    if gse == 'GSE248499':
        arm = 'IVF' if 'IVF' in t else 'SCNT_control' if 'cont' in t else 'SCNT_Kdm4d' if 'Kdm4d' in t else 'SCNT_Kdm3a'
        return 'L2C', ('control' if arm == 'IVF' else arm)
    if gse == 'GSE235547':
        for k, v in [('ICSI', 'control_ICSI'), ('SCNT EGFP', 'SCNT_EGFP'), ('SCNT Obox3', 'SCNT_Obox3'),
                     ('siCtrl', 'control_siCtrl'), ('siObox3', 'siObox3')]:
            if k in t:
                return 'L2C', v
    raise ValueError((gse, t))


def load_dataset(gse):
    S = SEL[SEL.gse == gse]
    cols, meta = {}, []
    for _, r in S.iterrows():
        q = f'{B}/data/perturb/{gse}/quant/{r.run_accession}'
        if not os.path.exists(f'{q}/DONE'):
            continue
        a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
        info = open(f'{q}/run_info.json').read()
        p_al = float(re.search(r'"p_pseudoaligned": ([\d.]+)', info).group(1))
        cols[r.run_accession] = a.est_counts.groupby(a.index.to_series().map(t2g).values).sum()
        stage, arm = annotate(gse, r.title_clean)
        meta.append(dict(run=r.run_accession, stage=stage, arm=arm, p_pseudoaligned=p_al, expected=len(S)))
    C = pd.DataFrame(cols).fillna(0)
    M = pd.DataFrame(meta).set_index('run')
    det = np.log10((C >= 1).sum(axis=0))
    mad = 1.4826 * np.median(np.abs(det - det.median()))
    M['detected_log10'] = det
    M['qc_pass'] = (M.p_pseudoaligned >= 30) & (det >= det.median() - 3 * mad)
    return C, M


def dynamic_genes(C, M):
    cpm = C / C.sum(axis=0) * 1e6
    e = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
    l = M.index[(M.arm == 'control') & (M.stage == 'L2C')]
    mat = ((cpm[e] >= 10).mean(axis=1) >= 0.8) & (cpm[l].median(axis=1) <= cpm[e].median(axis=1) / 8)
    zyg = ((cpm[e] < 1).mean(axis=1) >= 0.8) & ((cpm[l] >= 10).mean(axis=1) >= 0.5)
    return set(C.index[mat | zyg])


def score(C, ref):
    pp = tp.preprocess(C, reference_samples=list(ref))
    return tp.predict(MODEL, pp['scaled_diff'], species=None)


def boot_mean(v):
    return np.array([rng.choice(v, len(v)).mean() for _ in range(2000)])


rows_lib, rows_arm, rows_sec, decisions = [], [], [], {}
for gse, kind in [('GSE280522', 'P'), ('GSE221985', 'P'), ('GSE300734', 'P'),
                  ('GSE248499', 'S'), ('GSE235547', 'S'), ('GSE162345', 'S')]:
    C, M = load_dataset(gse)
    M = M[M.qc_pass]
    C = C[M.index]
    for variant in ['V0', 'V2']:
        if kind == 'P':
            ref = M.index[(M.arm == 'control') & (M.stage == 'E2C')]
            Cv = C if variant == 'V0' else C.loc[~C.index.isin(dynamic_genes(C, M))]
            v = score(Cv, ref)
            for run, x in v.items():
                rows_lib.append(dict(gse=gse, variant=variant, run=run, stage=M.loc[run, 'stage'], arm=M.loc[run, 'arm'], tAge=x))
            drops = {}
            for arm in M.arm.unique():
                e, l = v[(M.arm == arm) & (M.stage == 'E2C')].values, v[(M.arm == arm) & (M.stage == 'L2C')].values
                if len(e) == 0 or len(l) == 0:
                    continue
                bd = boot_mean(l) - boot_mean(e)
                drops[arm] = (l.mean() - e.mean(), bd)
                rows_arm.append(dict(gse=gse, variant=variant, arm=arm, n_E2C=len(e), n_L2C=len(l), drop=l.mean() - e.mean(),
                                     ci_lo=np.percentile(bd, 2.5), ci_hi=np.percentile(bd, 97.5)))
            for arm in drops:
                if arm == 'control':
                    continue
                bi = drops[arm][1] - drops['control'][1]
                rows_arm.append(dict(gse=gse, variant=variant, arm=f'INTERACTION {arm} - control', drop=drops[arm][0] - drops['control'][0],
                                     ci_lo=np.percentile(bi, 2.5), ci_hi=np.percentile(bi, 97.5)))
            if 'A485+Dux' in drops:
                br = drops['A485+Dux'][1] - drops['A485'][1]
                rows_arm.append(dict(gse=gse, variant=variant, arm='RESCUE A485+Dux - A485', drop=drops['A485+Dux'][0] - drops['A485'][0],
                                     ci_lo=np.percentile(br, 2.5), ci_hi=np.percentile(br, 97.5)))
        else:
            ctrl = M.index[M.arm.str.startswith('control')]
            if variant == 'V2':
                continue          # V2 needs E2C/L2C control arms; not defined for single-stage secondaries
            v = score(C, ctrl)
            for run, x in v.items():
                rows_lib.append(dict(gse=gse, variant=variant, run=run, stage=M.loc[run, 'stage'], arm=M.loc[run, 'arm'], tAge=x))
            for arm in M.arm.unique():
                if arm.startswith('control'):
                    continue
                x = v[M.arm == arm].values
                if gse == 'GSE162345':
                    ref_arm = 'control_' + arm.split('_')[-1]          # addendum 2: same collection window
                elif gse == 'GSE248499':
                    ref_arm = 'control'                                # IVF; the Kdm4d/Kdm3a arms are read against it
                else:                                                  # GSE235547
                    ref_arm = 'control_siCtrl' if arm == 'siObox3' else 'control_ICSI'
                c = v[M.arm == ref_arm].values
                if len(c) == 0:
                    raise ValueError(f'{gse}: reference arm {ref_arm} has no libraries for arm {arm}')
                bd = boot_mean(x) - boot_mean(c)
                rows_sec.append(dict(gse=gse, arm=arm, reference=ref_arm, n=len(x), n_ref=len(c), diff=x.mean() - c.mean(),
                                     ci_lo=np.percentile(bd, 2.5), ci_hi=np.percentile(bd, 97.5)))

L = pd.DataFrame(rows_lib)
A = pd.DataFrame(rows_arm)
S2 = pd.DataFrame(rows_sec)
L.to_csv(f'{B}/results/gate3_library_tage.tsv', sep='\t', index=False)
A.to_csv(f'{B}/results/gate3_arm_drops.tsv', sep='\t', index=False)
S2.to_csv(f'{B}/results/gate3_secondary_single_stage.tsv', sep='\t', index=False)

verdicts = {}
for variant in ['V0', 'V2']:
    Av = A[A.variant == variant]
    informative, pos = [], 0
    for gse in ['GSE280522', 'GSE221985', 'GSE300734']:      # addendum 2: P4 GSE162345 is no longer primary
        ctrl = Av[(Av.gse == gse) & (Av.arm == 'control')]
        if len(ctrl) and ctrl['drop'].iloc[0] < 0:
            informative.append(gse)
            inter = Av[(Av.gse == gse) & Av.arm.str.startswith('INTERACTION')]
            # P1: the ZGA-blocking arm is A485 (A485+Dux is the rescue)
            if gse == 'GSE280522':
                inter = inter[inter.arm == 'INTERACTION A485 - control']
            pos += int((inter['drop'] > 0).all())
    neg = len(informative) - pos
    rescue = Av[Av.arm == 'RESCUE A485+Dux - A485']
    rescue_ok = bool(len(rescue) and rescue['drop'].iloc[0] < 0)
    need = 3 if len(informative) >= 4 else len(informative)
    if len(informative) >= 2 and pos >= need and rescue_ok:
        v = 'ZGA-DEPENDENT'
    elif neg >= 2:
        v = 'NOT ZGA-DEPENDENT'
    else:
        v = 'MIXED'
    verdicts[variant] = dict(verdict=v, informative=informative, attenuated=pos, not_attenuated=neg, rescue_restores=rescue_ok)
overall = 'ZGA-DEPENDENT' if all(x['verdict'] == 'ZGA-DEPENDENT' for x in verdicts.values()) else \
          'NOT ZGA-DEPENDENT' if all(x['verdict'] == 'NOT ZGA-DEPENDENT' for x in verdicts.values()) else 'MIXED'

pd.set_option('display.width', 250)
print(A.round(4).to_string(index=False))
print(S2.round(4).to_string(index=False))
print(verdicts)
print('OVERALL (requires V0 and V2 agreement):', overall)
with open(f'{B}/results/gate3_decision.md', 'w', encoding='utf-8') as fh:
    fh.write(f'# Gate 3 decision\n\nOVERALL: **{overall}**\n\n{verdicts}\n\n')
    fh.write(A.round(4).to_string(index=False) + '\n\nSecondary (late 2-cell):\n' + S2.round(4).to_string(index=False) + '\n')
