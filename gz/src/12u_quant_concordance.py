"""Technical concordance of kallisto outputs for the same run (no clock score, no stage/arm labels used).
Usage: python 12u_quant_concordance.py <label> <quant_dir_A> <quant_dir_B> [<label> <A> <B> ...]
Compares transcript est_counts / TPM and Ensembl-gene summed counts; appends to results/ubuntu_quant_concordance.tsv."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import sys, re, os
import numpy as np, pandas as pd

B = f'{ANALYSIS}/gz'
t2g = pd.read_csv(f'{B}/../metadata/tx2gene_full.tsv', sep='\t', usecols=['tx', 'gene_id'])
t2g = dict(zip(t2g.tx, t2g.gene_id))
OUT = f'{B}/results/ubuntu_quant_concordance.tsv'


def load(q):
    a = pd.read_csv(f'{q}/abundance.tsv', sep='\t', index_col=0)
    # Windows run_info.json has unescaped backslashes in "call"; read fields by regex (as in 11_gate3_analysis.py)
    txt = open(f'{q}/run_info.json').read()
    info = {k: float(re.search(rf'"{k}": ([\d.]+)', txt).group(1)) for k in ['n_processed', 'p_pseudoaligned']}
    return a, info


rows = []
args = sys.argv[1:]
for i in range(0, len(args), 3):
    label, qa, qb = args[i:i + 3]
    a, ia = load(qa)
    b, ib = load(qb)
    assert a.index.equals(b.index), 'target order differs'
    ga = a.est_counts.groupby(a.index.map(t2g)).sum()
    gb = b.est_counts.groupby(b.index.map(t2g)).sum()
    expressed = ga.index[(ga >= 10) | (gb >= 10)]
    lfc = np.log2(ga[expressed] + 1) - np.log2(gb[expressed] + 1)
    rows.append(dict(
        comparison=label, run=os.path.basename(qa.rstrip('/')),
        n_processed_A=ia['n_processed'], n_processed_B=ib['n_processed'],
        p_pseudoaligned_A=ia['p_pseudoaligned'], p_pseudoaligned_B=ib['p_pseudoaligned'],
        tx_est_counts_identical=bool((a.est_counts == b.est_counts).all()),
        tx_max_abs_diff_counts=float((a.est_counts - b.est_counts).abs().max()),
        tx_pearson_tpm=float(np.corrcoef(a.tpm, b.tpm)[0, 1]),
        gene_pearson_log_counts=float(np.corrcoef(np.log1p(ga), np.log1p(gb))[0, 1]),
        gene_n_expressed_ge10=len(expressed),
        gene_max_abs_log2_diff=float(lfc.abs().max()),
        gene_frac_abs_log2_diff_gt_0p1=float((lfc.abs() > 0.1).mean()),
        eff_len_max_abs_diff=float((a.eff_length - b.eff_length).abs().max())))
R = pd.DataFrame(rows)
R.to_csv(OUT, sep='\t', index=False, mode='a', header=not os.path.exists(OUT))
pd.set_option('display.width', 250)
print(R.T.to_string())
