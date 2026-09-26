"""Source data for the junction-level validation figure (post hoc; plan POSTHOC_junction_psi_frozen.md).

  figJ_junction_vs_tpm.tsv  binned density of junction PSI against transcript-TPM PSI over all scored pairs
  figJ_sashimi.tsv          per-arm read coverage and junction counts for the three events the plan selects
"""
import os
import subprocess

import numpy as np
import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
GSE = 'GSE280522'
ALIGN = f'{B}/align/{GSE}'
SD = f'{B}/figures/source_data'
ARMS = ['control', 'A485', 'A485+Dux']
NBIN = 40

J = pd.read_csv(f'{B}/results/posthoc_junction_psi.tsv', sep='\t')
ok = J.dropna(subset=['psi_tpm'])

# ---- (a) 2-D density, so the figure ships a table rather than 59,482 points
h, xe, ye = np.histogram2d(ok.psi_tpm, ok.psi_junction, bins=NBIN, range=[[0, 1], [0, 1]])
dens = pd.DataFrame([(xe[i], xe[i + 1], ye[j], ye[j + 1], int(h[i, j]))
                     for i in range(NBIN) for j in range(NBIN) if h[i, j] > 0],
                    columns=['tpm_lo', 'tpm_hi', 'junction_lo', 'junction_hi', 'n'])
dens.to_csv(f'{SD}/figJ_junction_vs_tpm.tsv', sep='\t', index=False)
print(f'wrote figJ_junction_vs_tpm.tsv: {len(dens)} occupied bins, {int(h.sum()):,} pairs')

# ---- (b) coverage and junctions for the three events selected by the plan's rule
EVENTS = [l.split('`')[1] for l in open(f'{B}/results/posthoc_junction_psi.md') if l.startswith('- `')]
assert len(EVENTS) == 3, EVENTS
Mc = pd.read_csv(f'{SD}/fig4a_heatmap_columns.tsv', sep='\t')[['run', 'arm', 'stage']].drop_duplicates()
late = Mc[Mc.stage == 'L2C']

rows = []
for ev in EVENTS:
    gene, rest = ev.split(';')
    _, chrom, j1, j2, strand = rest.split(':')
    e1, s2 = (int(v) for v in j1.split('-'))
    e2, s3 = (int(v) for v in j2.split('-'))
    lo, hi = min(e1, s3) - 200, max(e1, s3) + 200
    for arm in ARMS:
        runs = list(late.run[late.arm == arm])
        cov = np.zeros(hi - lo + 1)
        for r in runs:
            bam = f'{ALIGN}/{r}/Aligned.sortedByCoord.out.bam'
            if not os.path.exists(f'{bam}.bai'):
                subprocess.run(['samtools', 'index', bam], check=True)
            out = subprocess.run(['samtools', 'depth', '-a', '-Q', '255', '-r', f'{chrom}:{lo}-{hi}', bam],
                                 capture_output=True, text=True, check=True).stdout
            if out:
                d = np.loadtxt(out.splitlines(), dtype=int, usecols=(1, 2), ndmin=2)
                cov[d[:, 0] - lo] += d[:, 1]
        cov /= max(1, len(runs))                       # mean depth per library in the arm
        sj = {}
        for r in runs:
            t = pd.read_csv(f'{ALIGN}/{r}/SJ.out.tab', sep='\t', header=None, dtype={0: str},
                            names=['chrom', 'start', 'end', 'strand', 'motif', 'annot', 'nu', 'nm', 'oh'])
            t = t[t.chrom == chrom]
            for k, v in zip(zip(t.start, t.end), t.nu):
                sj[k] = sj.get(k, 0) + v
        n = len(runs)
        junc = {'inc1': sj.get((e1 + 1, s2 - 1), 0) / n, 'inc2': sj.get((e2 + 1, s3 - 1), 0) / n,
                'skip': sj.get((e1 + 1, s3 - 1), 0) / n}
        step = max(1, (hi - lo) // 1200)               # thin the coverage to keep the table small
        for p in range(0, hi - lo + 1, step):
            rows.append((ev, gene, chrom, strand, arm, e1, s2, e2, s3, lo + p, cov[p],
                         junc['inc1'], junc['inc2'], junc['skip']))
        print(f'{gene} {arm}: inc1 {junc["inc1"]:.0f}  inc2 {junc["inc2"]:.0f}  skip {junc["skip"]:.0f}'
              f'  (mean over {n} late two-cell libraries)')

S = pd.DataFrame(rows, columns=['event', 'gene', 'chrom', 'strand', 'arm', 'e1', 's2', 'e2', 's3',
                                'pos', 'coverage', 'inc1', 'inc2', 'skip'])
S.to_csv(f'{SD}/figJ_sashimi.tsv', sep='\t', index=False)
print(f'wrote figJ_sashimi.tsv: {len(S):,} rows, {S.event.nunique()} events x {S.arm.nunique()} arms')
