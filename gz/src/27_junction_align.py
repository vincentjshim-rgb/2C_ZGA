"""POST HOC (plan: plan/POSTHOC_junction_psi_frozen.md). Align the GSE280522 libraries with STAR so that splice
junctions can be counted from reads, independently of the kallisto/SUPPA2 quantification used for Gate 4.

Only the libraries that passed the prespecified Gate 3 quality control are aligned — the same set that was scored.
Two-pass mapping per sample; only SJ.out.tab (uniquely mapped junction reads) and the sorted BAM are kept.
"""
import os
import subprocess
import sys

import pandas as pd

ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
B = f'{ANALYSIS}/gz'
GSE = 'GSE280522'
IDX = f'{B}/ref/star_idx'
OUT = f'{B}/align/{GSE}'
FQ = f'{B}/data/perturb/{GSE}/fastq'
THREADS = 20

Q = pd.read_csv(f'{B}/figures/source_data/figS7_library_qc.tsv', sep='\t')
runs = sorted(Q.loc[(Q.gse == GSE) & Q.qc_pass, 'run'])
assert len(runs) == 23, f'expected the 23 scored libraries, got {len(runs)}'

os.makedirs(OUT, exist_ok=True)
for i, run in enumerate(runs, 1):
    done = f'{OUT}/{run}/SJ.out.tab'
    if os.path.exists(done):
        print(f'[{i}/{len(runs)}] {run} already aligned', flush=True)
        continue
    os.makedirs(f'{OUT}/{run}', exist_ok=True)
    r1, r2 = f'{FQ}/{run}_1.fastq.gz', f'{FQ}/{run}_2.fastq.gz'
    assert os.path.exists(r1) and os.path.exists(r2), run
    cmd = ['STAR', '--runThreadN', str(THREADS), '--genomeDir', IDX,
           '--readFilesIn', r1, r2, '--readFilesCommand', 'zcat',
           '--outFileNamePrefix', f'{OUT}/{run}/', '--twopassMode', 'Basic',
           '--outSAMtype', 'BAM', 'SortedByCoordinate',
           '--outSJfilterReads', 'Unique', '--outSAMattrIHstart', '0',
           '--outBAMsortingThreadN', '6', '--limitBAMsortRAM', '20000000000']
    print(f'[{i}/{len(runs)}] {run}', flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f'STAR failed on {run}:\n{r.stderr[-2000:]}')
    # the two-pass working directories are large and not needed afterwards
    for junk in ('_STARpass1', '_STARgenome', '_STARtmp'):
        subprocess.run(['rm', '-rf', f'{OUT}/{run}/{junk}'])

n = sum(os.path.exists(f'{OUT}/{r}/SJ.out.tab') for r in runs)
print(f'\naligned {n} of {len(runs)} libraries; junctions in {OUT}/<run>/SJ.out.tab')
