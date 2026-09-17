"""Ubuntu worker lists for Gate 3 downloads (same 76 runs as gate3_selected_runs.tsv; execution order only).
Adds ENA fastq_bytes/fastq_md5 (metadata/perturb/gate3_ena_md5.tsv, fetched 2026-09-15) so 10u can verify each file.
Order: the 3 runs quantified on Windows first (cross-index concordance check), then P2, P1, P3, S2, S1, P4 (P4 on hold,
see audit log 2026-09-15). Runs are dealt round-robin to 3 workers so datasets complete roughly in this order."""
import os
ANALYSIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # .../analysis
import pandas as pd

B = f'{ANALYSIS}/gz'
SEL = pd.read_csv(f'{B}/metadata/perturb/gate3_selected_runs.tsv', sep='\t')
MD5 = pd.read_csv(f'{B}/metadata/perturb/gate3_ena_md5.tsv', sep='\t')
WINDOWS_QUANT = ['SRR30827777', 'SRR16201705', 'SRR16201707']
PRIORITY = ['GSE221985', 'GSE280522', 'GSE300734', 'GSE235547', 'GSE248499', 'GSE162345']

x = SEL.merge(MD5, on='run_accession', suffixes=('', '_ena'), validate='one_to_one')
assert len(x) == 76 and (x.library_layout == x.library_layout_ena).all() and (x.read_count == x.read_count_ena).all()
x['first'] = ~x.run_accession.isin(WINDOWS_QUANT)
x['prio'] = x.gse.map(PRIORITY.index)
x = x.sort_values(['first', 'prio', 'title_clean']).reset_index(drop=True)
x = x.rename(columns={'run_accession': 'run', 'library_layout': 'layout'})
cols = ['gse', 'run', 'layout', 'fastq_ftp', 'fastq_bytes', 'fastq_md5']
for w in range(3):
    x.iloc[w::3][cols].to_csv(f'{B}/metadata/perturb/gate3_ubuntu_worker{w}.tsv', sep='\t', index=False)
    print(w, len(x.iloc[w::3]))
