"""Rebuild the transcript-to-gene map from BOTH Ensembl r112 FASTA files.

The original metadata/tx2gene.tsv was built from the cDNA FASTA only, so the
29,930 ncRNA transcripts that are in the kallisto index (mm_full.idx, built from
cDNA + ncRNA) were silently dropped at gene aggregation. This writes
metadata/tx2gene_full.tsv and leaves the old file untouched.
"""
import gzip, re, pandas as pd

rows = []
for fa, src in [('data/raw/mm_cdna.fa.gz', 'cdna'), ('data/raw/mm_ncrna.fa.gz', 'ncrna')]:
    with gzip.open(fa, 'rt') as f:
        for line in f:
            if not line.startswith('>'):
                continue
            h = line[1:].rstrip('\n')
            tx = h.split()[0]
            gid = re.search(r'gene:(\S+)', h)
            sym = re.search(r'gene_symbol:(\S+)', h)
            gbt = re.search(r'gene_biotype:(\S+)', h)
            tbt = re.search(r'transcript_biotype:(\S+)', h)
            loc = re.search(r'chromosome:[^:]+:([^:]+):(\d+):(\d+)', h)
            gene_id = gid.group(1).split('.')[0] if gid else tx
            rows.append(dict(tx=tx, gene_id=gene_id,
                             symbol=sym.group(1) if sym else gene_id,
                             gene_biotype=gbt.group(1) if gbt else 'NA',
                             transcript_biotype=tbt.group(1) if tbt else 'NA',
                             chrom=loc.group(1) if loc else 'NA', source=src))
T = pd.DataFrame(rows)
T.to_csv('metadata/tx2gene_full.tsv', sep='\t', index=False)
old = pd.read_csv('metadata/tx2gene.tsv', sep='\t')
print('transcripts: %d (cdna %d, ncrna %d)' % (len(T), (T.source == 'cdna').sum(), (T.source == 'ncrna').sum()))
print('genes: %d; old map genes: %d' % (T.gene_id.nunique(), old.gene_id.nunique()))
print('old map transcripts found in new map: %d / %d' % (old.tx.isin(T.tx).sum(), len(old)))
print('symbols mapping to >1 gene_id: %d' % (T.groupby('symbol').gene_id.nunique() > 1).sum())
print(T[T.source == 'ncrna'].gene_biotype.value_counts().head(10).to_string())
for g in ['Xist', 'Tsix', 'H19', 'Meg3', 'Airn', 'Kcnq1ot1', 'Gm29013']:
    print('  %-9s in map: %s (%s)' % (g, g in set(T.symbol), ','.join(sorted(set(T.source[T.symbol == g])))))
