#!/usr/bin/env bash
# Gate 2b dataset R2 (GSE66582, Wu 2016): download stage RNA-seq runs from ENA and quantify with kallisto
# against the project GRCm39 r112 cDNA+ncRNA index. FASTQ deleted after a successful quant. Resumable.
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IDX="${IDX:-$B/../data/processed/mm_full.idx}"   # GRCm39 r112 cDNA+ncRNA kallisto index
OUT=$B/data/GSE66582
mkdir -p "$OUT/fastq" "$OUT/quant"
LOG=$OUT/quant.log
RUNS="SRR2927026 SRR2927027 SRR1840514 SRR1840515 SRR2927028 SRR2927029 SRR1840516 SRR1840517 SRR1840518 SRR1840519 SRR1840520 SRR1840521 SRR1840522 SRR1840523 SRR1840526"
curl -s -m 120 "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP055882&result=read_run&fields=run_accession,library_layout,fastq_ftp&format=tsv" | tr -d '\r' > "$OUT/ena_fastq.tsv"
for r in $RUNS; do
  if [ -f "$OUT/quant/$r/DONE" ]; then echo "$r done" >> "$LOG"; continue; fi
  line=$(awk -F"	" -v r="$r" '$1==r' "$OUT/ena_fastq.tsv")
  layout=$(echo "$line" | cut -f2)
  urls=$(echo "$line" | cut -f3 | tr ';' ' ')
  files=""
  for u in $urls; do
    f="$OUT/fastq/$(basename "$u")"
    curl -s -C - --retry 5 -m 7200 -o "$f" "https://$u" || { echo "$r download failed $u" >> "$LOG"; continue 2; }
    gzip -t "$f" || { echo "$r corrupt $f" >> "$LOG"; rm -f "$f"; continue 2; }
    files="$files $f"
  done
  mkdir -p "$OUT/quant/$r"
  if [ "$layout" = "PAIRED" ]; then
    "${KALLISTO:-kallisto}" quant -i "$IDX" -o "$OUT/quant/$r" -t 6 $files >> "$LOG" 2>&1
  else
    "${KALLISTO:-kallisto}" quant -i "$IDX" -o "$OUT/quant/$r" -t 6 --single -l 200 -s 30 $files >> "$LOG" 2>&1
  fi
  if [ -s "$OUT/quant/$r/abundance.tsv" ]; then
    touch "$OUT/quant/$r/DONE"; rm -f $files; echo "$r quant ok ($layout)" >> "$LOG"
  else
    echo "$r quant failed" >> "$LOG"
  fi
done
touch "$OUT/ALL_DONE"
