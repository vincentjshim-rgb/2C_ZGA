#!/usr/bin/env bash
# Gate 3 downloads + kallisto quant for the perturbation datasets. FASTQ are KEPT (Gate 4 splicing).
# Usage: bash 10_gate3_download_quant.sh <runlist.tsv>   (columns: gse, run, layout, fastq_ftp)
# Resumable (curl -C -, DONE markers). One worker per call; launch several calls on disjoint run lists.
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IDX="${IDX:-$B/../data/processed/mm_full.idx}"   # GRCm39 r112 cDNA+ncRNA kallisto index
KAL="${KALLISTO:-kallisto}"
LIST=$1
while IFS=$'\t' read -r gse run layout ftp; do
  [ "$gse" = "gse" ] && continue
  OUT=$B/data/perturb/$gse
  mkdir -p "$OUT/fastq" "$OUT/quant"
  LOG=$OUT/quant.log
  if [ -f "$OUT/quant/$run/DONE" ]; then continue; fi
  files=""
  ok=1
  for u in $(echo "$ftp" | tr ';' ' '); do
    f="$OUT/fastq/$(basename "$u")"
    if [ ! -f "$f.ok" ]; then
      curl -s -C - --retry 8 --retry-delay 20 -m 14400 -o "$f" "https://$u" && gzip -t "$f" && touch "$f.ok" || { echo "$run download/corrupt $u" >> "$LOG"; ok=0; }
    fi
    files="$files $f"
  done
  [ $ok -eq 1 ] || continue
  mkdir -p "$OUT/quant/$run"
  if [ "$layout" = "PAIRED" ]; then
    "$KAL" quant -i "$IDX" -o "$OUT/quant/$run" -t 4 $files > "$OUT/quant/$run/kallisto.log" 2>&1
  else
    "$KAL" quant -i "$IDX" -o "$OUT/quant/$run" -t 4 --single -l 200 -s 30 $files > "$OUT/quant/$run/kallisto.log" 2>&1
  fi
  if [ -s "$OUT/quant/$run/abundance.tsv" ]; then touch "$OUT/quant/$run/DONE"; echo "$run quant ok" >> "$LOG"; else echo "$run quant failed" >> "$LOG"; fi
done < <(tr -d '\r' < "$LIST")
