#!/usr/bin/env bash
# Gate 3 downloads + kallisto quant, Ubuntu version of 10_gate3_download_quant.sh (kallisto parameters unchanged).
# Changes vs 10_ (audit log 2026-09-15): each FASTQ verified against ENA size + md5 (not only gzip -t); a file that
# reaches its declared size but fails md5 is deleted and fetched again; interrupted files resume (curl -C -) in up to
# PASSES passes over the list; index = Ubuntu rebuild of mm_full.idx. FASTQ are KEPT (Gate 4).
# Usage: bash 10u_gate3_download_quant.sh <worklist.tsv> [passes]
#        (columns: gse, run, layout, fastq_ftp, fastq_bytes, fastq_md5; from src/10u_make_worklists.py)
set -u
B="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IDX="${IDX:-$B/../data/processed/mm_full.idx}"   # GRCm39 r112 cDNA+ncRNA kallisto index
KAL="${KALLISTO:-kallisto}"
LIST=$1
PASSES=${2:-30}
for pass in $(seq 1 "$PASSES"); do
  pending=0
  while IFS=$'\t' read -r gse run layout ftp bytes md5s; do
    [ "$gse" = "gse" ] && continue
    OUT=$B/data/perturb/$gse
    mkdir -p "$OUT/fastq" "$OUT/quant"
    LOG=$OUT/quant.log
    [ -f "$OUT/quant/$run/DONE" ] && continue
    IFS=';' read -r -a urls <<< "$ftp"
    IFS=';' read -r -a sizes <<< "$bytes"
    IFS=';' read -r -a sums <<< "$md5s"
    files=()
    ok=1
    for i in "${!urls[@]}"; do
      f="$OUT/fastq/$(basename "${urls[$i]}")"
      files+=("$f")
      [ -f "$f.ok" ] && continue
      have=$(stat -c %s "$f" 2>/dev/null || echo 0)
      if [ "$have" -lt "${sizes[$i]}" ]; then
        curl -s -C - --speed-limit 2000 --speed-time 600 -o "$f" "https://${urls[$i]}"
        have=$(stat -c %s "$f" 2>/dev/null || echo 0)
      fi
      if [ "$have" -lt "${sizes[$i]}" ]; then
        echo "$(date '+%F %T') $run incomplete $(basename "$f") $have/${sizes[$i]} bytes (pass $pass)" >> "$LOG"
        ok=0
        continue
      fi
      got=$(md5sum "$f" | cut -d' ' -f1)
      if [ "$got" = "${sums[$i]}" ]; then
        echo "$got" > "$f.ok"
        echo "$(date '+%F %T') $run md5 ok $(basename "$f") $have bytes" >> "$LOG"
      else
        echo "$(date '+%F %T') $run md5 MISMATCH $(basename "$f") got $got expected ${sums[$i]} ($have bytes); deleted" >> "$LOG"
        rm -f "$f"
        ok=0
      fi
    done
    if [ $ok -eq 0 ]; then pending=1; sleep 30; continue; fi
    mkdir -p "$OUT/quant/$run"
    if [ "$layout" = "PAIRED" ]; then
      "$KAL" quant -i "$IDX" -o "$OUT/quant/$run" -t 4 "${files[@]}" > "$OUT/quant/$run/kallisto.log" 2>&1
    else
      "$KAL" quant -i "$IDX" -o "$OUT/quant/$run" -t 4 --single -l 200 -s 30 "${files[@]}" > "$OUT/quant/$run/kallisto.log" 2>&1
    fi
    if [ -s "$OUT/quant/$run/abundance.tsv" ]; then
      touch "$OUT/quant/$run/DONE"; echo "$(date '+%F %T') $run quant ok" >> "$LOG"
    else
      echo "$(date '+%F %T') $run quant failed" >> "$LOG"; pending=1
    fi
  done < <(tr -d '\r' < "$LIST")
  [ $pending -eq 0 ] && { echo "$(date '+%F %T') $LIST all runs DONE after pass $pass"; exit 0; }
  echo "$(date '+%F %T') $LIST pass $pass finished with pending runs"
done
echo "$(date '+%F %T') $LIST stopped after $PASSES passes with pending runs"
exit 1
