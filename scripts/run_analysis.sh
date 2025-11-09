#!/usr/bin/env bash
export MPLBACKEND=Agg
# scripts/run_analysis.sh — robust batch runner for bioseq.analyzer
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 <fasta1> [fasta2 ...]" >&2
  exit 2
fi

OUTDIR="results"
mkdir -p "$OUTDIR"

echo "🚀 BioSeq Toolkit batch run: $(date)"

# enable nullglob so patterns that have no matches expand to nothing
shopt -s nullglob

for pat in "$@"; do
  # If the caller passed a pattern that contains a wildcard, expand it here.
  files=()
  if [[ "$pat" == *\** ]]; then
    # unquoted expansion so shell will expand wildcards
    for f in $pat; do
      files+=("$f")
    done
  else
    # treat as a literal filename
    files+=("$pat")
  fi

  # If no files found yet, try common FASTA extensions for the given base
  if [ "${#files[@]}" -eq 0 ]; then
    base="${pat%.*}"
    for ext in fasta fa fna fas; do
      for f in "${base}.${ext}"; do
        [ -e "$f" ] && files+=("$f") || true
      done
    done
  fi

  # If still empty, warn and continue
  if [ "${#files[@]}" -eq 0 ]; then
    echo "⚠️  Skipping (not found): $pat"
    continue
  fi

  # Process each resolved file
  for fasta in "${files[@]}"; do
    if [ ! -f "$fasta" ]; then
      echo "⚠️  Skipping (not a file): $fasta"
      continue
    fi

    base=$(basename "$fasta")
    name="${base%.*}"
    csv="$OUTDIR/${name}_summary.csv"
    log="$OUTDIR/${name}.log"

    echo "🔹 Processing $fasta → $csv (log: $log)"
    python -m bioseq.analyzer "$fasta" \
      -o "$csv" \
      --orf \
      --plot-stats \
      --logfile "$log" \
      --quiet

    echo "✅ Finished $fasta"
  done
done

# turn off nullglob to avoid side effects in interactive shells sourced after this
shopt -u nullglob

echo "🎉 Batch complete. Results in: $OUTDIR"

