#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
INPUT="/Users/darveshgorhe/columbia/phd/marko/structural_spidr_analysis/results/replicate_comparison/LAS_AK/bedgraphs/peaks/LAS/Ctr_LARP1.bed"
GTFDB="$ROOT/benchmark/gencode.v50.primary_assembly.annotation.gtf.db"
INDEX="${GTFDB}.annotator.idx"
OUT_DIR="$ROOT/benchmark"
VENV="$OUT_DIR/modern-venv"
ANNOTATOR="$VENV/bin/annotator"
PRIOR="$OUT_DIR/Ctr_LARP1.modern_c1.bed"

mkdir -p "$OUT_DIR"

echo "=== Building annotator index (one-time) ==="
if [[ ! -f "$INDEX" ]]; then
  /usr/bin/time -p "$VENV/bin/python" -c "
from annotator.index_store import build_annotator_index
build_annotator_index('${GTFDB}', force=True, gtf_format='gencode')
print('Index written to ${INDEX}')
" 2> "$OUT_DIR/build_index.time"
else
  echo "Index already exists: $INDEX"
fi

run_case() {
  local label="$1"
  shift
  local out="$OUT_DIR/Ctr_LARP1.${label}.bed"
  local meta="$OUT_DIR/Ctr_LARP1.${label}.meta.txt"

  echo "=== Running $label ==="
  if (($#)); then
    /usr/bin/time -p "$ANNOTATOR" \
      --input "$INPUT" \
      --output "$out" \
      --gtfdb "$GTFDB" \
      --cores 1 \
      "$@" \
      2> "$meta.time" || true
  else
    /usr/bin/time -p "$ANNOTATOR" \
      --input "$INPUT" \
      --output "$out" \
      --gtfdb "$GTFDB" \
      --cores 1 \
      2> "$meta.time" || true
  fi

  {
    echo "label=$label"
    echo "input_sha256=$(shasum -a 256 "$INPUT" | awk '{print $1}')"
    echo "output_sha256=$(shasum -a 256 "$out" | awk '{print $1}')"
    echo "row_count=$(wc -l < "$out" | tr -d ' ')"
    echo "gtfdb=$GTFDB"
    echo "index=$INDEX"
    cat "$meta.time"
  } > "$meta"
}

run_case "index_c1"
run_case "noindex_c1" --no-index

if [[ -f "$PRIOR" ]]; then
  diff -q "$OUT_DIR/Ctr_LARP1.index_c1.bed" "$PRIOR" > "$OUT_DIR/index_vs_prior_modern.diff" 2>&1 || true
fi
diff -q "$OUT_DIR/Ctr_LARP1.index_c1.bed" "$OUT_DIR/Ctr_LARP1.noindex_c1.bed" \
  > "$OUT_DIR/index_vs_noindex.diff" 2>&1 || true

echo "=== Summary ==="
for f in "$OUT_DIR"/Ctr_LARP1.{index_c1,noindex_c1}.meta.txt; do
  echo "--- $(basename "$f") ---"
  grep -E '^(label|output_sha256|row_count|real)' "$f" || true
done
echo "--- diffs ---"
cat "$OUT_DIR/index_vs_noindex.diff" 2>/dev/null || echo "(no diff file)"
cat "$OUT_DIR/index_vs_prior_modern.diff" 2>/dev/null || true
