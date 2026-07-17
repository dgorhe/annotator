#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
INPUT="/Users/darveshgorhe/columbia/phd/marko/structural_spidr_analysis/results/replicate_comparison/LAS_AK/bedgraphs/peaks/LAS/Ctr_LARP1.bed"
GTFDB="$ROOT/benchmark/gencode.v50.primary_assembly.annotation.gtf.db"
OUT_DIR="$ROOT/benchmark"
LEGACY_VENV="$OUT_DIR/baseline-venv"
MODERN_VENV="$OUT_DIR/modern-venv"

mkdir -p "$OUT_DIR"

run_case() {
  local label="$1"
  local venv="$2"
  local cores="$3"
  local out="$OUT_DIR/Ctr_LARP1.${label}.bed"
  local meta="$OUT_DIR/Ctr_LARP1.${label}.meta.txt"
  local -a extra_args=()
  if [[ "$label" != legacy_* ]]; then
    extra_args=(--cores "$cores")
  fi

  echo "Running $label"
  if ((${#extra_args[@]})); then
    /usr/bin/time -p "$venv/bin/annotator" \
      --input "$INPUT" \
      --output "$out" \
      --gtfdb "$GTFDB" \
      "${extra_args[@]}" \
      2> "$meta.time" || true
  else
    /usr/bin/time -p "$venv/bin/annotator" \
      --input "$INPUT" \
      --output "$out" \
      --gtfdb "$GTFDB" \
      2> "$meta.time" || true
  fi

  {
    echo "label=$label"
    echo "cores=$cores"
    echo "venv=$venv"
    echo "input_sha256=$(shasum -a 256 "$INPUT" | awk '{print $1}')"
    echo "output_sha256=$(shasum -a 256 "$out" | awk '{print $1}')"
    echo "row_count=$(wc -l < "$out" | tr -d ' ')"
    echo "gtfdb=$GTFDB"
    echo "gtfdb_bytes=$(stat -f%z "$GTFDB")"
    cat "$meta.time"
  } > "$meta"
}

if [[ ! -f "$GTFDB" ]]; then
  echo "Missing database: $GTFDB" >&2
  exit 1
fi

run_case "legacy_c1" "$LEGACY_VENV" 1
run_case "modern_c1" "$MODERN_VENV" 1
run_case "modern_c4" "$MODERN_VENV" 4

diff -q "$OUT_DIR/Ctr_LARP1.legacy_c1.bed" "$OUT_DIR/Ctr_LARP1.modern_c1.bed" \
  > "$OUT_DIR/legacy_vs_modern_c1.diff" 2>&1 || true
diff -q "$OUT_DIR/Ctr_LARP1.modern_c1.bed" "$OUT_DIR/Ctr_LARP1.modern_c4.bed" \
  > "$OUT_DIR/modern_c1_vs_c4.diff" 2>&1 || true

echo "Benchmark artifacts written to $OUT_DIR"
