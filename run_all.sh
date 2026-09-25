#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--quick" ) ]]; then
  echo "Usage: ./run_all.sh [--quick]" >&2
  exit 2
fi

ROOT="$(cd "$(dirname "$0")" && pwd)"
PYTHON="${PYTHON:-/opt/homebrew/bin/python3.11}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

step() {
  local label="$1"
  shift
  local started="$SECONDS"
  if ! "$@" >"$WORK/$label.log" 2>&1; then
    echo "FAIL $label" >&2
    tail -n 8 "$WORK/$label.log" >&2
    exit 1
  fi
  echo "PASS $label ($((SECONDS-started))s)"
}

step catalog "$PYTHON" -B "$ROOT/scripts/round9_matroid_backsolve/run.py" \
  --stage catalog --output "$WORK/catalog.json"
step family "$PYTHON" -B "$ROOT/scripts/round9_matroid_backsolve/run.py" \
  --stage family --output "$WORK/family.json"
if [[ "${1:-}" != "--quick" ]]; then
  step matroid_graphs "$PYTHON" -B "$ROOT/scripts/round9_matroid_backsolve/run.py" \
    --stage graphs --max-vertices 8 --max-edges 14 --output "$WORK/matroid_graphs.json"
fi
step exchange "$PYTHON" -B "$ROOT/scripts/round11_h1_proof/verify_exchange.py" \
  --output "$WORK/exchange.json"
step kernel "$PYTHON" -B "$ROOT/scripts/round12_kernel_review/verify_kernel.py" \
  --output "$WORK/kernel" --seconds 600
step cover "$PYTHON" -B "$ROOT/scripts/round8_review/t21_check.py" \
  --output "$WORK/cover.json"
step graph_small "$PYTHON" -B "$ROOT/scripts/round7/cover_census.py" small
step paper_numbers "$PYTHON" -B "$ROOT/scripts/verify_summaries.py" "$WORK"
tail -n 4 "$WORK/paper_numbers.log"
if [[ "${1:-}" == "--quick" ]]; then
  echo "PASS quick suite ($SECONDS s)"
else
  echo "PASS standard suite ($SECONDS s)"
fi
