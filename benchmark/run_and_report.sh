#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS_FILE="${1:-$ROOT_DIR/benchmark/results.csv}"
REPORT_FILE="${2:-$ROOT_DIR/benchmark/REPORT.md}"

"$ROOT_DIR/benchmark/run.sh" "$RESULTS_FILE"
python3 "$ROOT_DIR/benchmark/generate_report.py" "$RESULTS_FILE" "$REPORT_FILE"
