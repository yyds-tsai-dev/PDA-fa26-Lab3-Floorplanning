#!/usr/bin/env bash
# Build the dev tuning set: 5 presets x seeds 1001-1004 -> cases/dev/ (20 cases).
# Graded cases use different seeds (9000s); never tune on them.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$ROOT/cases/dev"
for preset in public hidden1 hidden2 hidden3 hidden4; do
  for seed in 1001 1002 1003 1004; do
    python3 "$ROOT/generator/gen.py" --preset "$preset" --seed "$seed" \
      --out "$ROOT/cases/dev/${preset}_s${seed}"
  done
done
