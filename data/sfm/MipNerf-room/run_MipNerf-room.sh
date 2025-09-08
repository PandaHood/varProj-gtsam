#!/usr/bin/env bash
set -euo pipefail

# Paths
REPO_ROOT="$(git rev-parse --show-toplevel)"

BIN="$REPO_ROOT/cmake-build-default/bin/BA_example"
DIM=3
DATA_DIR="$REPO_ROOT/data/sfm/MipNerf-room"
PYFG="$DATA_DIR/MipNerf-room.pyfg"
INITS_DIR="$DATA_DIR/inits"

# Checks
[[ -x "$BIN" ]] || { echo "ERROR: Binary not found/executable: $BIN"; exit 1; }
[[ -f "$PYFG" ]] || { echo "ERROR: .pyfg not found: $PYFG"; exit 1; }
[[ -d "$INITS_DIR" ]] || { echo "ERROR: inits/ not found: $INITS_DIR"; exit 1; }

shopt -s nullglob
for INIT_PATH in "$INITS_DIR"/rank*_init*.txt; do
  INIT_FILE="$(basename "$INIT_PATH")"

  # Expected pattern: rankP_initN.txt (e.g., rank6_init7.txt)
  if [[ "$INIT_FILE" =~ ^rank([0-9]+)_init([0-9]+)\.txt$ ]]; then
    RANK="${BASH_REMATCH[1]}"
    INIT_NUM="${BASH_REMATCH[2]}"

    echo "=== Running: DIM=$DIM  RANK=$RANK  INIT=$INIT_NUM ==="
    "$BIN" "$DIM" "$RANK" "$PYFG" "$INIT_PATH"
  else
    echo "Skipping unrecognized file name: $INIT_FILE" >&2
  fi
done

echo "All runs complete."