#!/usr/bin/env bash
set -uo pipefail

# Base directory for MRCLAM datasets
# You can override with: MRCLAM_BASE=/custom/path ./run_mrclam.sh
BASE_DIR="${MRCLAM_BASE:-${RASLAM_BASE:-/home/nikolas/varProj-gtsam/data/raslam}}/mrclam"

# Default targets; can be overridden by CLI args
DEFAULT_TARGETS=(mrclam2 mrclam4 mrclam6 mrclam7)

# Build target list from args (accept "2 4" or "mrclam2 mrclam4")
TARGETS=()
if [[ $# -gt 0 ]]; then
  for arg in "$@"; do
    if [[ "$arg" =~ ^[0-9]+$ ]]; then
      TARGETS+=("mrclam$arg")
    else
      TARGETS+=("$arg")
    fi
  done
else
  TARGETS=("${DEFAULT_TARGETS[@]}")
fi

failures=()

for sub in "${TARGETS[@]}"; do
  dir="$BASE_DIR/$sub"
  script="run_${sub}.sh"
  path="$dir/$script"

  if [[ ! -d "$dir" ]]; then
    echo "⚠️  Skipping '$sub': missing directory $dir" >&2
    failures+=("$sub (missing dir)")
    continue
  fi

  if [[ ! -f "$path" ]]; then
    echo "⚠️  Skipping '$sub': missing script $path" >&2
    failures+=("$sub (missing script)")
    continue
  fi

  # Ensure executable (best-effort)
  [[ -x "$path" ]] || chmod +x "$path" 2>/dev/null || true

  echo "▶️  Running $script in $dir"
  ( cd "$dir" && "./$script" )
  rc=$?

  if [[ $rc -ne 0 ]]; then
    echo "❌ $script failed with exit code $rc"
    failures+=("$sub (exit $rc)")
  else
    echo "✅ $script completed"
  fi
  echo
done

if [[ ${#failures[@]} -gt 0 ]]; then
  echo "Finished with failures:"
  printf ' - %s\n' "${failures[@]}"
  exit 1
fi

echo "🎉 All MRCLAM scripts completed successfully."
