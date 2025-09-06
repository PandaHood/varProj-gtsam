#!/usr/bin/env bash
set -uo pipefail

# Root containing the SNL dataset folders
# Override with: SNL_BASE=/custom/path ./run_snl_all.sh
BASE_DIR="${SNL_BASE:-/home/nikolas/StiefelManifold/data/snl}"

# Default folders. If you pass args, they'll replace this list (e.g., "./run_snl_all.sh MIT_snl sphere2500_snl").
FOLDERS=(intel_snl M3500_snl MIT_snl parking-garage_snl sphere2500_snl)
if [[ $# -gt 0 ]]; then
  FOLDERS=("$@")
fi

failures=()

for folder in "${FOLDERS[@]}"; do
  dir="$BASE_DIR/$folder"
  script="run_${folder}.sh"
  path="$dir/$script"

  if [[ ! -d "$dir" ]]; then
    echo "⚠️  Skipping '$folder': missing directory $dir" >&2
    failures+=("$folder (missing dir)")
    continue
  fi
  if [[ ! -f "$path" ]]; then
    echo "⚠️  Skipping '$folder': missing script $path" >&2
    failures+=("$folder (missing script)")
    continue
  fi

  # Ensure it's executable; if not, try to fix it
  [[ -x "$path" ]] || chmod +x "$path" 2>/dev/null || true

  echo "▶️  Running $script in $dir"
  ( cd "$dir" && "./$script" )
  rc=$?

  if [[ $rc -ne 0 ]]; then
    echo "❌ $script failed with exit code $rc"
    failures+=("$folder (exit $rc)")
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

echo "🎉 All SNL scripts completed successfully."
