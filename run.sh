#!/bin/bash
# run.sh — Invoke ST3R scripts using the correct conda environment.
#
# Usage: ./run.sh <script_name.py> [args...]
# Example: ./run.sh 01_gef2h5ad.py -C rawData/sample.tsv -BT cell_bins -O Output/01
#
# Why this wrapper? The shipped scripts use hardcoded shebangs pointing to a
# placeholder path (e.g. /path/to/envs/Stereopy/bin/python). This wrapper
# dispatches each script to its matching conda environment via `conda run`,
# so you do not need to edit any shebang.

set -e

if [ $# -lt 1 ]; then
    echo "Usage: ./run.sh <script_name.py> [args...]"
    echo "Example: ./run.sh 01_gef2h5ad.py -C rawData/sample.tsv -BT cell_bins -O Output/01"
    exit 1
fi

SCRIPT="$1"
shift

# Map script name → conda environment
case "$(basename "$SCRIPT" .py)" in
    01_gef2h5ad)
        ENV=stereopy
        ;;
    02_concat|03_preprocess|04_squidpy|05_dataConvert|11_report|11_report_subset)
        ENV=scanpy
        ;;
    06_align|07_tdr|08_backbone|09_morph|10_interpolation)
        ENV=spateo_env
        ;;
    *)
        echo "[run.sh] Unknown script: $SCRIPT"
        echo "Valid scripts: 01_gef2h5ad, 02_concat, 03_preprocess, 04_squidpy, 05_dataConvert,"
        echo "                06_align, 07_tdr, 08_backbone, 09_morph, 10_interpolation,"
        echo "                11_report, 11_report_subset"
        exit 1
        ;;
esac

# Invoke script via conda run, with the repo's Scripts/ dir on PYTHONPATH
# (needed for Step 11 to find Report_config/lib/* modules).
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec conda run -n "$ENV" \
    env PYTHONPATH="$SCRIPT_DIR/Scripts/Report_config/lib:${PYTHONPATH:-}" \
    python "$SCRIPT_DIR/Scripts/$SCRIPT" "$@"