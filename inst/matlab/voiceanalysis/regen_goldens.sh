#!/usr/bin/env bash
# Regenerate MATLAB golden reference outputs and convert to RDS.
# Run from voiceanalysis/ directory.
set -euo pipefail

MATLAB=/Volumes/Fredrik_Nylen_0705853304/Applications/MATLAB_R2025b.app/bin/matlab

echo "[1/2] Running MATLAB pipeline on a1.wav..."
"$MATLAB" -nosplash -nodesktop -nojvm -noawt \
  -r "run('inst/matlab/gen_goldens.m'); exit" 2>&1 | tail -30

echo "[2/2] Converting .mat -> .rds..."
Rscript tools/mat_to_rds.R

echo "Done. Goldens at tests/testthat/golden/a1.rds and inst/extdata/creak_ann.rds"
