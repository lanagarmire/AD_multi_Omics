#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON_BIN:-python3}"
rscript_bin="${RSCRIPT_BIN:-Rscript}"
mkdir -p "$root/results" "$root/figures" "$root/panels" "$root/docs"
"$python_bin" "$root/code/preprocess_and_transfer.py"
"${TRAJECTORY_RSCRIPT:-$rscript_bin}" "$root/code/fit_trajectory.R" "$root"
"$rscript_bin" "$root/code/analyze_DE.R" "$root"
"$rscript_bin" "$root/code/pair_subset_gam.R" "$root"
"$python_bin" - "$root" <<'PYCODE'
from pathlib import Path
import sys, shutil
root=Path(sys.argv[1])
for group in ['EMCI','LMCI']:
    out=root/'results/metaboanalyst'/group
    out.mkdir(parents=True,exist_ok=True)
    for name in ['compound_db.qs','hsa.qs']:
        shutil.copy2(root/'resources/metaboanalyst'/group/name,out/name)
PYCODE
"${PATHWAY_RSCRIPT:-$rscript_bin}" "$root/code/run_metaboanalyst_pathways.R" "$root"
"$python_bin" "$root/code/export_full_metabolite_names.py"
"$python_bin" "$root/code/build_figure.py"
