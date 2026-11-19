"""Zone ranking of the fixation field L against the terrace suitability field T (canonical run S1, seed 42). No model runs.

For the reference K-means partitions with 20, 30 and 40 zones: Spearman correlation between the zone means of L and of
T, whether the two zone rankings are identical, and Spearman correlation between the zone means of L and of the traffic
field F. The last column is the Spearman correlation between L and T over all valid cells.
Output: out/final/informe/diagnostico_L_vs_T.csv

Added in release v2.0.0 so that this file has its own script in the public repository. The computation is copied
without changes from the report script used during the revision (b8_informe.py), whose other parts read review files
that are not part of this repository. Running this script reproduces the stored CSV (verify_repository.py, step 7)."""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import final_common as fc
import final_runs as fr

INF = _Path("out/final/informe"); INF.mkdir(parents=True, exist_ok=True)

st = fc.state_for(fr.cfg())
z = np.load("out/final/runs/S1/canon_s42.npz"); L, F = z["L"], z["F"]; V = st.valid_mask
rows = []
for n in (20, 30, 40):
    part = fc.partition(st, n)
    zl, zt, zf = (fc.zone_frame(x, part).L_mean for x in (L, st.T_grid, F))
    rows.append(dict(n=n, spearman_zone_L_T=spearmanr(zl, zt)[0], identical_zone_ranking=bool((zl.rank() == zt.rank()).all()),
                     spearman_zone_L_F=spearmanr(zl, zf)[0]))
LT = pd.DataFrame(rows); LT["spearman_cell_L_T"] = spearmanr(L[V], st.T_grid[V])[0]
LT.to_csv(INF / "diagnostico_L_vs_T.csv", index=False)
print(LT.to_string(index=False))
