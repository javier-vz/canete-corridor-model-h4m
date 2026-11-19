"""JCAA-80, Phase B9 (report layer; reads the frozen S0 and v2 states, no runs): affinity at attractors
and at site cells, before (v2: E, 125 m buffer) and after (H4m: E') - Reviewer 2, comment on E.
Also distance of attractors and site cells to the nearest HydroRIVERS reach (dist_river_m of the grid).
Outputs: out/final/informe/diagnostico_E_atractores_sitios.csv, rangos_campos_figura4.csv"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import pickle
import numpy as np
import pandas as pd
import final_common as fc
import final_runs as fr

st = fc.state_for(fr.cfg())                                   # H4m (S0)
v2 = pickle.load(open("out/state_v2_repro.pkl", "rb"))        # published v2 state
_, _, d0 = fc.masks_and_dist(st)
S = st.site_mask
rows = []
for model, field, attrs in (("v2 (E, 125 m buffer)", v2.E_grid, v2.attractors), ("H4m (E')", st.E_grid, st.attractors)):
    a = np.array([field[i, j] for i, j in attrs]); s = field[S]
    da = np.array([d0[i, j] for i, j in attrs]); ds = d0[S]
    for grp, x, d in (("attractors", a, da), ("site cells", s, ds)):
        rows.append(dict(model=model, group=grp, n=len(x), affinity_mean=x.mean(), affinity_median=np.median(x),
                         affinity_min=x.min(), affinity_max=x.max(), dist_river_median_m=np.median(d),
                         dist_river_min_m=d.min(), dist_river_max_m=d.max()))
out = pd.DataFrame(rows)
_Path("out/final/informe").mkdir(parents=True, exist_ok=True)
out.to_csv("out/final/informe/diagnostico_E_atractores_sitios.csv", index=False)
print(out.round(3).to_string(index=False))

# value ranges of the fields drawn in Figure 4 (caption states a common [0, 1] display range)
V = st.valid_mask
rng =pd.DataFrame([dict(field=k, min=float(a[V].min()), max=float(a[V].max())) for k, a in
                    (("E' (H4m)", st.E_grid), ("c' (H4m)", st.c_grid), ("W", st.w_grid), ("T", st.T_grid))])
rng.to_csv("out/final/informe/rangos_campos_figura4.csv", index=False)
print(rng.to_string(index=False))
