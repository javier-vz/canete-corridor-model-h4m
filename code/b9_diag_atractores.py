"""JCAA-80, Phase B9 (report layer; no simulation, no change to any frozen output).

1. Attractor sets of the sensitivity configurations (Tables 3-5): the deterministic states are rebuilt with
   the frozen code (final_common.state_for) and each one is checked against the state_hash recorded in its
   S3/S4/S5 run JSON, so the counts refer to exactly the states used in the runs.
   Output: out/final/informe/atractores_sensibilidades.csv
2. Attractor count near the 8 priority lines (correction note to the signed checklist): per attractor of S0,
   distance to the 8 lines (same geometry and FIDs as final_evaluate.py) and whether it disappears in
   H4m-excl8 (frozen file out/b5/t5/attractors_changed.csv).
   Output: out/final/informe/conteo_atractores_8_lineas.csv"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import json
import numpy as np
import pandas as pd
import final_common as fc
import final_runs as fr
import final_h4m_state as fs

OUT = _Path("out/final/informe"); OUT.mkdir(parents=True, exist_ok=True)
base = fc.state_for(fr.cfg())
B = set(map(tuple, base.attractors))

# ---------------------------------------------------------------- 1. sensitivity states
rows = []
fam = [("Table 3 (weights)", "S3", 3000, lab, fr.cfg(a=list(a), b=list(b))) for lab, a, b in fs.weight_configs()]
fam += [("Table 4 (T)", "S4", 5000, lab, fr.cfg(T=[r, s, w])) for lab, r, s, w in fs.T_configs()]
fam += [("Table 5 (quantile)", "S5", 4000, f"q={q}", fr.cfg(q=q)) for q in fs.QUANTILES]
for table, block, seed, lab, conf in fam:
    st = fc.state_for(conf)
    jid = f"q{lab[2:]}_s{seed}" if block == "S5" else f"{fr.safe(lab)}_s{seed}"
    j = json.load(open(f"out/final/runs/{block}/{jid}.json"))
    h = fc.state_key_hash(st)
    assert h == j["state_hash"], (block, jid)
    A = set(map(tuple, st.attractors))
    rows.append(dict(table=table, config=lab, run_json=f"{block}/{jid}.json", state_hash_matches_run=True,
                     attractors=len(A), pairs=len(st.pair_candidates), shared_with_S0=len(A & B),
                     only_in_config=len(A - B), only_in_S0=len(B - A)))
t = pd.DataFrame(rows)
t.to_csv(OUT / "atractores_sensibilidades.csv", index=False)
print(t.to_string(index=False))

# ---------------------------------------------------------------- 2. 8 priority lines
import geopandas as gpd
from shapely.geometry import Point
cl = gpd.read_file("out/hydro/hydrorivers_model_lines_classified.gpkg")                 # as in final_evaluate.py
fid = pd.read_csv("out/b4/t1/hydrorivers_official_match.csv").query("dataset == 'C87'").set_index("line").shp_fid
cl["shp_fid"] = [int(fid[i]) for i in range(len(cl))]
EXCL_FID = [1131506, 1130585, 1129684, 1130586, 1127717, 1130415, 1126308, 1126640]   # = b5_exclusion_runs.EXCL_FID
L8 = cl[cl.shp_fid.isin(EXCL_FID)].geometry.union_all()
ix2x = {v: k for k, v in base.lattice["x_to_ix"].items()}; iy2y = {v: k for k, v in base.lattice["y_to_iy"].items()}
ch = pd.read_csv("out/b5/t5/attractors_changed.csv")
gone = set(zip(ch[ch.change == "only in H4m"].iy, ch[ch.change == "only in H4m"].ix))
new = ch[ch.change == "only in H4m-excl8"]
rows = []
for i, j in base.attractors:
    d = Point(ix2x[j], iy2y[i]).distance(L8)
    rows.append(dict(set="S0 (H4m)", iy=i, ix=j, X=ix2x[j], Y=iy2y[i], dist_to_8_lines_m=round(d, 1),
                     within_1km=d <= 1000, disappears_in_excl8=(i, j) in gone))
for r in new.itertuples():
    rows.append(dict(set="new in H4m-excl8", iy=r.iy, ix=r.ix, X=r.X, Y=r.Y, dist_to_8_lines_m=float(r.dist_to_8_lines_m),
                     within_1km=r.dist_to_8_lines_m <= 1000, disappears_in_excl8=False))
c = pd.DataFrame(rows)
c.to_csv(OUT / "conteo_atractores_8_lineas.csv", index=False)
s0 = c[c.set == "S0 (H4m)"]
t11 = pd.read_csv("out/final/eval/tabla11_hidrologia.csv").set_index("scenario")
summ = dict(S0_attractors=len(s0), within_1km=int(s0.within_1km.sum()),
            within_1km_min_m=s0[s0.within_1km].dist_to_8_lines_m.min(), within_1km_max_m=s0[s0.within_1km].dist_to_8_lines_m.max(),
            within_1km_and_disappear=int((s0.within_1km & s0.disappears_in_excl8).sum()),
            disappear_total=int(s0.disappears_in_excl8.sum()), appear=len(new),
            excl8_attractors=len(s0) - int(s0.disappears_in_excl8.sum()) + len(new),
            excl8_attractors_table11=float(t11.loc["H4m-excl8", "attractors_mean (B5)"]),
            excl8_pairs_table11=float(t11.loc["H4m-excl8", "pairs_mean (B5)"]))
assert summ["excl8_attractors"] == summ["excl8_attractors_table11"]
print(summ)
