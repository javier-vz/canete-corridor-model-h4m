"""Attractor/pair counts of every sensitivity configuration under H4m (state construction only,
no model runs). Output: out/b6/design_counts_H4m.csv"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import pandas as pd
import final_h4m_state as fs
CSV, RIV = "data/grilla_con_W_hidrica.csv", "data/celdas_fluviales_1020.csv"
rows = []
for lab, a, b in fs.weight_configs():
    st = fs.build_state_h4m(CSV, RIV, a_weights=a, b_weights=b)
    rows.append(dict(family="weights", config=lab, attractors=len(st.attractors), pairs=len(st.pair_candidates)))
for lab, r, s, w in fs.T_configs():
    st = fs.build_state_h4m(CSV, RIV, T_r=r, T_s=s, T_window=w)
    rows.append(dict(family="T", config=lab, attractors=len(st.attractors), pairs=len(st.pair_candidates)))
for q in fs.QUANTILES:
    st = fs.build_state_h4m(CSV, RIV, peak_quantile=q)
    rows.append(dict(family="quantile", config=f"q={q}", attractors=len(st.attractors), pairs=len(st.pair_candidates)))
df = pd.DataFrame(rows); df.to_csv("out/b6/design_counts_H4m.csv", index=False); print(df.to_string(index=False))
