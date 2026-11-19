"""Checks for final_h4m_state.py (S0 of the recompute matrix). No model runs.

1. The river-cell CSV reproduces out/river_cell_mask.npy (1,020 cells) and river_cells_main.npy (218).
2. build_state_h4m() with nominal parameters is identical to the Phase B4 H4m state
   (out/b4/runs/def_H4m.pkl and h_scenarios.make_state("H4m")): attractors, pairs, E', c',
   and every edge weight of the routing graph.
3. cost_mask="main" reproduces the B4 H4m-main state.
4. rules="v2" reproduces the attractor/pair counts of the published Tables 4 and 5 (design check).
Writes out/b6/fingerprint_H4m.json (used by the Ryzen environment test).
"""
import os as _os, sys as _sys, json
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import pickle, numpy as np
import final_h4m_state as fs
from h_scenarios import make_state

CSV = "data/grilla_con_W_hidrica.csv"
RIV = "data/celdas_fluviales_1020.csv"
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"[{'OK ' if cond else 'FAIL'}] {name} {detail}")


st_v2 = pickle.load(open("out/state_v2_repro.pkl", "rb"))
R0 = np.load("out/river_cell_mask.npy"); Rm0 = np.load("out/hydro/river_cells_main.npy")
if not _Path(RIV).exists():
    fs.write_river_cells_csv(st_v2, R0, Rm0, RIV)
R, Rm = fs.river_masks(st_v2, RIV)
check("river cells CSV = river_cell_mask.npy", np.array_equal(R, R0 & st_v2.valid_mask), f"({int(R.sum())} cells)")
check("main cells CSV = river_cells_main.npy", np.array_equal(Rm, Rm0 & st_v2.valid_mask), f"({int(Rm.sum())} cells)")
d_csv = fs.dist_grid(st_v2); d_npy = np.load("out/dist_river_grid.npy")
check("dist grid from CSV = dist_river_grid.npy", np.array_equal(np.nan_to_num(d_csv, nan=-1), np.nan_to_num(d_npy, nan=-1)))


def edge_diff(G1, G2):
    assert G1.number_of_edges() == G2.number_of_edges()
    return max(max(abs(G1.edges[e]["weight_base"] - G2.edges[e]["weight_base"]),
                   abs(G1.edges[e]["weight"] - G2.edges[e]["weight"])) for e in G1.edges)


for scen, kw in [("H4m", {}), ("H4m-main", {"cost_mask": "main"})]:
    st = fs.build_state_h4m(CSV, RIV, **kw)
    ref = make_state(scen)
    V = st.valid_mask
    check(f"{scen}: attractors identical to B4", st.attractors == ref.attractors, f"({len(st.attractors)})")
    check(f"{scen}: pairs identical to B4", st.pair_candidates == ref.pair_candidates, f"({len(st.pair_candidates)})")
    check(f"{scen}: E' identical", np.abs(st.E_grid - ref.E_grid)[V].max() == 0)
    check(f"{scen}: c' identical", np.abs(st.c_grid - ref.c_grid)[V].max() == 0, f"(max c' = {st.c_grid[V].max():.4f})")
    check(f"{scen}: T identical", np.abs(st.T_grid - ref.T_grid)[V].max() == 0)
    check(f"{scen}: graph edge weights identical", edge_diff(st.G_template, ref.G_template) == 0)
    if scen == "H4m":
        b4 = pickle.load(open("out/b4/runs/def_H4m.pkl", "rb"))
        check("H4m: same as def_H4m.pkl (B4 runs)", st.attractors == b4["attractors"] and st.pair_candidates == b4["pairs"]
              and np.abs(st.E_grid - b4["E"])[V].max() == 0 and np.abs(st.c_grid - b4["c"])[V].max() == 0)
        fp = fs.fingerprint(st)
        _Path("out/b6").mkdir(parents=True, exist_ok=True)
        json.dump(fp, open("out/b6/fingerprint_H4m.json", "w"), indent=1)
        print("   fingerprint:", fp["n_attractors"], fp["n_pairs"], fp["sha256_E"][:12], fp["sha256_c"][:12])

# Design check against the published v2 tables (rules="v2" = 125 m buffer, as in v8_11)
TAB4 = {"baseline": 25, "r=0.80": 24, "r=0.50": 25, "s=0.40": 25, "s=0.00": 24, "w=7": 25, "w=3": 24}
for lab, r, s, w in fs.T_configs():
    st = fs.build_state_h4m(CSV, RIV, T_r=r, T_s=s, T_window=w, rules="v2")
    check(f"v2 Table 4 {lab}: attractors", len(st.attractors) == TAB4[lab], f"({len(st.attractors)} vs {TAB4[lab]})")
TAB5 = {0.90: (26, 93), 0.92: (25, 91), 0.95: (25, 91), 0.975: (22, 78)}
for q in fs.QUANTILES:
    st = fs.build_state_h4m(CSV, RIV, peak_quantile=q, rules="v2")
    got = (len(st.attractors), len(st.pair_candidates))
    check(f"v2 Table 5 q={q}: attractors/pairs", got == TAB5[q], f"({got} vs {TAB5[q]})")
w = fs.weight_configs()
check("weight design: 17 configs, unit sums", len(w) == 17 and all(abs(sum(a) - 1) < 1e-12 and abs(sum(b) - 1) < 1e-12 for _, a, b in w))
print("ALL OK" if ok else "SOME CHECKS FAILED")
_sys.exit(0 if ok else 1)
