"""Hydrological scenarios H1..H4 (pre-registered in PRERREGISTRO_FASE_B.md).

Usage: python code/h_scenarios.py SCEN:SEED [SCEN:SEED ...]
Saves out/h/<scen>_s<seed>.npz (F, L, auc, lcc, m2) and out/h/def_<scen>.pkl
(attractors, pairs, E, c) for the route diagnostics.
"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import sys, time, pickle, numpy as np
from dataclasses import replace
import ca_kernel_v2 as ca
from scenario_eval import with_attractors

SCEN = {  # name: (delta for E, gamma for c)
    "H1": (0.0, 0.0), "H1+H2m": (0.25, 0.0), "H1+H2s": (0.50, 0.0),
    "H1+H3m": (0.0, 0.10), "H1+H3s": (0.0, 0.30), "H4m": (0.25, 0.10), "H4s": (0.50, 0.30),
    "H4m-main": (0.25, 0.10)}          # post hoc (ADENDA_REGISTRO_FASE_B2.md): cost penalty on main channel only
COST_MASK = {"H4m-main": "out/hydro/river_cells_main.npy"}   # default: all river cells
LAMBDA_M = 250.0


def make_state(name):
    st0 = pickle.load(open("out/state_v2_repro.pkl", "rb"))
    R = np.load("out/river_cell_mask.npy"); dist = np.load("out/dist_river_grid.npy")
    delta, gamma = SCEN[name]
    d0 = np.where(np.isnan(dist), 1e9, dist)
    E = st0.E_grid * (1.0 - delta * np.exp(-d0 / LAMBDA_M)) if delta > 0 else st0.E_grid
    st = with_attractors(st0, R | (d0 <= 125), E_grid=E)
    if gamma > 0:
        Rc = np.load(COST_MASK[name]) if name in COST_MASK else R
        c = st0.c_grid + gamma * Rc
        st = replace(st, c_grid=c, G_template=ca.build_graph(c, st0.valid_mask))
    assert not any(R[i, j] for i, j in st.attractors)
    return st


if __name__ == "__main__":
    cache = {}
    for job in sys.argv[1:]:
        name, seed = job.split(":"); seed = int(seed)
        if name not in cache:
            cache[name] = make_state(name)
            s = cache[name]
            pickle.dump(dict(attractors=s.attractors, pairs=s.pair_candidates, E=s.E_grid, c=s.c_grid),
                        open(f"out/h/def_{name}.pkl", "wb"))
        st = cache[name]
        t0 = time.time()
        r = ca.run_one(st, eps=0.15, kappa=0.03, seed=seed, return_fields=True)
        np.savez_compressed(f"out/h/{name}_s{seed}.npz", F=r["F_norm"], L=r["L"], auc=r["auc"],
                            lcc=r["lcc_size"], m2=r["m2"])
        print(f"{name} seed={seed} attr={len(st.attractors)} pairs={len(st.pair_candidates)} "
              f"AUC={r['auc']:.4f} ({time.time()-t0:.0f}s)", flush=True)
