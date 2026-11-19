"""JCAA-80 final recompute: diagnostics for the response letter that need routing but no kernel runs.
D1  E2-A under the final configuration: the 90 H4m pairs routed on c' at 250 m and 50 m, 8 and 16
    directions (same method as e2a_constant_cost.py; c' replicated 5x5 to 50 m).
D9  pairs with long detours: sinuosity of the deterministic least-cost path (A1 = kernel graph) and
    each pair's share of the canonical seed-42 traffic (from the S1 path records).
Outputs: out/final/eval/resolucion_E2A_pares.csv, resolucion_E2A_resumen.csv, pares_rutas.csv"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import time
import numpy as np
import pandas as pd
from shapely import points as shp_points

import final_common as fc
import final_runs as fr
from routing import build_graph, route, polyline, sample_pts, mean_sep, angular
from river_geometry import load_rivers


def river_metrics(line, union):
    pts = sample_pts(line, 25.0)
    d = union.distance(shp_points([(p.x, p.y) for p in pts]))
    near = d <= 125.0
    runs, cur = [], 0
    for flag in near:
        if flag:
            cur += 1
        else:
            if cur:
                runs.append(cur); cur = 0
    if cur:
        runs.append(cur)
    step = line.length / max(len(pts) - 1, 1)
    follow = sum(r for r in runs if r * step >= 750) / len(near)
    inter = line.intersection(union)
    ncross = 0 if inter.is_empty else (len(inter.geoms) if hasattr(inter, "geoms") else 1)
    return dict(frac_within_125m=float(near.mean()), frac_following_ge750m=float(follow),
                river_crossings=ncross, crossings_per_10km=ncross / (line.length / 1e4))


def e2a(st):
    union = load_rivers(st.df).geometry.union_all()
    xs250, ys250 = np.array(st.lattice["xs"]), np.array(st.lattice["ys"])
    xs50 = 342700 + (np.arange(144 * 5) + 0.5) * 50; ys50 = 8547500 + (np.arange(93 * 5) + 0.5) * 50
    c50 = np.kron(st.c_grid, np.ones((5, 5))); v50 = np.kron(st.valid_mask, np.ones((5, 5), dtype=bool))
    CFG = {"A1": (st.c_grid, st.valid_mask, 250.0, False, 1), "A4": (st.c_grid, st.valid_mask, 250.0, True, 1),
           "A2": (c50, v50, 50.0, False, 5), "A3": (c50, v50, 50.0, True, 5)}
    lines, rows = {}, []
    for name, (c, V, res, kn, f) in CFG.items():
        t0 = time.time(); G, idx = build_graph(c, V, res, kn)
        xs, ys = (xs250, ys250) if f == 1 else (xs50, ys50)
        for k, (a, b) in enumerate(st.pair_candidates):
            (i1, j1), (i2, j2) = st.attractors[a], st.attractors[b]
            s, t = (i1 * f + f // 2, j1 * f + f // 2), (i2 * f + f // 2, j2 * f + f // 2)
            cells, cost = route(G, idx, s, t)
            ln = polyline(cells, xs, ys); lines[(name, k)] = ln
            straight = np.hypot(xs250[j2] - xs250[j1], ys250[i2] - ys250[i1])
            rows.append(dict(config=name, pair_idx=k, attr_a=a, attr_b=b, length_m=ln.length, straight_m=straight,
                             sinuosity=ln.length / straight, cost_map_units=cost, frac_off_45deg=1 - angular(cells),
                             **river_metrics(ln, union)))
        print(f"E2-A {name}: {time.time() - t0:.0f}s", flush=True)
    D = pd.DataFrame(rows)
    D["mean_sep_vs_A3_m"] = [mean_sep(lines[(r.config, r.pair_idx)], lines[("A3", r.pair_idx)]) for r in D.itertuples()]
    D["hausdorff_vs_A3_m"] = [lines[(r.config, r.pair_idx)].hausdorff_distance(lines[("A3", r.pair_idx)]) for r in D.itertuples()]
    S = D.groupby("config").agg(mean_sep_median=("mean_sep_vs_A3_m", "median"), mean_sep_p90=("mean_sep_vs_A3_m", lambda x: x.quantile(.9)),
                                hausdorff_median=("hausdorff_vs_A3_m", "median"), pairs_sep_lt_250=("mean_sep_vs_A3_m", lambda x: int((x < 250).sum())),
                                n_pairs=("pair_idx", "nunique"), sinuosity_median=("sinuosity", "median"),
                                within125_median=("frac_within_125m", "median"), following_median=("frac_following_ge750m", "median"),
                                crossings_per_10km=("crossings_per_10km", "median"), off45_median=("frac_off_45deg", "median")).reset_index()
    Lp = D.pivot(index="pair_idx", columns="config", values="length_m")
    S["len_ratio_8_16_at_50m_median"] = (Lp.A2 / Lp.A3).median(); S["len_ratio_8_16_at_50m_max"] = (Lp.A2 / Lp.A3).max()
    S["len_ratio_8_16_at_250m_median"] = (Lp.A1 / Lp.A4).median(); S["len_ratio_8_16_at_250m_max"] = (Lp.A1 / Lp.A4).max()
    return D, S


def pair_traffic(st, records, beta=1.5):
    """Share of the canonical traffic mass carried by each pair (kernel weighting, as traffic_from_records)."""
    idx = {(tuple(st.attractors[a]), tuple(st.attractors[b])): k for k, (a, b) in enumerate(st.pair_candidates)}
    mass = np.zeros(len(st.pair_candidates))
    for r in records:
        c0 = r["costs"].min(); p = np.exp(-beta * (r["costs"] - c0)); p /= p.sum() + 1e-12
        (i1, j1), (i2, j2) = r["source"], r["target"]
        pw = np.exp(2.0 * st.E_grid[i1, j1]) * np.exp(2.0 * st.E_grid[i2, j2]) * np.exp(-3.0 * c0)
        k = idx.get((tuple(r["source"]), tuple(r["target"])), idx.get((tuple(r["target"]), tuple(r["source"]))))
        mass[k] += sum(pw * pp * len(path) for path, pp in zip(r["paths"], p))
    return 100 * mass / mass.sum()


if __name__ == "__main__":
    from final_evaluate import recget, check_block
    st = fc.state_for(fr.cfg())
    D, S = e2a(st)
    fc.EVAL.mkdir(parents=True, exist_ok=True)
    D.to_csv(fc.EVAL / "resolucion_E2A_pares.csv", index=False); S.to_csv(fc.EVAL / "resolucion_E2A_resumen.csv", index=False)
    check_block("S1")
    share = pair_traffic(st, recget("S1", "canon_s42"))
    A1 = D[D.config == "A1"].set_index("pair_idx")
    P = A1[["attr_a", "attr_b", "length_m", "straight_m", "sinuosity", "cost_map_units"]].copy()
    P["traffic_share_pct"] = share
    P.to_csv(fc.EVAL / "pares_rutas.csv")
    print(S.round(3).to_string(index=False))
    print("pairs with sinuosity > 2:", int((P.sinuosity > 2).sum()), "| their traffic share %:", round(P.loc[P.sinuosity > 2, "traffic_share_pct"].sum(), 2))
