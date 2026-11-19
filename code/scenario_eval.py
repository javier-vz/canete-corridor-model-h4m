"""Shared evaluator for E1/E6 scenarios. Re-uses ca_kernel_v2 for everything
algorithmic (no duplicated corridor code)."""
import sys, time, numpy as np, pandas as pd
from dataclasses import replace
from scipy.stats import spearmanr
import ca_kernel_v2 as ca
from paths import OUT as _OUT

OUT = str(_OUT)
BANDS = [(0, 125), (125, 250), (250, 500), (500, 1000), (1000, np.inf)]
BAND_LABELS = ["0-125", "125-250", "250-500", "500-1000", ">=1000"]


def with_attractors(state, exclusion_mask, E_grid=None):
    """New State whose attractors respect `exclusion_mask`; route graph untouched."""
    E = state.E_grid if E_grid is None else E_grid
    A, _, attr = ca.extract_local_maxima_adaptive(E, state.valid_mask, attractor_exclusion_mask=exclusion_mask,
                                                  peak_quantile=0.92, min_dist=5, n_bands=4, min_per_band=4)
    pairs = ca.build_local_pairs(attr, E, k_neighbors=6)
    return replace(state, A_mask=A, attractors=attr, pair_candidates=pairs, E_grid=E)


def zonal_rho(field, labels, iyix, sites):
    v = field[iyix[:, 0], iyix[:, 1]]
    out = []
    for lab in labels:
        g = pd.DataFrame({"z": lab, "v": v, "s": sites}).groupby("z")
        out.append(spearmanr(g.v.mean(), g.s.sum())[0])
    return np.array(out)


def evaluate(name, state, res, river_mask, dist, ref=None, extra=None):
    V = state.valid_mask
    F, L = res["F_norm"], res["L"]
    labels = np.load(f"{OUT}/kmeans_labels_100.npy"); iyix = np.load(f"{OUT}/rowmajor_iy_ix.npy")
    sites = np.load(f"{OUT}/rowmajor_sites.npy")
    rho_all = zonal_rho(L, labels, iyix, sites)
    ad = np.array([dist[i, j] for i, j in state.attractors])
    r = {"scenario": name, "n_attractors": len(state.attractors), "n_pairs": len(state.pair_candidates),
         "attractor_cells_on_river": int(sum(river_mask[i, j] for i, j in state.attractors)),
         "attr_dist_min": ad.min(), "attr_dist_median": np.median(ad), "attr_dist_max": ad.max(),
         "E_mean_attractors": float(np.mean([state.E_grid[i, j] for i, j in state.attractors])),
         "AUC_L": res["auc"], "rho20_L_canonical_P2": rho_all[42],
         "rho20_L_median_100km": np.median(rho_all), "rho20_L_p05": np.percentile(rho_all, 5),
         "rho20_L_p95": np.percentile(rho_all, 95),
         "F_on_river_cells": float(F[river_mask].sum() / F[V].sum()),
         "L_on_river_cells": float(L[river_mask].sum() / L[V].sum()),
         "lcc": res["lcc_size"], "m2": res["m2"]}
    for (lo, hi), lab in zip(BANDS, BAND_LABELS):
        m = V & (dist >= lo) & (dist < hi)
        r[f"F_{lab}"] = float(F[m].sum() / F[V].sum()); r[f"L_{lab}"] = float(L[m].sum() / L[V].sum())
        r[f"cells_{lab}"] = float(m.sum() / V.sum())
    if ref is not None:
        r["spearman_F_vs_H0"] = spearmanr(F[V], ref["F_norm"][V])[0]
        r["spearman_L_vs_H0"] = spearmanr(L[V], ref["L"][V])[0]
    if extra:
        r.update(extra)
    return r
