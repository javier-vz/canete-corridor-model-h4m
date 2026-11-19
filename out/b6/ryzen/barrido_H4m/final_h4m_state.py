"""JCAA-80, final recompute (ADR-002 rev. 5): frozen construction of the H4m model state.

Uses the UNMODIFIED kernel (ca_kernel_v2.py, sha256 15d1409a...) and adds the approved
hydrological rules on top of ca.build_state, exactly as h_scenarios.make_state("H4m") did in
Phases B-B5:

  H1  attractor exclusion: no attractor on a river cell R (cell polygon intersects one of the
      HydroRIVERS reaches of the model window; 1,020 cells) nor at dist_river_m <= 125 m
      (the 829 cells at <= 125 m are a subset of R).
  H2m E' = E * (1 - 0.25 * exp(-dist_river_m / 250 m)); E' is used for attractor extraction and
      pair priorities (the kernel reads state.E_grid).
  H3  c' = c + 0.10 * R on all 1,020 river cells (cost_mask="all"); the routing graph is rebuilt
      from c' with the kernel's build_graph.

Sensitivity families pass their perturbations to ca.build_state (weights a, b; T parameters;
attractor quantile q) and the H4m rules are applied on top, so every configuration is "H4m with
that perturbation". rules="v2" rebuilds the published v2 state (125 m buffer, no H2m/H3); it is
only used to check that the perturbation designs reproduce the attractor and pair counts of the
published tables.

Inputs: the canonical grid CSV (dist_river_m is read from it) and a CSV with the coordinates of
the river cells (X, Y, main), written by write_river_cells_csv() from out/river_cell_mask.npy.
"""
from __future__ import annotations

import hashlib
from dataclasses import replace

import numpy as np
import pandas as pd

import ca_kernel_v2 as ca

DELTA_H2M = 0.25        # H2m amplitude
LAMBDA_M = 250.0        # H2m decay length (m)
GAMMA_H3 = 0.10         # H3 cost penalty
EXCL_D_M = 125.0        # H1 distance rule (m)
NMS_MIN_DIST, N_BANDS, MIN_PER_BAND, K_NEIGHBORS = 5, 4, 4, 6   # kernel defaults (paper Sec. 3.2)

A_NOMINAL = (0.10, 0.20, 0.50, 0.20)
B_NOMINAL = (0.18, 0.36, 0.26, 0.20)
T_NOMINAL = (0.65, 0.20, 5)
Q_NOMINAL = 0.92


# ----------------------------------------------------------------------------- designs
def perturb(vec, idx, sign, delta=0.20):
    """Paper Eq. (redist); identical to rerun_weights_ablation_v2.perturb."""
    vec = list(vec)
    orig_x = vec[idx]
    new_x = (1 + sign * delta) * orig_x
    scale = (1.0 - new_x) / (1.0 - orig_x) if orig_x < 1 else 1.0
    return tuple(vec[j] * scale if j != idx else new_x for j in range(len(vec)))


def weight_configs():
    """17 configurations, same labels and order as rerun_weights_ablation_v2.build_configs."""
    names = ["h", "p", "w", "t"]
    cfgs = [("baseline", A_NOMINAL, B_NOMINAL)]
    for idx, name in enumerate(names):
        for sign, s_lbl in [(+1, "+20"), (-1, "-20")]:
            cfgs.append((f"a_{name}_{s_lbl}", perturb(A_NOMINAL, idx, sign), B_NOMINAL))
    for idx, name in enumerate(names):
        for sign, s_lbl in [(+1, "+20"), (-1, "-20")]:
            cfgs.append((f"b_{name}_{s_lbl}", A_NOMINAL, perturb(B_NOMINAL, idx, sign)))
    return cfgs


def T_configs():
    """7 configurations, same labels and order as rerun_T_ablation_v2.build_configs."""
    r0, s0, w0 = T_NOMINAL
    cfgs = [("baseline", r0, s0, w0)]
    for r in (0.50, 0.80):
        cfgs.append((f"r={r:.2f}", r, s0, w0))
    for s in (0.00, 0.40):
        cfgs.append((f"s={s:.2f}", r0, s, w0))
    for w in (3, 7):
        cfgs.append((f"w={w}", r0, s0, w))
    return cfgs


QUANTILES = (0.90, 0.92, 0.95, 0.975)


# ----------------------------------------------------------------------------- inputs
def write_river_cells_csv(state, river_mask, main_mask, path):
    """Coordinates of the river cells (robust to array layout across machines)."""
    lat = state.lattice
    ix_to_x = {v: k for k, v in lat["x_to_ix"].items()}
    iy_to_y = {v: k for k, v in lat["y_to_iy"].items()}
    iy, ix = np.nonzero(river_mask & state.valid_mask)
    pd.DataFrame({"X": [ix_to_x[i] for i in ix], "Y": [iy_to_y[i] for i in iy],
                  "main": main_mask[iy, ix].astype(int)}).to_csv(path, index=False)


def river_masks(state, river_csv):
    rc = pd.read_csv(river_csv)
    lat = state.lattice
    R = np.zeros(state.valid_mask.shape, bool)
    Rm = np.zeros(state.valid_mask.shape, bool)
    iy = np.array([lat["y_to_iy"][y] for y in rc.Y]); ix = np.array([lat["x_to_ix"][x] for x in rc.X])
    R[iy, ix] = True
    Rm[iy[rc.main.values == 1], ix[rc.main.values == 1]] = True
    return R, Rm


def dist_grid(state):
    lat = state.lattice
    d = np.full(state.valid_mask.shape, np.nan)
    for x, y, v in zip(state.df.X, state.df.Y, state.df.dist_river_m):
        d[lat["y_to_iy"][y], lat["x_to_ix"][x]] = v
    return d


# ----------------------------------------------------------------------------- state
def build_state_h4m(csv_path, river_csv, a_weights=A_NOMINAL, b_weights=B_NOMINAL,
                    peak_quantile=Q_NOMINAL, T_r=T_NOMINAL[0], T_s=T_NOMINAL[1], T_window=T_NOMINAL[2],
                    rules="H4m", cost_mask="all"):
    """Return the kernel State for one configuration.

    rules="H4m": H1 + H2m + H3 as described in the module docstring (final configuration).
    rules="v2":  the published v2 state (ca.build_state with river_buffer_m=125); checks only.
    cost_mask="main" gives H4m-main (penalty on the 218 main-channel cells); sensitivity only.
    """
    if rules == "v2":
        return ca.build_state(csv_path, a_weights=a_weights, b_weights=b_weights,
                              peak_quantile=peak_quantile, T_r=T_r, T_s=T_s, T_window=T_window,
                              river_buffer_m=125)
    if rules != "H4m":
        raise ValueError(rules)
    st0 = ca.build_state(csv_path, a_weights=a_weights, b_weights=b_weights,
                         peak_quantile=peak_quantile, T_r=T_r, T_s=T_s, T_window=T_window,
                         river_buffer_m=0)
    R, Rmain = river_masks(st0, river_csv)
    d = dist_grid(st0)
    d0 = np.where(np.isnan(d), 1e9, d)
    E1 = st0.E_grid * (1.0 - DELTA_H2M * np.exp(-d0 / LAMBDA_M))
    excl = R | (d0 <= EXCL_D_M)
    A, _, attr = ca.extract_local_maxima_adaptive(E1, st0.valid_mask, attractor_exclusion_mask=excl,
                                                  peak_quantile=peak_quantile, min_dist=NMS_MIN_DIST,
                                                  n_bands=N_BANDS, min_per_band=MIN_PER_BAND)
    pairs = ca.build_local_pairs(attr, E1, k_neighbors=K_NEIGHBORS)
    Rc = R if cost_mask == "all" else Rmain
    c1 = st0.c_grid + GAMMA_H3 * Rc
    st = replace(st0, E_grid=E1, A_mask=A, attractors=attr, pair_candidates=pairs,
                 c_grid=c1, G_template=ca.build_graph(c1, st0.valid_mask))
    assert not any(R[i, j] for i, j in st.attractors), "H1 violated"
    return st


def fingerprint(state):
    """Compact identity of a state, for cross-machine checks."""
    V = state.valid_mask
    h = lambda a: hashlib.sha256(np.round(np.asarray(a, float)[V], 12).tobytes()).hexdigest()
    return {"n_attractors": len(state.attractors), "n_pairs": len(state.pair_candidates),
            "attractors": [list(map(int, a)) for a in state.attractors],
            "pairs": [list(map(int, p)) for p in state.pair_candidates],
            "sha256_E": h(state.E_grid), "sha256_c": h(state.c_grid), "sha256_T": h(state.T_grid)}
