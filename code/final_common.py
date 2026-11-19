"""JCAA-80 final recompute (ADR-002 rev. 5; matrix validated 05/11): shared helpers, frozen.

Paths are relative to the project root (callers chdir there). Nothing here runs the model.
Procedures implemented here are those of out/final/PROCEDIMIENTOS.md.
"""
from __future__ import annotations

import hashlib
import json
import os
import pickle
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import ca_kernel_v2 as ca
import final_h4m_state as fs
from canonical_partition import canonical_partition

CSV = "data/grilla_con_W_hidrica.csv"
RIV = "data/celdas_fluviales_1020.csv"
E4_FIELDS = "out/e4/scenario_fields.pkl"          # Phase B slope variants (p, W, T, E, c), frozen
OUT = Path("out/final")
EVAL = OUT / "eval"
N_ATTR_S0, N_PAIRS_S0 = 25, 90
CODE_FILES = ["code/final_common.py", "code/final_runs.py", "code/final_h4m_state.py", "code/ca_kernel_v2.py",
              "code/canonical_partition.py"]
DATA_FILES = [CSV, RIV, E4_FIELDS]


# ----------------------------------------------------------------------------- files
def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def firma_calculo():
    """Signature of the code and data that produce the runs (every run JSON records it)."""
    parts = {p: sha256(p) for p in CODE_FILES + DATA_FILES}
    return hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest(), parts


def atomic_bytes(path, data: bytes):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_bytes(data)
    os.replace(tmp, path)


def atomic_json(path, obj):
    atomic_bytes(path, json.dumps(obj, indent=1, sort_keys=True).encode())


def atomic_npz(path, **arrays):
    import io
    buf = io.BytesIO()
    np.savez_compressed(buf, **arrays)
    atomic_bytes(path, buf.getvalue())


# ----------------------------------------------------------------------------- states
_CACHE = {}


def state_for(cfg: dict):
    """Model state of one configuration (cached). cfg keys: rules ('H4m'|'v2'), a, b, q, T (r, s, w),
    slope (None or a Phase B variant name)."""
    key = json.dumps(cfg, sort_keys=True)
    if key in _CACHE:
        return _CACHE[key]
    a, b, q = tuple(cfg["a"]), tuple(cfg["b"]), cfg["q"]
    r, s, w = cfg["T"]
    st = fs.build_state_h4m(CSV, RIV, a_weights=a, b_weights=b, peak_quantile=q, T_r=r, T_s=s, T_window=int(w),
                            rules=cfg["rules"])
    if cfg.get("slope"):
        st = slope_variant_state(st, cfg["slope"])
    _CACHE[key] = st
    return st


def masks_and_dist(st):
    R, Rm = fs.river_masks(st, RIV)
    d = fs.dist_grid(st)
    return R, Rm, np.where(np.isnan(d), 1e9, d)


def slope_variant_state(st_h4m, variant):
    """S9: H4m rules applied to a Phase B slope variant (p, W, T, E, c computed there with the nominal
    weights). Same steps and constants as final_h4m_state.build_state_h4m."""
    F = pickle.load(open(E4_FIELDS, "rb"))[variant]
    R, _, d0 = masks_and_dist(st_h4m)
    V = st_h4m.valid_mask
    E1 = F["E"] * (1.0 - fs.DELTA_H2M * np.exp(-d0 / fs.LAMBDA_M))
    A, _, attr = ca.extract_local_maxima_adaptive(E1, V, attractor_exclusion_mask=R | (d0 <= fs.EXCL_D_M),
                                                  peak_quantile=fs.Q_NOMINAL, min_dist=fs.NMS_MIN_DIST,
                                                  n_bands=fs.N_BANDS, min_per_band=fs.MIN_PER_BAND)
    pairs = ca.build_local_pairs(attr, E1, k_neighbors=fs.K_NEIGHBORS)
    c1 = F["c"] + fs.GAMMA_H3 * R
    return replace(st_h4m, p_grid=F["p"], w_grid=F["W"], T_grid=F["T"], E_grid=E1, c_grid=c1, A_mask=A,
                   attractors=attr, pair_candidates=pairs, G_template=ca.build_graph(c1, V))


def state_key_hash(st):
    fp = fs.fingerprint(st)
    return hashlib.sha256(json.dumps(fp, sort_keys=True).encode()).hexdigest()


# ----------------------------------------------------------------------------- partitions and zonal stats
_PART = {}


def partition(st, n=20):
    """P2 procedure (row-major cells, K-means random_state 42, n_init 10) for n zones; cached."""
    if n not in _PART:
        _PART[n] = canonical_partition(st, n_zones=n)
    return _PART[n]


def zone_frame(field, part):
    v = field[part.iy.values, part.ix.values]
    g = pd.DataFrame({"zone": part.zone.values, "v": v, "s": part.presencia_sitio_nuevo.values}).groupby("zone")
    return pd.DataFrame({"L_mean": g.v.mean(), "L_max": g.v.max(), "n_sites": g.s.sum(), "n_cells": g.v.size()})


def zonal_stats(field, part):
    """rho (mean and max), p, leave-one-out without the zone with most sites (PROCEDIMIENTOS 3)."""
    z = zone_frame(field, part)
    r, p = spearmanr(z.L_mean, z.n_sites); rm, pm = spearmanr(z.L_max, z.n_sites)
    dom = z.n_sites.idxmax(); loo = z.drop(dom)
    rl, pl = spearmanr(loo.L_mean, loo.n_sites)
    return dict(n_zones=len(z), rho=r, p=p, rho_max=rm, p_max=pm, rho_loo=rl, p_loo=pl,
                zones_with_sites=int((z.n_sites > 0).sum()), dominant_zone=int(dom), dominant_sites=int(z.n_sites.max()))


def zone_prefix(field, part):
    """Complete zones ranked by mean field, cumulative area/sites/lift (no interpolation)."""
    z = zone_frame(field, part).sort_values("L_mean", ascending=False)
    N, S = z.n_cells.sum(), z.n_sites.sum()
    z["rank"] = np.arange(1, len(z) + 1)
    z["cum_area_pct"] = 100 * z.n_cells.cumsum() / N
    z["cum_sites"] = z.n_sites.cumsum()
    z["cum_sites_pct"] = 100 * z.cum_sites / S
    z["lift"] = (z.cum_sites_pct / 100) / (z.cum_area_pct / 100)
    return z.reset_index()


def rho20(field, st):
    return float(zonal_stats(field, partition(st, 20))["rho"])
