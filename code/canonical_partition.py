"""Single source of truth for the paper's 20-zone K-means partition.

Paper (Sec. 3.5): features (column index, row index, h_norm, p_norm),
standardised; valid cells ordered by row and then column; K-means with
random_state=42 and n_init=10.  Row order matters: K-means++ seeding
with a fixed random_state draws different initial centres when the
rows are permuted, so every script must call this function instead of
re-implementing the partition.
"""
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


def canonical_partition(state, n_zones=20, random_state=42, n_init=10):
    """Return a DataFrame (row-major order) with iy, ix, X, Y, site flag and zone."""
    lat = state.lattice
    d = state.df[["X", "Y", "h_norm", "p_norm", "presencia_sitio_nuevo"]].copy()
    d["iy"] = [lat["y_to_iy"][y] for y in d.Y]
    d["ix"] = [lat["x_to_ix"][x] for x in d.X]
    d = d.sort_values(["iy", "ix"]).reset_index(drop=True)
    X = StandardScaler().fit_transform(d[["ix", "iy", "h_norm", "p_norm"]].to_numpy())
    d["zone"] = KMeans(n_clusters=n_zones, random_state=random_state,
                       n_init=n_init).fit(X).labels_
    return d


def zone_table(part, field):
    """Zone means of a 2-D field, site counts, cell counts."""
    v = field[part.iy.values, part.ix.values]
    g = pd.DataFrame({"zone": part.zone, "v": v, "s": part.presencia_sitio_nuevo}).groupby("zone")
    return pd.DataFrame({"mean": g.v.mean(), "max": g.v.max(), "n_sites": g.s.sum(), "n_cells": g.v.size()})
