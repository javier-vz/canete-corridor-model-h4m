"""Geometric river masks on the analysis lattice (E1).

Uses exactly the HydroRIVERS selection rule of original_scripts/CA1.py:
all HydroRIVERS v1.0 lines clipped to the grid bbox + 5 km (no stream-order
filter). Verified: centre distances reproduce dist_river_m to 0.0000 m.
"""
import numpy as np
from paths import DATA
import geopandas as gpd
from shapely.geometry import box
from shapely.strtree import STRtree

RIVERS_GPKG = str(DATA / "HydroRIVERS_canete_UTM18S.gpkg")


def load_rivers(df, buffer_m=5000):
    riv = gpd.read_file(RIVERS_GPKG)
    xmin, ymin, xmax, ymax = df.X.min(), df.Y.min(), df.X.max(), df.Y.max()
    return gpd.clip(riv, box(xmin - buffer_m, ymin - buffer_m, xmax + buffer_m, ymax + buffer_m))


def river_cell_mask(state, rivers, cell_m=250.0):
    """True where the cell polygon intersects at least one river line."""
    lat = state.lattice
    h = cell_m / 2
    geoms = list(rivers.geometry)
    tree = STRtree(geoms)
    M = np.zeros(lat["shape"], dtype=bool)
    for x, y in zip(state.df.X, state.df.Y):
        cell = box(x - h, y - h, x + h, y + h)
        if len(tree.query(cell, predicate="intersects")):
            M[lat["y_to_iy"][y], lat["x_to_ix"][x]] = True
    return M & state.valid_mask
