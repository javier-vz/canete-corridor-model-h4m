"""JCAA-80 final recompute: files for the QGIS maps (Figures 4-7), EPSG:32718, 250 m grid.
Orientation: the kernel lattice has row 0 = south; GeoTIFF row 0 = north, so arrays are flipped
once and every raster is read back at a known cell (same check as fig6_files_P2.py).
Fig 4: h_norm, p_norm, W, T, E' (E_prime), c' (c_prime)      Fig 5: E', attractors, sites
Fig 6: F deterministic and canonical (seed 42), L canonical, routes of both runs
Fig 7: mean L per P2 zone, zone ids, centroids with site counts
Output: out/final/qgis/"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import LineString, Point

import final_common as fc
import final_runs as fr

Q = fc.OUT / "qgis"
TR = from_origin(342700.0, 8570750.0, 250, 250)


def write_tif(arr, name, st, dtype="float32", nodata=np.nan):
    a = np.where(st.valid_mask, arr, nodata).astype(dtype)
    with rasterio.open(Q / name, "w", driver="GTiff", height=a.shape[0], width=a.shape[1], count=1, dtype=dtype,
                       crs="EPSG:32718", transform=TR, nodata=nodata) as d:
        d.write(np.flipud(a), 1)
    with rasterio.open(Q / name) as d:                       # read-back orientation check
        iy, ix = map(int, np.argwhere(st.valid_mask)[len(np.argwhere(st.valid_mask)) // 3])
        x = st.lattice["xs"][ix]; y = st.lattice["ys"][iy]
        r, c = d.index(x, y)
        assert np.isclose(d.read(1)[r, c], a[iy, ix]), f"orientation error in {name}"


def routes_gdf(st, records, beta=1.5):
    xs, ys = np.array(st.lattice["xs"]), np.array(st.lattice["ys"])
    rows = []
    for k, r in enumerate(records):
        c = np.asarray(r["costs"]); p = np.exp(-beta * (c - c.min())); p /= p.sum()
        for m, (path, pp, cc) in enumerate(zip(r["paths"], p, c)):
            rows.append(dict(pair=k, alternative=m + 1, probability=float(pp), cost_perturbed=float(cc), n_cells=len(path),
                             geometry=LineString(np.c_[xs[path[:, 1]], ys[path[:, 0]]])))
    return gpd.GeoDataFrame(rows, crs="EPSG:32718")


if __name__ == "__main__":
    from final_evaluate import fget, recget, check_block
    Q.mkdir(parents=True, exist_ok=True)
    check_block("S1")
    st = fc.state_for(fr.cfg())
    Fc, Lc = fget("S1", "canon_s42"); Fd, Ld = fget("S1", "det_s42")
    for name, arr in [("h_norm.tif", st.h_grid), ("p_norm.tif", st.p_grid), ("W.tif", st.w_grid), ("T.tif", st.T_grid),
                      ("E_prime.tif", st.E_grid), ("c_prime.tif", st.c_grid), ("F_det_s42.tif", Fd), ("F_canon_s42.tif", Fc),
                      ("L_canon_s42.tif", Lc)]:
        write_tif(arr, name, st)
    xs, ys = np.array(st.lattice["xs"]), np.array(st.lattice["ys"])
    R, _, d0 = fc.masks_and_dist(st)
    A = pd.DataFrame([dict(attractor=k, X=xs[j], Y=ys[i], E_prime=st.E_grid[i, j], dist_river_m=d0[i, j])
                      for k, (i, j) in enumerate(st.attractors)])
    A.to_csv(Q / "atractores_H4m.csv", index=False)
    gpd.GeoDataFrame(A, geometry=[Point(x, y) for x, y in zip(A.X, A.Y)], crs="EPSG:32718").to_file(Q / "atractores_H4m.gpkg")
    iy, ix = np.nonzero(st.site_mask)
    Sx = pd.DataFrame(dict(X=xs[ix], Y=ys[iy]))
    gpd.GeoDataFrame(Sx, geometry=[Point(x, y) for x, y in zip(Sx.X, Sx.Y)], crs="EPSG:32718").to_file(Q / "sitios_84.gpkg")
    routes_gdf(st, recget("S1", "canon_s42")).to_file(Q / "rutas_s42.gpkg", layer="rutas_canonica")
    routes_gdf(st, recget("S1", "det_s42")).to_file(Q / "rutas_s42.gpkg", layer="rutas_determinista")
    part = fc.partition(st, 20)
    z = fc.zone_frame(Lc, part).sort_values("L_mean", ascending=False); z["rank_by_L"] = np.arange(1, len(z) + 1)
    zm = np.full(st.valid_mask.shape, np.nan); zid = np.full(st.valid_mask.shape, -1)
    zm[part.iy, part.ix] = z.loc[part.zone, "L_mean"].values; zid[part.iy, part.ix] = part.zone.values
    write_tif(zm, "zonas_Lmedio_P2.tif", st); write_tif(zid, "zonas_id_P2.tif", st, dtype="int16", nodata=-1)
    cent = part.groupby("zone").agg(X=("X", "mean"), Y=("Y", "mean")).join(z).reset_index().sort_values("rank_by_L")
    cent.to_csv(Q / "zonas_centroides_P2.csv", index=False)
    gpd.GeoDataFrame(cent, geometry=[Point(x, y) for x, y in zip(cent.X, cent.Y)], crs="EPSG:32718").to_file(Q / "zonas_centroides_P2.gpkg")
    print("QGIS files:", sorted(p.name for p in Q.iterdir()))
