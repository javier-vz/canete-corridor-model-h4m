"""Read-only checks of the analysis grid (data/grilla_con_W_hidrica.csv) against its sources.

Added in release v2.0.0 to document how each column of the grid was obtained. It writes nothing.

1. Lattice: 12,826 cells whose centres lie on a regular 250 m lattice (EPSG:32718).
2. altitud is the value of the DEM pixel (data/DEM_vallebajo.tif, about 31 m) that contains the cell centre.
3. pendiente is the slope, in degrees, of that same pixel, computed by central differences (numpy.gradient) on the
   DEM with no-data pixels filled with the DEM mean (the method identified in the provenance check of the revision).
4. dist_river_m is the distance from the cell centre to the nearest HydroRIVERS v1.0 reach of
   data/HydroRIVERS_canete_UTM18S.gpkg clipped to the grid extent plus 5 km (code/river_geometry.load_rivers).
5. h_norm, p_norm and dist_river_norm are min-max normalisations, (x - min) / (max - min + 1e-12), of altitud,
   pendiente and dist_river_m; W_river = 1 - dist_river_norm; W_corridor = W_river * (1 - p_norm).
6. presencia_sitio_nuevo flags 84 cells. In those cells distancia_sitios_nuevos (distance from the cell centre to the
   nearest point of the archaeological reference set) is smaller than half the cell diagonal (176.8 m); in all the other
   cells it is larger than half the cell width (125 m).

Usage: python code/check_grid_columns.py   (exit status 0 if every check passes)."""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import numpy as np
import pandas as pd
import rasterio
import shapely

from river_geometry import load_rivers

GRID, DEM = "data/grilla_con_W_hidrica.csv", "data/DEM_vallebajo.tif"
CELL = 250.0
ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"[{'OK  ' if cond else 'FAIL'}] {name}{'  ' + detail if detail else ''}")


def minmax(s):
    s = s.astype(float)
    return (s - s.min()) / (s.max() - s.min() + 1e-12)


df = pd.read_csv(GRID)

# 1. lattice
on_lattice = ((df.X - df.X.min()) % CELL == 0).all() and ((df.Y - df.Y.min()) % CELL == 0).all()
check("1. 12,826 cells on a regular 250 m lattice", len(df) == 12826 and on_lattice and not df.duplicated(["X", "Y"]).any(),
      f"({len(df)} cells; X {df.X.min():.0f}-{df.X.max():.0f} m, Y {df.Y.min():.0f}-{df.Y.max():.0f} m)")

# 2-3. elevation and slope at the DEM pixel of the cell centre
with rasterio.open(DEM) as src:
    epsg, res = src.crs.to_epsg(), src.transform.a
    dem = src.read(1).astype(float)
    dem[dem == src.nodata] = np.nan
    r, c = map(np.array, rasterio.transform.rowcol(src.transform, df.X.values, df.Y.values))
check("   DEM in EPSG:32718", epsg == 32718, f"(pixel {res:.2f} m)")
z = dem[r, c]
check("2. altitud = DEM pixel at the cell centre", np.array_equal(z, df.altitud.values.astype(float)),
      f"({np.mean(z == df.altitud.values) * 100:.1f}% of cells)")
fill = np.where(np.isnan(dem), np.nanmean(dem), dem)
slope = np.degrees(np.arctan(np.hypot(*np.gradient(fill, res))))[r, c]
d_slope = np.abs(slope - df.pendiente.values).max()
check("3. pendiente = central-difference slope of that pixel", d_slope < 1e-4, f"(max |difference| {d_slope:.1e} degrees)")

# 4. distance to the nearest reach
riv = load_rivers(df)
union = shapely.union_all(riv.geometry.values)
dist = shapely.distance(shapely.points(df.X.values, df.Y.values), union)
d_dist = np.abs(dist - df.dist_river_m.values).max()
check("4. dist_river_m = distance to the nearest HydroRIVERS reach (grid extent + 5 km)", d_dist < 1e-6,
      f"({len(riv)} reaches; max |difference| {d_dist:.1e} m)")

# 5. derived columns
pairs = [("h_norm", minmax(df.altitud)), ("p_norm", minmax(df.pendiente)), ("dist_river_norm", minmax(df.dist_river_m)),
         ("W_river", 1.0 - df.dist_river_norm), ("W_corridor", df.W_river * (1.0 - df.p_norm))]
d_der = max(np.abs(df[col].values - v.values).max() for col, v in pairs)
check("5. h_norm, p_norm, dist_river_norm, W_river, W_corridor as documented", d_der < 1e-12, f"(max |difference| {d_der:.1e})")

# 6. site cells
s = df.presencia_sitio_nuevo == 1
half_diag = CELL / 2 * np.sqrt(2)
d_in, d_out = df.distancia_sitios_nuevos[s].max(), df.distancia_sitios_nuevos[~s].min()
check("6. 84 site cells, consistent with the distances to the reference sites",
      int(s.sum()) == 84 and set(df.presencia_sitio_nuevo.unique()) <= {0, 1} and d_in < half_diag and d_out > CELL / 2,
      f"({int(s.sum())} cells; max distance in site cells {d_in:.1f} m, min elsewhere {d_out:.1f} m)")

print("GRID CHECKS PASSED" if ok else "GRID CHECKS FAILED")
_sys.exit(0 if ok else 1)
