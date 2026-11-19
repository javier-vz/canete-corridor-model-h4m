"""New Figure 7: Ward-4 regionalisation (Ward clustering with rook adjacency, 20 zones), mean simulated traffic F per zone
(canonical seed 42) and the number of site cells per zone. Stored partition and field only; zone boundaries are the cell
edges of the stored labels (no smoothing). Sites are not drawn: only counts per zone.
Panel a: zones filled by the rank of their mean F (1 = highest), zone identifiers. Panel b: site cells per zone, zones in the
same rank order, bars in the colour of the zone in panel a. Outputs (revision/out/figure7/): fig7_ward4_F_sitios.pdf/png,
fig7_datos.csv (zone, rank, mean F, area, site cells, density)."""
import os, sys
from pathlib import Path
W = Path(__file__).resolve().parents[2]; os.chdir(W); sys.path.insert(0, str(W / "code")); sys.dont_write_bytecode = True
import numpy as np, pandas as pd
from scipy import ndimage
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm, colors
import final_common as fc, final_runs as fr

OUT = W / "revision/out/figure7"; OUT.mkdir(parents=True, exist_ok=True)
st = fc.state_for(fr.cfg()); V, S = st.valid_mask, st.site_mask
F = np.load("out/final/runs/S1/canon_s42.npz")["F"]
Z = np.load("out/e3/partitions.npz"); g = np.full(V.shape, -1, int); g[Z["iy"], Z["ix"]] = Z["Ward_rook_n20"]
assert (g[V] >= 0).all() and (g[~V] == -1).all()
zones = np.unique(g[V])
d = pd.DataFrame({"zone": zones,
                  "mean_F": [F[g == z].mean() for z in zones],
                  "cells": [int((g == z).sum()) for z in zones],
                  "site_cells": [int(S[g == z].sum()) for z in zones]})
d["area_km2"] = d.cells * 0.0625; d["density_per_km2"] = d.site_cells / d.area_km2
d = d.sort_values(["mean_F", "zone"], ascending=[False, True]).reset_index(drop=True)
d["rank_F"] = d.mean_F.rank(ascending=False, method="min").astype(int)   # zones without traffic (F = 0) share the last rank
rho_c = spearmanr(d.mean_F, d.site_cells).statistic; rho_d = spearmanr(d.mean_F, d.density_per_km2).statistic
assert abs(rho_c - 0.6128) < 5e-4, rho_c                                   # = ablation, seed 42, Ward-4
assert d.site_cells.sum() == 84 and d.cells.sum() == int(V.sum())
d.to_csv(OUT / "fig7_datos.csv", index=False)

SURF, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#6b6a66", "#e1e0d9"
cmap = colors.LinearSegmentedColormap.from_list("blues", ["#0b3c7d", "#2a78d6", "#9ec5ee", "#e8f1fb"])   # one hue, dark = high
norm = colors.Normalize(vmin=1, vmax=20)
col = {r.zone: cmap(norm(r.rank_F)) for r in d.itertuples()}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "text.color": INK, "axes.edgecolor": MUTED,
                     "xtick.color": INK, "ytick.color": INK})
fig = plt.figure(figsize=(11.2, 4.6), facecolor=SURF)
ax = fig.add_axes([0.01, 0.19, 0.56, 0.74]); ax.set_facecolor(SURF)
rank = np.full(V.shape, np.nan)
for r in d.itertuples():
    rank[g == r.zone] = r.rank_F
ax.imshow(rank, origin="lower", cmap=cmap, norm=norm, interpolation="nearest")
segs = []
for i, j in zip(*np.nonzero((g[:, 1:] != g[:, :-1]) & ((g[:, 1:] >= 0) | (g[:, :-1] >= 0)))):
    segs.append([(j + 0.5, i - 0.5), (j + 0.5, i + 0.5)])
for i, j in zip(*np.nonzero((g[1:, :] != g[:-1, :]) & ((g[1:, :] >= 0) | (g[:-1, :] >= 0)))):
    segs.append([(j - 0.5, i + 0.5), (j + 0.5, i + 0.5)])
ax.add_collection(LineCollection(segs, colors="#ffffff", linewidths=0.9))
for r in d.itertuples():
    dist = ndimage.distance_transform_edt(np.pad(g == r.zone, 1))[1:-1, 1:-1]
    yy, xx = np.unravel_index(np.argmax(dist), dist.shape)
    ax.text(xx, yy, str(r.zone), ha="center", va="center", fontsize=7, color=INK,
            bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.9))
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.plot([3, 23], [2, 2], color=INK, lw=2); ax.text(13, 3.6, "5 km", ha="center", fontsize=7)
ax.set_title("a  Ward-4 zones by mean simulated traffic F (seed 42)", fontsize=9, loc="left")
cax = fig.add_axes([0.06, 0.125, 0.30, 0.022])
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_ticks([1, 5, 10, 15, 20]); cb.ax.tick_params(length=0, labelsize=7)
cb.set_label("zone rank by mean F (1 = highest; zones with F = 0 share the last rank)", fontsize=7.5)

bx = fig.add_axes([0.62, 0.22, 0.36, 0.66]); bx.set_facecolor(SURF)
x = np.arange(len(d))
bx.bar(x, d.site_cells, width=0.72, color=[col[z] for z in d.zone], edgecolor=SURF, linewidth=1.0, zorder=3)
for xi, r in zip(x, d.itertuples()):
    if r.site_cells > 0:
        bx.text(xi, r.site_cells + 0.8, str(r.site_cells), ha="center", va="bottom", fontsize=6.8, color=INK)
bx.set_xticks(x, [str(z) for z in d.zone], fontsize=6.8)
bx.set_xlabel("zone, ordered by mean F (highest at left)", fontsize=7.5)
bx.set_ylabel("site cells in the zone", fontsize=7.5)
bx.set_ylim(0, max(d.site_cells) * 1.12)
bx.grid(axis="y", color=GRID, lw=0.6, zorder=0); bx.tick_params(length=0)
for side in ("top", "right"):
    bx.spines[side].set_visible(False)
zero = [xi for xi, r in zip(x, d.itertuples()) if r.mean_F == 0]
if zero:
    hz = max(d.site_cells) * 0.10
    bx.plot([min(zero) - 0.4, max(zero) + 0.4], [hz, hz], color=MUTED, lw=0.8)
    bx.text((min(zero) + max(zero)) / 2, hz + 0.8, "no simulated traffic", ha="center", va="bottom", fontsize=6.6, color=MUTED)
bx.set_title("b  Site cells per zone", fontsize=9, loc="left")
bx.text(0.98, 0.95, f"Spearman ρ, mean F vs site cells = {rho_c:.3f}\nwith site density per km² = {rho_d:.3f}",
        transform=bx.transAxes, ha="right", va="top", fontsize=7, color=INK)
fig.text(0.01, 0.012, "Ward clustering with rook adjacency, 20 zones; 250 m cells, EPSG:32718; north is up. "
         "Zone identifiers match in both panels. Zones are units of comparison, not survey sectors or routes; site locations are not shown.", fontsize=6.8, color=MUTED)
fig.savefig(OUT / "fig7_ward4_F_sitios.pdf", dpi=600)   # raster of the zone map at 600 dpi; fig.savefig(OUT / "fig7_ward4_F_sitios.png", dpi=220)
print(d.round(4).to_string(index=False)); print("rho counts", round(rho_c, 4), "rho density", round(rho_d, 4))
