"""S3: Ward with rook (4-neighbour) adjacency as a candidate contiguous regionalisation (no partition is recomputed).

Reads the stored partitions (out/e3/partitions.npz), the H4m state (valid cells, site cells) and the canonical L of seed 42.
For Ward rook, Ward queen, SKATER rook and the reference K-means (20 zones; 30 and 40 for diagnostics): 4-neighbour
connectivity of each zone, holes, area, perimeter, Polsby-Popper compactness, elongation and size dispersion. The choice
criteria are geographic and fixed before looking at sites: every zone connected under 4 neighbours, no zone below 1 % of
the domain, size dispersion and compactness. The zonal rho of Table 12 is reported only as a secondary column.
Outputs (04_S3_WARD4/): S3_WARD4_ZONAS.csv (20 zones), S3_COMPARACION_PARTICIONES.csv, S3_MAPA.pdf/png, S3_log.txt."""
import os, sys
from pathlib import Path
W = Path(__file__).resolve().parents[2]; os.chdir(W); sys.path.insert(0, str(W / "code")); sys.dont_write_bytecode = True
import numpy as np, pandas as pd
from scipy import ndimage
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import final_common as fc, final_runs as fr

OUT = W / "revision/out/S5_ward4"; OUT.mkdir(parents=True, exist_ok=True)
st = fc.state_for(fr.cfg()); V, S = st.valid_mask, st.site_mask
L = np.load("out/final/runs/S1/canon_s42.npz")["L"]
Z = np.load("out/e3/partitions.npz")
N, CELL_KM2, EDGE_KM = int(V.sum()), 0.0625, 0.25
ROOK = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])
t12 = pd.read_csv("out/final/eval/tabla12_particiones.csv")


def grid_of(key):
    g = np.full(V.shape, -1, int); g[Z["iy"], Z["ix"]] = Z[key]
    assert (g[V] >= 0).all() and (g[~V] == -1).all(), "labels must cover exactly the valid cells"
    return g


def zone_stats(g, z):
    m = g == z
    lab, ncomp = ndimage.label(m, structure=ROOK)
    comp_sizes = np.bincount(lab.ravel())[1:]
    cl, nc = ndimage.label(~m, structure=ROOK)                       # holes: complement components not touching the border
    border = set(np.unique(np.r_[cl[0], cl[-1], cl[:, 0], cl[:, -1]])) - {0}
    holes = len(set(range(1, nc + 1)) - border)
    pad = np.pad(m, 1)
    edges = sum(int((pad & ~np.roll(pad, s, axis=a)).sum()) for a in (0, 1) for s in (1, -1))
    area, per = m.sum() * CELL_KM2, edges * EDGE_KM
    iy, ix = np.nonzero(m)
    ev = np.sort(np.linalg.eigvalsh(np.cov(np.c_[ix, iy].T))) if m.sum() > 2 else np.array([1, 1])
    return dict(cells=int(m.sum()), area_km2=area, pct_domain=100 * m.sum() / N, components_4n=int(ncomp),
                largest_component_pct=100 * comp_sizes.max() / m.sum(), holes=holes, perimeter_km=per,
                polsby_popper=4 * np.pi * area / per ** 2, elongation=float(np.sqrt(ev[1] / max(ev[0], 1e-9))),
                sites=int(S[m].sum()), mean_L=float(L[m].mean()), cy=float(iy.mean()), cx=float(ix.mean()))


rows, comp = [], []
for key in [k for k in Z.files if k not in ("iy", "ix")]:
    alg, n = key.rsplit("_n", 1); g = grid_of(key)
    zs = pd.DataFrame([dict(zone=int(z), **zone_stats(g, z)) for z in np.unique(g[V])])
    r12 = t12[(t12.algorithm == alg) & (t12.n == int(n))]
    comp.append(dict(partition=alg, n=int(n), zones=len(zs), zones_not_connected_4n=int((zs.components_4n > 1).sum()),
                     zones_with_holes=int((zs.holes > 0).sum()), zones_below_1pct=int((zs.pct_domain < 1).sum()),
                     area_min_km2=zs.area_km2.min(), area_median_km2=zs.area_km2.median(), area_max_km2=zs.area_km2.max(),
                     area_cv=zs.area_km2.std(ddof=0) / zs.area_km2.mean(), polsby_popper_median=zs.polsby_popper.median(),
                     elongation_median=zs.elongation.median(), rho_secondary=float(r12.rho.iloc[0]) if len(r12) else np.nan))
    if key == "Ward_rook_n20":
        w4 = zs
C = pd.DataFrame(comp)
C.to_csv(OUT / "S3_COMPARACION_PARTICIONES.csv", index=False)

# table of the 20 Ward-4 zones: ranked by mean L, cumulative complete zones (no interpolation)
w4 = w4.sort_values("mean_L", ascending=False).reset_index(drop=True)
w4["rank_by_mean_L"] = np.arange(1, len(w4) + 1)
w4["site_density_per_km2"] = w4.sites / w4.area_km2
w4["cum_area_pct"] = w4.pct_domain.cumsum(); w4["cum_sites_pct"] = 100 * w4.sites.cumsum() / int(S.sum())
cols = ["zone", "rank_by_mean_L", "mean_L", "cells", "area_km2", "pct_domain", "sites", "site_density_per_km2", "cum_area_pct", "cum_sites_pct",
        "components_4n", "holes", "perimeter_km", "polsby_popper", "elongation"]
w4[cols].to_csv(OUT / "S3_WARD4_ZONAS.csv", index=False)
assert w4.cells.sum() == N and w4.sites.sum() == int(S.sum())

# map: zones filled by rank of mean L (one sequential hue), boundaries, zone numbers; no sites
g = grid_of("Ward_rook_n20")
rank = np.full(V.shape, np.nan); rmap = dict(zip(w4.zone, w4.rank_by_mean_L))
for z, rk in rmap.items():
    rank[g == z] = rk
SURF, TXT, MUTED = "#fcfcfb", "#0b0b0b", "#898781"
fig, ax = plt.subplots(figsize=(9.5, 6.4), facecolor=SURF); ax.set_facecolor(SURF)
im = ax.imshow(rank, origin="lower", cmap="Blues_r", vmin=1, vmax=24, interpolation="nearest")
from matplotlib.collections import LineCollection
segs = []
for i, j in zip(*np.nonzero((g[:, 1:] != g[:, :-1]) & (g[:, 1:] >= 0) & (g[:, :-1] >= 0))):
    segs.append([(j + 0.5, i - 0.5), (j + 0.5, i + 0.5)])           # vertical edge between (i, j) and (i, j + 1)
for i, j in zip(*np.nonzero((g[1:, :] != g[:-1, :]) & (g[1:, :] >= 0) & (g[:-1, :] >= 0))):
    segs.append([(j - 0.5, i + 0.5), (j + 0.5, i + 0.5)])           # horizontal edge between (i, j) and (i + 1, j)
ax.add_collection(LineCollection(segs, colors=TXT, linewidths=0.6))
for r in w4.itertuples():                                            # label at the cell farthest from the zone boundary
    dist = ndimage.distance_transform_edt(np.pad(g == r.zone, 1))[1:-1, 1:-1]
    yy, xx = np.unravel_index(np.argmax(dist), dist.shape)
    ax.text(xx, yy, str(r.zone), ha="center", va="center", fontsize=7, color=TXT,
            bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.85))
ax.set_xticks([]); ax.set_yticks([])
ax.plot([4, 24], [3, 3], color=TXT, lw=2); ax.text(14, 4.5, "5 km", ha="center", fontsize=7)
cb = fig.colorbar(im, ax=ax, shrink=0.6); cb.set_label("rank of the zone by mean L (1 = highest)", fontsize=8); cb.ax.invert_yaxis()
ax.set_title("Candidate contiguous landscape strata: Ward clustering with rook adjacency, 20 zones (not survey units)", fontsize=9, loc="left")
fig.text(0.01, 0.01, "Numbers: zone identifiers. 250 m cells, EPSG:32718; north is up. Archaeological sites are not shown.", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(OUT / "S3_MAPA.pdf"); fig.savefig(OUT / "S3_MAPA.png", dpi=200)
log = [C[C.n == 20].round(3).to_string(index=False), "", w4[cols].round(3).to_string(index=False)]
(OUT / "S3_log.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
print(C.round(3).to_string(index=False))
print(w4[["zone", "rank_by_mean_L", "area_km2", "pct_domain", "sites", "components_4n", "holes", "polsby_popper", "elongation"]].round(2).to_string(index=False))
