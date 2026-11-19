"""Zonal diagnostic: site counts versus site densities, for the four stored 20-zone partitions.

For each partition (K-means P2 and the three contiguity-constrained ones) and each field (1 - h, T, E', F, L; F and L
for seed 42 and for each of the 30 canonical realisations 2000-2029), the zone mean of the field is compared with
(a) the number of site cells in the zone and (b) the site density (site cells per km2). The relation between zone area
and site count is reported too. Spearman rank correlations only, n = 20 zones; no p-values are reported because zones are
spatially dependent and the 84 sites are the same in every comparison. Stored fields and partitions only; nothing is simulated.
Outputs (revision/out/S5_counts_densities/): ZONAS_AREA_SITIOS.csv, RHO_CONTEOS_DENSIDADES.csv, diagnostico_log.txt."""
import glob, os, sys
from pathlib import Path
W = Path(__file__).resolve().parents[2]; os.chdir(W); sys.path.insert(0, str(W / "code")); sys.dont_write_bytecode = True
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import final_common as fc, final_runs as fr
import ca_kernel_v2 as ca

OUT = W / "revision/out/S5_counts_densities"; OUT.mkdir(parents=True, exist_ok=True)
st = fc.state_for(fr.cfg()); V, S = st.valid_mask, st.site_mask
Z = np.load("out/e3/partitions.npz")
PART = {"K-means P2": "KMeans_P2_rs42_n20", "Ward-4 (rook)": "Ward_rook_n20", "Ward-8 (queen)": "Ward_queen_n20", "SKATER (rook)": "SKATER_rook_n20"}
fixed = {"1 - h": 1 - st.h_grid, "T": st.T_grid, "E'": st.E_grid}
assert abs(ca.auc_score(st.E_grid, S, V) - 0.808870) < 1e-5 and abs(ca.auc_score(st.T_grid, S, V) - 0.817657) < 1e-5   # = Table 10
runs = {42: "out/final/runs/S1/canon_s42.npz", **{int(p.split("_s")[-1][:-4]): p for p in sorted(glob.glob("out/final/runs/S2/canon_s*.npz"))}}
assert len(runs) == 31
CELL_KM2 = 0.0625
zrows, rrows, lrows = [], [], []
for pname, key in PART.items():
    lab = Z[key]; iy, ix = Z["iy"], Z["ix"]
    sites = S[iy, ix].astype(int)
    zt = pd.DataFrame({"zone": lab, "site": sites}).groupby("zone").agg(cells=("site", "size"), sites=("site", "sum"))
    zt["area_km2"] = zt.cells * CELL_KM2; zt["density_per_km2"] = zt.sites / zt.area_km2
    rho_area_count = spearmanr(zt.area_km2, zt.sites).statistic
    rho_area_density = spearmanr(zt.area_km2, zt.density_per_km2).statistic
    def zone_mean(field):
        return pd.Series(field[iy, ix]).groupby(lab).mean().reindex(zt.index)
    for fname, field in fixed.items():
        m = zone_mean(field)
        rrows.append(dict(partition=pname, field=fname, realisation="deterministic", rho_count=spearmanr(m, zt.sites).statistic,
                          rho_density=spearmanr(m, zt.density_per_km2).statistic, rho_field_area=spearmanr(m, zt.area_km2).statistic))
        zt[f"mean_{fname}"] = m.values
        dom = zt.sites.idxmax()
        lrows.append(dict(partition=pname, field=fname, realisation="deterministic", zone_left_out=int(dom), sites_left_out=int(zt.sites.max()),
                          rho_count_19=spearmanr(m.drop(dom), zt.sites.drop(dom)).statistic))
    for seed, path in runs.items():
        npz = np.load(path)
        for fname in ("F", "L"):
            m = zone_mean(npz[fname])
            rrows.append(dict(partition=pname, field=fname, realisation=f"seed {seed}", rho_count=spearmanr(m, zt.sites).statistic,
                              rho_density=spearmanr(m, zt.density_per_km2).statistic, rho_field_area=spearmanr(m, zt.area_km2).statistic))
            if seed == 42:
                zt[f"mean_{fname}_s42"] = m.values
            dom = zt.sites.idxmax()
            lrows.append(dict(partition=pname, field=fname, realisation=f"seed {seed}", zone_left_out=int(dom), sites_left_out=int(zt.sites.max()),
                              rho_count_19=spearmanr(m.drop(dom), zt.sites.drop(dom)).statistic))
    zt = zt.reset_index(); zt.insert(0, "partition", pname)
    zt["rho_area_count"] = rho_area_count; zt["rho_area_density"] = rho_area_density
    zrows.append(zt)
ZT = pd.concat(zrows, ignore_index=True); RT = pd.DataFrame(rrows)
ZT.to_csv(OUT / "ZONAS_AREA_SITIOS.csv", index=False); RT.to_csv(OUT / "RHO_CONTEOS_DENSIDADES.csv", index=False)
# without the zone with most site cells (19 zones): seed 42 and range over the 31 realisations for F and L
LT = pd.DataFrame(lrows)
LS = LT.assign(field=LT.field.map(lambda f: f)).groupby(["partition", "field"]).agg(zone_left_out=("zone_left_out", "first"), sites_left_out=("sites_left_out", "first"),
                                                    rho_count_19_min=("rho_count_19", "min"), rho_count_19_max=("rho_count_19", "max")).reset_index()
LS.to_csv(OUT / "RHO_SIN_ZONA_DOMINANTE.csv", index=False)

# controls against frozen results: rho with counts for L seed 42 and fixed fields under P2 = Table 10; F, L Ward-4 = ablation
t10 = pd.read_csv("out/final/eval/tabla10_predictores.csv").set_index("predictor").rho20
chk = {("K-means P2", "1 - h"): t10["1 - h"], ("K-means P2", "T"): t10["T"], ("K-means P2", "E'"): t10["E' (final affinity)"]}
for (p, f), v in chk.items():
    got = RT[(RT.partition == p) & (RT.field == f)].rho_count.iloc[0]; assert abs(got - v) < 1e-9, (p, f, got, v)
for f, v in (("F", t10["F"]), ("L", t10["L"])):
    got = RT[(RT.partition == "K-means P2") & (RT.field == f) & (RT.realisation == "seed 42")].rho_count.iloc[0]; assert abs(got - v) < 1e-9, (f, got, v)
ab = pd.read_csv("revision/out/S3_ablation/P1A_ABLACION.csv")
for f, sc in (("F", "F"), ("L", "L_canon")):
    for p, col in (("Ward-4 (rook)", "rho_Ward_rook"), ("Ward-8 (queen)", "rho_Ward_queen"), ("SKATER (rook)", "rho_SKATER_rook")):
        got = RT[(RT.partition == p) & (RT.field == f) & (RT.realisation == "seed 42")].rho_count.iloc[0]
        ref = ab[(ab.block == "S1") & (ab.scenario == sc)][col].iloc[0]; assert abs(got - ref) < 1e-9, (f, p, got, ref)

# summary: deterministic fields as is; F and L as seed 42 and range over the 30 canonical realisations
def summ(g):
    d = g[g.realisation.str.startswith("seed 2")]
    s42 = g[g.realisation.isin(["seed 42", "deterministic"])].iloc[0]
    return pd.Series(dict(rho_count_s42=s42.rho_count, rho_density_s42=s42.rho_density, rho_field_area_s42=s42.rho_field_area,
                          rho_count_30_min=d.rho_count.min() if len(d) else np.nan, rho_count_30_max=d.rho_count.max() if len(d) else np.nan,
                          rho_density_30_min=d.rho_density.min() if len(d) else np.nan, rho_density_30_max=d.rho_density.max() if len(d) else np.nan))
SU = RT.groupby(["partition", "field"], sort=False).apply(summ, include_groups=False).reset_index()
SU.to_csv(OUT / "RHO_RESUMEN.csv", index=False)
area = ZT.groupby("partition", sort=False)[["rho_area_count", "rho_area_density"]].first()
area["zones_with_sites"] = ZT.groupby("partition", sort=False).sites.apply(lambda s: int((s > 0).sum()))
area["area_km2_min_max"] = ZT.groupby("partition", sort=False).area_km2.apply(lambda a: f"{a.min():.1f}-{a.max():.1f}")
txt = SU.round(3).to_string(index=False) + "\n\n" + area.round(3).to_string() + "\n\nControls against Table 10 and the ablation: OK\n"
(OUT / "diagnostico_log.txt").write_text(txt, encoding="utf-8"); print(txt)
