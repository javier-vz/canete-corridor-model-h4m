"""JCAA-80 final recompute: evaluation = the single numeric source of tables and text (frozen before
any final result exists; procedures in out/final/PROCEDIMIENTOS.md).

Reads: out/final/runs (final_runs.py, validated again here), the frozen hydrology sensitivities
(out/b4/runs, out/b5/runs, out/b5/t5/decision_metrics.csv, out/b5/t5/attractors_changed.csv), the
frozen partitions (out/e3/partitions.npz, out/kmeans_labels_100.npy, out/rowmajor_iy_ix.npy,
out/e3/regionalization_comparison.csv) and, when present, the S6 sweep (out/final/S6).
Writes: out/final/eval/*.csv and out/final/eval/indice.json (sha256 of every output).
Every pure metric below is regression-tested on v1/v2 fields in test_final_pipeline.py.
"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import importlib.util, json, pickle, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation, gaussian_filter, uniform_filter
from scipy.stats import ks_2samp, rankdata, spearmanr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

import ca_kernel_v2 as ca
import final_common as fc
import final_h4m_state as fs
import final_runs as fr

RUNS = fc.OUT / "runs"
EV = fc.EVAL
S6_DIR = fc.OUT / "S6"
BUD = (0.05, 0.10, 0.20, 0.30)
KMULT = (0.5, 1, 2, 5, 10, 15)


# =============================================================================== pure metrics
def auc(field, st):
    return float(ca.auc_score(field, st.site_mask, st.valid_mask))


def dominant_zone_sites(st):
    part = fc.partition(st, 20)
    z = fc.zone_frame(np.zeros(st.valid_mask.shape), part)
    dom = z.n_sites.idxmax()
    m = np.zeros(st.valid_mask.shape, bool)
    sel = (part.zone.values == dom) & (part.presencia_sitio_nuevo.values == 1)
    m[part.iy.values[sel], part.ix.values[sel]] = True
    return m


def ks(field, st, sites_mask):
    """PROCEDIMIENTOS 6: two-sample, two-sided, exact; background = valid cells without sites."""
    bg = st.valid_mask & ~st.site_mask
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        r = ks_2samp(field[sites_mask], field[bg], alternative="two-sided", method="exact")
    fb = any("exact" in str(x.message).lower() for x in w)
    return dict(D=float(r.statistic), p=float(r.pvalue), method_used="asymp (fallback)" if fb else "exact",
                n_sites=int(sites_mask.sum()), n_background=int(bg.sum()))


def overlap(field, st, K, b):
    V, S = st.valid_mask, st.site_mask
    thr = np.sort(field[V])[::-1][K - 1]
    T = V & (field >= thr)
    B = binary_dilation(S, structure=np.ones((2 * b + 1, 2 * b + 1))) if b > 0 else S
    return float((T & B).sum() / T.sum())


def lcc_top(field, st):
    S_top, _ = ca.threshold_topK(field, st.valid_mask, K=2 * st.n_sites)
    sizes, _ = ca.cluster_sizes(S_top, st.valid_mask)
    return int(sizes.max())


def cheapest_prob(records, beta=1.5):
    pmax = []
    for r in records:
        c = np.asarray(r["costs"], float)
        p = np.exp(-beta * (c - c.min())); p /= p.sum()
        pmax.append(p.max())
    return float(min(pmax)), float(max(pmax)), len(pmax)


class Capture:
    """PROCEDIMIENTOS 4 (identical rules to capture_no_interpolation.py, Phase B4)."""

    def __init__(self, st):
        self.iy, self.ix = np.nonzero(st.valid_mask)
        self.N = len(self.iy)
        self.sv = st.site_mask[self.iy, self.ix].astype(int)
        self.NS = int(self.sv.sum())
        self.TIE = np.random.default_rng(0).permutation(self.N)
        self.NB = {b: int(round(b * self.N)) for b in BUD}
        Z = np.load("out/e3/partitions.npz")
        pos = {(a, b): k for k, (a, b) in enumerate(zip(self.iy, self.ix))}
        self.perm_e3 = np.array([pos[(a, b)] for a, b in zip(Z["iy"], Z["ix"])])
        rm = np.load("out/rowmajor_iy_ix.npy")
        self.perm_rm = np.array([pos[(a, b)] for a, b in rm])

    def cells(self, field):
        v = field[self.iy, self.ix]
        order = np.lexsort((self.TIE, -v))
        out = {}
        for b, n in self.NB.items():
            vstar = v[order[n - 1]]
            above, tied = v > vstar, v == vstar
            a, t = int(above.sum()), int(tied.sum()); sa, stt = int(self.sv[above].sum()), int(self.sv[tied].sum())
            take = n - a
            out[b] = dict(n_cells=n, sites=int(self.sv[order[:n]].sum()), tie_block_cells=t,
                          tie_min=sa + max(0, take - (t - stt)), tie_max=sa + min(stt, take))
        return out

    def zones(self, lab_e3, field):
        lab = np.empty(self.N, int); lab[self.perm_e3] = lab_e3
        g = pd.DataFrame({"z": lab, "v": field[self.iy, self.ix], "s": self.sv}).groupby("z").agg(
            v=("v", "mean"), s=("s", "sum"), n=("v", "size")).sort_values("v", ascending=False)
        cn, cs = g.n.cumsum().values, g.s.cumsum().values
        out = {}
        for b in BUD:
            k = int(np.searchsorted(cn, b * self.N, side="right"))
            out[b] = dict(zones_used=k, area_consumed_pct=100 * (cn[k - 1] if k else 0) / self.N,
                          sites=int(cs[k - 1]) if k else 0, sites_pct=100 * (cs[k - 1] if k else 0) / self.NS,
                          first_zone_not_fitting_area_pct=100 * g.n.values[k] / self.N if k < len(g) else np.nan)
        return out


def zonal_rho_many(field, labels_rm, iyix, sites):
    v = field[iyix[:, 0], iyix[:, 1]]
    out = []
    for lab in labels_rm:
        g = pd.DataFrame({"z": lab, "v": v, "s": sites}).groupby("z")
        z = pd.DataFrame({"m": g.v.mean(), "s": g.s.sum()})
        r = spearmanr(z.m, z.s)[0]
        loo = z.drop(z.s.idxmax())
        out.append((r, spearmanr(loo.m, loo.s)[0]))
    return np.array(out)


def kmeans_labels(st, n, seeds=range(100)):
    """K-means initialisations (n zones) with the P2 features, row-major order (= rowmajor_iy_ix)."""
    part = fc.partition(st, 20)
    X = StandardScaler().fit_transform(part[["ix", "iy"]].assign(h=part.h_norm, p=part.p_norm).to_numpy())
    return np.array([KMeans(n_clusters=n, random_state=rs, n_init=10).fit(X).labels_ for rs in seeds])


def mm(a, V):
    o = np.zeros_like(a); v = a[V]; o[V] = (v - v.min()) / (v.max() - v.min()); return o


def prank(a, V):
    o = np.zeros_like(a); o[V] = (rankdata(a[V]) - 1) / (V.sum() - 1); return o


def normalization(F, L, st, labels100, iyix, sites, R):
    """E5 under the final fields (same definitions as e5_normalization.py)."""
    V = st.valid_mask
    F0 = np.where(V, F, 0); C = mm(gaussian_filter(F0, 1.0, mode="nearest"), V); D = mm(uniform_filter(F0, 5, mode="nearest"), V); T = st.T_grid
    assert np.abs(mm(0.55 * F + 0.2 * T + 0.15 * C + 0.1 * D, V) - L)[V].max() < 1e-9, "L is not the kernel combination"
    S = {"N1_canonical": L, "N2_F_minmax": mm(0.55 * mm(F, V) + 0.2 * T + 0.15 * C + 0.1 * D, V),
         "N3_all_percentile_rank": mm(0.55 * prank(F, V) + 0.2 * prank(T, V) + 0.15 * prank(C, V) + 0.1 * prank(D, V), V)}
    rows, comp = [], []
    K = 2 * st.n_sites
    for name, Lx in S.items():
        rr = zonal_rho_many(Lx, labels100, iyix, sites)[:, 0]
        top = Lx[V] >= np.sort(Lx[V])[-K]
        rows.append(dict(surface=name, AUC=auc(Lx, st), rho20_P2=fc.rho20(Lx, st), rho20_median_100=np.median(rr),
                         rho20_p05=np.percentile(rr, 5), rho20_p95=np.percentile(rr, 95),
                         spearman_vs_N1=spearmanr(Lx[V], L[V])[0], sites_in_top_2n_cells=int(st.site_mask[V][top].sum()),
                         L_on_river_cells=float(Lx[R & V].sum() / Lx[V].sum())))
    for name, (fF, fT, fC, fD) in {"N1/N2 kernel scales": (F, T, C, D), "N2 F min-max": (mm(F, V), T, C, D),
                                    "N3 percentile ranks": (prank(F, V), prank(T, V), prank(C, V), prank(D, V))}.items():
        for lab, w, a in [("F", .55, fF), ("T", .2, fT), ("C", .15, fC), ("D", .1, fD)]:
            x = (w * a)[V]; comp.append(dict(scaling=name, term=f"{w}*{lab}", mean=x.mean(), max=x.max(), var=x.var()))
    cd = pd.DataFrame(comp)
    cd["share_of_mean_sum"] = cd.groupby("scaling")["mean"].transform(lambda m: m / m.sum())
    return pd.DataFrame(rows), cd


# =============================================================================== readers
JOBS = {(j["block"], j["id"]): j for j in fr.jobs()}


def check_block(block):
    sig, _ = fc.firma_calculo()
    J = [j for j in JOBS.values() if j["block"] == block]
    bad = [(j["id"], why) for j in J for ok, why in [fr.validate(j, sig, RUNS)] if not ok]
    if bad:
        raise SystemExit(f"block {block}: {len(bad)} invalid or missing jobs, e.g. {bad[:3]}")
    return J


def jget(block, jid):
    return json.loads(fr.paths_for(JOBS[(block, jid)], RUNS)[0].read_text())


def fget(block, jid):
    z = np.load(fr.paths_for(JOBS[(block, jid)], RUNS)[1])
    return z["F"], z["L"]


def recget(block, jid):
    return pickle.load(open(fr.paths_for(JOBS[(block, jid)], RUNS)[2], "rb"))


def msd(x):
    x = np.asarray(x, float)
    return x.mean(), x.std(ddof=1) if len(x) > 1 else np.nan


# =============================================================================== tables
def sensitivity_table(block, st0, label_of):
    J = check_block(block)
    rows = []
    for j in J:
        r = jget(block, j["id"]); _, L = fget(block, j["id"])
        rows.append(dict(config=label_of(j), seed=j["seed"], attractors=r["n_attractors"], pairs=r["n_pairs"],
                         auc=r["metrics"]["auc"], rho20=fc.rho20(L, st0)))
    D = pd.DataFrame(rows)
    T = D.groupby("config", sort=False).agg(attractors=("attractors", "first"), attr_unique=("attractors", "nunique"),
                                            pairs=("pairs", "first"), n_seeds=("seed", "size"),
                                            auc_mean=("auc", "mean"), auc_sd=("auc", "std"),
                                            rho20_mean=("rho20", "mean"), rho20_sd=("rho20", "std")).reset_index()
    assert (T.attr_unique == 1).all()
    return T.drop(columns="attr_unique"), D


def main():
    EV.mkdir(parents=True, exist_ok=True)
    out = {}
    st = fc.state_for(fr.cfg())                       # final H4m state (S0)
    st_E = fc.state_for(fr.cfg(rules="v2"))           # only for E without H2m (Eq. E) in Table 10
    V, S = st.valid_mask, st.site_mask
    R, Rm, d0 = fc.masks_and_dist(st)
    cap = Capture(st)
    iyix = np.load("out/rowmajor_iy_ix.npy"); sites_rm = np.load("out/rowmajor_sites.npy")
    km100 = np.load("out/kmeans_labels_100.npy")

    # ---------------- S1
    check_block("S1")
    Fc, Lc = fget("S1", "canon_s42"); Fd, Ld = fget("S1", "det_s42")
    jc, jd = jget("S1", "canon_s42"), jget("S1", "det_s42")
    rec = recget("S1", "canon_s42")
    pmin, pmax, npairs = cheapest_prob(rec)

    # ---------------- Tables 3-5
    lab3 = {fr.safe(l): l for l, _, _ in fs.weight_configs()}
    T3, D3 = sensitivity_table("S3", st, lambda j: lab3[j["id"].rsplit("_s", 1)[0]])
    lab4 = {fr.safe(l): l for l, *_ in fs.T_configs()}
    T4, D4 = sensitivity_table("S4", st, lambda j: lab4[j["id"].rsplit("_s", 1)[0]])
    T5, D5 = sensitivity_table("S5", st, lambda j: f"q={j['cfg']['q']}")
    out.update(tabla3_pesos=T3, tabla4_T=T4, tabla5_cuantil=T5,
               sensibilidades_por_semilla=pd.concat([D3.assign(bloque="S3"), D4.assign(bloque="S4"), D5.assign(bloque="S5")]))

    # ---------------- Table 6 (KS)
    res_mask = S & ~dominant_zone_sites(st)
    T6 = pd.DataFrame([dict(site_set=nm, field=f, **ks(fld, st, msk))
                       for nm, msk in [("all 84 sites", S), ("residual (outside dominant P2 zone)", res_mask)]
                       for f, fld in [("E'", st.E_grid), ("L", Lc)]])
    out["tabla6_ks"] = T6

    # ---------------- Tables 8-9 (FETE)
    check_block("S2"); check_block("S7"); check_block("S8")
    s2 = pd.DataFrame([jget("S2", f"canon_s{s}")["metrics"] for s in range(2000, 2030)])
    s7 = pd.DataFrame([dict(n=n, **jget("S7", f"fete_c_n{n}_s{s}")["metrics"]) for n in (50, 90, 200, 500, 1000) for s in range(2000, 2030)])
    s8 = pd.DataFrame([jget("S8", f"fete_1mE_n90_s{s}")["metrics"] for s in range(2000, 2030)])

    def summ(df, label, n):
        a, b = msd(df.auc); c, d = msd(df.lcc_size); e, f = msd(df.m2)
        ok = bool((df.n_successful_paths == df.n_random_pairs).all()) if "n_successful_paths" in df else True
        return dict(method=label, n_pairs=n, n_seeds=len(df), auc_mean=a, auc_sd=b, lcc_mean=c, lcc_sd=d, m2_mean=e, m2_sd=f, all_paths_found=ok)
    out["tabla8_fete"] = pd.DataFrame([summ(s2, "corridor model (H4m)", 90), summ(s7[s7.n == 90], "FETE, cost c'", 90),
                                       summ(s8, "FETE, cost 1 - E'", 90)])
    out["tabla9_fete_presupuestos"] = pd.DataFrame([summ(s7[s7.n == n], "FETE, cost c'", n) for n in (50, 90, 200, 500, 1000)])

    # ---------------- Table 10
    preds = {"1 - h": 1 - st.h_grid, "T": st.T_grid, "E (Eq. E, without H2m)": st_E.E_grid, "E' (final affinity)": st.E_grid,
             "F": Fc, "L": Lc}
    rows = []
    for nm, f in preds.items():
        zs = fc.zonal_stats(f, fc.partition(st, 20))
        rows.append(dict(predictor=nm, auc=auc(f, st), rho20=zs["rho"], p20=zs["p"]))
    out["tabla10_predictores"] = pd.DataFrame(rows)

    # ---------------- zonal (Sec. 4.4), Fig 8-9
    zrows, prefix = [], []
    for n in (20, 30, 40):
        part = fc.partition(st, n)
        zrows.append(dict(n=n, **fc.zonal_stats(Lc, part)))
        prefix.append(fc.zone_prefix(Lc, part).assign(n=n))
    out["zonal_P2"] = pd.DataFrame(zrows)
    out["zonas_prefijo"] = pd.concat(prefix)

    # ---------------- Table 12 (partitions) and 100 K-means
    Z = np.load("out/e3/partitions.npz")
    spatial = pd.read_csv("out/e3/regionalization_comparison.csv")
    prow = []
    for key in [k for k in Z.files if k not in ("iy", "ix")]:
        alg, n_ = key.rsplit("_n", 1)
        lab_rm = np.empty(cap.N, int); lab_rm[cap.perm_e3] = Z[key]
        part = pd.DataFrame({"iy": cap.iy, "ix": cap.ix, "zone": lab_rm, "presencia_sitio_nuevo": cap.sv})
        zs = fc.zonal_stats(Lc, part); zp = fc.zone_prefix(Lc, part)
        prow.append(dict(algorithm=alg, n=int(n_), rho=zs["rho"], p=zs["p"], rho_loo=zs["rho_loo"], p_loo=zs["p_loo"],
                         zones_with_sites=zs["zones_with_sites"], top1_sites=int(zp.cum_sites.iloc[0]),
                         top1_area_pct=zp.cum_area_pct.iloc[0], top3_sites=int(zp.cum_sites.iloc[2]), top3_area_pct=zp.cum_area_pct.iloc[2]))
    T12 = pd.DataFrame(prow)
    keep = ["fragments_total", "zones_fragmented", "max_fragments", "polsby_popper_median", "size_min", "size_max",
            "largest_zone_area_share", "within_sd_alt_m", "within_sd_slope_deg"]
    T12 = T12.merge(spatial[["algorithm", "n"] + keep], on=["algorithm", "n"], how="left", validate="one_to_one")
    assert T12[keep].notna().all().all(), "spatial metrics missing for some partition"
    out["tabla12_particiones"] = T12
    krows = []
    for n in (20, 30, 40):
        labs = km100 if n == 20 else kmeans_labels(st, n)
        rr = zonal_rho_many(Lc, labs, iyix, sites_rm)
        krows.append(dict(n=n, n_inits=len(labs), rho_median=np.median(rr[:, 0]), rho_p05=np.percentile(rr[:, 0], 5),
                          rho_p95=np.percentile(rr[:, 0], 95), rho_loo_median=np.median(rr[:, 1]), rho_loo_p05=np.percentile(rr[:, 1], 5)))
    out["kmeans100_resumen"] = pd.DataFrame(krows)

    # ---------------- capture (decision 2)
    Ls2 = {s: fget("S2", f"canon_s{s}")[1] for s in range(2000, 2030)}
    Lmean = np.mean(list(Ls2.values()), axis=0)
    crow = []
    for s, f in Ls2.items():
        for b, r in cap.cells(f).items():
            crow.append(dict(predictor="L", field=f"seed {s}", budget_pct=int(b * 100), **r))
    for nm, f in {"L seed 42": Lc, "L mean of S2 (30 seeds)": Lmean, "1 - h": 1 - st.h_grid, "T": st.T_grid, "E'": st.E_grid}.items():
        for b, r in cap.cells(f).items():
            crow.append(dict(predictor=nm.split()[0] if nm.startswith("L") else nm, field=nm, budget_pct=int(b * 100), **r))
    C = pd.DataFrame(crow)
    for c in ("sites", "tie_min", "tie_max"):
        C[c + "_pct"] = 100 * C[c] / cap.NS
    out["captura_celda"] = C
    ps = C[C.field.str.startswith("seed")]
    out["captura_celda_resumen"] = ps.groupby("budget_pct").agg(n_seeds=("field", "nunique"), sites_pct_mean=("sites_pct", "mean"),
                                                               sites_pct_sd=("sites_pct", "std"), sites_pct_min=("sites_pct", "min"),
                                                               sites_pct_max=("sites_pct", "max"),
                                                               tie_width_max_pct=("tie_max_pct", lambda x: (x - ps.loc[x.index, "tie_min_pct"]).max())).reset_index()
    zc = []
    for key, lab in [("K-means P2", "KMeans_P2_rs42_n20"), ("Ward 4-neighbour", "Ward_rook_n20"), ("Ward 8-neighbour", "Ward_queen_n20"),
                     ("SKATER 4-neighbour", "SKATER_rook_n20")]:
        for nm, f in {"L seed 42": Lc, "L mean of S2 (30 seeds)": Lmean, "1 - h": 1 - st.h_grid, "T": st.T_grid, "E'": st.E_grid}.items():
            for b, r in cap.zones(Z[lab], f).items():
                zc.append(dict(partition=key, ranking_field=nm, budget_pct=int(b * 100), **r))
    out["captura_zonas"] = pd.DataFrame(zc)

    # ---------------- Table 11 (hydrology sensitivities, frozen runs)
    DM = pd.read_csv("out/b5/t5/decision_metrics.csv").set_index("scenario")
    hrows = []
    for scen, folder in [("H4m", "out/b4/runs"), ("H4m-excl8", "out/b5/runs"), ("H4m-main", "out/b4/runs"), ("H1+H2m", "out/b4/runs")]:
        for s in (42, 43, 44):
            z = np.load(f"{folder}/{scen}_s{s}.npz")
            hrows.append(dict(scenario=scen, seed=s, auc=float(z["auc"]), rho20=fc.rho20(z["L"], st),
                              capture10_pct=100 * cap.cells(z["L"])[0.10]["sites"] / cap.NS))
    H = pd.DataFrame(hrows)
    T11 = H.groupby("scenario", sort=False).agg(auc_mean=("auc", "mean"), auc_min=("auc", "min"), auc_max=("auc", "max"),
                                                rho20_mean=("rho20", "mean"), capture10_mean=("capture10_pct", "mean"),
                                                capture10_min=("capture10_pct", "min"), capture10_max=("capture10_pct", "max")).reset_index()
    for c in ("attractors_mean", "pairs_mean", "F_main_pct_mean", "F_secondary_pct_mean", "F_8lines_cells_pct_mean",
              "F_other_secondary_pct_mean", "C3_angular_pair_equal_mean", "C3_angular_traffic_mean",
              "AUC_mean", "rho20_P2_mean", "capture10_cell_pct_mean"):
        T11[c + " (B5)"] = T11.scenario.map(DM[c])
    T11["check_vs_B5_maxdiff"] = np.nanmax(np.abs(np.c_[T11.auc_mean - T11["AUC_mean (B5)"], T11.rho20_mean - T11["rho20_P2_mean (B5)"],
                                                        (T11.capture10_mean - T11["capture10_cell_pct_mean (B5)"]) / 100]), axis=1)
    out["tabla11_hidrologia"] = T11
    out["hidrologia_por_semilla"] = H

    # ---------------- S9 slope variants under H4m (E4)
    check_block("S9")
    s9 = []
    for v in fr.SLOPE_VARIANTS:
        stv = fc.state_for(fr.cfg(slope=v)); F9, L9 = fget("S9", f"{v}_s42"); r9 = jget("S9", f"{v}_s42")
        rr = zonal_rho_many(L9, km100, iyix, sites_rm)[:, 0]
        s9.append(dict(variant=v, attractors=len(stv.attractors), pairs=len(stv.pair_candidates),
                       attractors_shared_with_S1=len(set(map(tuple, stv.attractors)) & set(map(tuple, st.attractors))),
                       spearman_c_vs_S1=spearmanr(stv.c_grid[V], st.c_grid[V])[0], spearman_F_vs_S1=spearmanr(F9[V], Fc[V])[0],
                       spearman_L_vs_S1=spearmanr(L9[V], Lc[V])[0], F_on_river_cells_pct=100 * F9[R & V].sum() / F9[V].sum(),
                       auc=r9["metrics"]["auc"], rho20_P2=fc.rho20(L9, st), rho20_median_100=np.median(rr)))
    s9.append(dict(variant="S1 (min-max, final)", attractors=len(st.attractors), pairs=len(st.pair_candidates),
                   attractors_shared_with_S1=len(st.attractors), spearman_c_vs_S1=1.0, spearman_F_vs_S1=1.0, spearman_L_vs_S1=1.0,
                   F_on_river_cells_pct=100 * Fc[R & V].sum() / Fc[V].sum(), auc=jc["metrics"]["auc"], rho20_P2=fc.rho20(Lc, st),
                   rho20_median_100=np.median(zonal_rho_many(Lc, km100, iyix, sites_rm)[:, 0])))
    out["pendiente_E4_H4m"] = pd.DataFrame(s9)

    # ---------------- E5 normalisation and L decomposition
    N5, C5 = normalization(Fc, Lc, st, km100, iyix, sites_rm, R)
    out["normalizacion_E5"] = N5; out["descomposicion_L"] = C5

    # ---------------- overlap
    out["solapamiento"] = pd.DataFrame([dict(b=b, K=int(k * st.n_sites), K_mult=k, overlap=overlap(Lc, st, int(k * st.n_sites), b))
                                        for b in (0, 1, 2) for k in KMULT])

    # ---------------- S6 sweep (if available)
    sweep = sweep_summary()
    if sweep is not None:
        out["barrido_S6"] = sweep

    # ---------------- numbers for the text (single source)
    A = pd.DataFrame(st.attractors, columns=["iy", "ix"]); A["d"] = d0[A.iy, A.ix]
    ch = pd.read_csv("out/b5/t5/attractors_changed.csv")
    import geopandas as gpd
    from shapely.geometry import Point
    cl = gpd.read_file("out/hydro/hydrorivers_model_lines_classified.gpkg")             # same loading as b5_exclusion_runs
    fid = pd.read_csv("out/b4/t1/hydrorivers_official_match.csv").query("dataset == 'C87'").set_index("line").shp_fid
    cl["shp_fid"] = [int(fid[i]) for i in range(len(cl))]
    EXCL_FID = [1131506, 1130585, 1129684, 1130586, 1127717, 1130415, 1126308, 1126640]   # = b5_exclusion_runs.EXCL_FID
    assert set(EXCL_FID) <= set(cl.shp_fid)
    L8 = cl[cl.shp_fid.isin(EXCL_FID)].geometry.union_all()
    ix2x = {v: k for k, v in st.lattice["x_to_ix"].items()}; iy2y = {v: k for k, v in st.lattice["y_to_iy"].items()}
    d8 = np.array([Point(ix2x[j], iy2y[i]).distance(L8) for i, j in st.attractors])
    zr = out["zonal_P2"].set_index("n"); zp20 = out["zonas_prefijo"][out["zonas_prefijo"].n == 20].reset_index(drop=True)
    k100 = int((zp20.cum_sites >= cap.NS).idxmax()) + 1
    t3, t4, t5 = out["tabla3_pesos"], out["tabla4_T"], out["tabla5_cuantil"]
    num = [
        ("X1", "attractors (S0)", len(st.attractors)), ("X5", "attractor pairs (S0)", len(st.pair_candidates)),
        ("X3a", "river cells R (H1 exclusion)", int((R & V).sum())), ("X3b", "river cells, % of Omega", 100 * (R & V).sum() / V.sum()),
        ("X3c", "sites in river cells", int((R & S).sum())), ("X3d", "sites in river cells, % of 84", 100 * (R & S).sum() / S.sum()),
        ("X3e", "attractor distance to river, min (m)", A.d.min()), ("X3f", "attractor distance to river, median (m)", A.d.median()),
        ("X3g", "attractor distance to river, max (m)", A.d.max()),
        ("X3h", "attractors shared with v2", len(set(map(tuple, st.attractors)) & set(map(tuple, pickle.load(open("out/state_v2_repro.pkl", "rb")).attractors)))),
        ("X6a", "cheapest-path probability, min over pairs (seed 42)", pmin), ("X6b", "cheapest-path probability, max over pairs", pmax),
        ("X6c", "pairs routed (seed 42)", npairs), ("X16", "path records (seed 42, canonical)", jc["n_path_alternatives"]),
        ("X16b", "path records (deterministic)", jd["n_path_alternatives"]),
        ("X8a", "weights: mean AUC min", t3.auc_mean.min()), ("X8b", "weights: mean AUC max", t3.auc_mean.max()),
        ("X8c", "weights: mean rho20 min", t3.rho20_mean.min()), ("X8d", "weights: mean rho20 max", t3.rho20_mean.max()),
        ("X8e", "T: mean AUC min", t4.auc_mean.min()), ("X8f", "T: mean AUC max", t4.auc_mean.max()),
        ("X8g", "T: mean rho20 min", t4.rho20_mean.min()), ("X8h", "T: mean rho20 max", t4.rho20_mean.max()),
        ("X8i", "T: attractors min", t4.attractors.min()), ("X8j", "T: attractors max", t4.attractors.max()),
        ("X8k", "quantile: attractors min", t5.attractors.min()), ("X8l", "quantile: attractors max", t5.attractors.max()),
        ("X8m", "quantile: pairs min", t5.pairs.min()), ("X8n", "quantile: pairs max", t5.pairs.max()),
        ("X8o", "quantile: mean AUC min", t5.auc_mean.min()), ("X8p", "quantile: mean AUC max", t5.auc_mean.max()),
        ("X9", "rho20 seed 42 (P2)", zr.loc[20, "rho"]), ("X12", "median run time S2 here (s)", float(np.median([jget("S2", f"canon_s{s}")["runtime_s"] for s in range(2000, 2030)]))),
        ("X15a", "LCC top-168 L, deterministic", jd["metrics"]["lcc_size"]), ("X15b", "LCC top-168 L, canonical seed 42", jc["metrics"]["lcc_size"]),
        ("X17a", "AUC deterministic", jd["metrics"]["auc"]), ("X17b", "AUC canonical seed 42", jc["metrics"]["auc"]),
        ("X17c", "AUC S2 mean (30 seeds)", s2.auc.mean()), ("X17d", "AUC S2 sd", s2.auc.std(ddof=1)),
        ("X18a", "overlap b=1, K=42", overlap(Lc, st, 42, 1)), ("X18b", "overlap b=1, K=168", overlap(Lc, st, 168, 1)),
        ("X18c", "overlap b=1, K=1260", overlap(Lc, st, 1260, 1)),
        ("X20a", "zones with sites (n=20)", zr.loc[20, "zones_with_sites"]), ("X20b", "rho n=20", zr.loc[20, "rho"]), ("X20c", "p n=20", zr.loc[20, "p"]),
        ("X20d", "rho Lmax n=20", zr.loc[20, "rho_max"]), ("X20e", "p Lmax n=20", zr.loc[20, "p_max"]),
        ("X22a", "rho n=30", zr.loc[30, "rho"]), ("X22b", "p n=30", zr.loc[30, "p"]), ("X22c", "rho n=40", zr.loc[40, "rho"]), ("X22d", "p n=40", zr.loc[40, "p"]),
        ("X22e", "rho_loo n=20", zr.loc[20, "rho_loo"]), ("X22f", "p_loo n=20", zr.loc[20, "p_loo"]),
        ("X22g", "rho_loo n=30", zr.loc[30, "rho_loo"]), ("X22h", "p_loo n=30", zr.loc[30, "p_loo"]),
        ("X22i", "rho_loo n=40", zr.loc[40, "rho_loo"]), ("X22j", "p_loo n=40", zr.loc[40, "p_loo"]),
        ("X21a", "top-1 zone area %", zp20.cum_area_pct[0]), ("X21b", "top-1 zone sites", zp20.cum_sites[0]), ("X21c", "top-1 lift", zp20.lift[0]),
        ("X21d", "top-3 zones area %", zp20.cum_area_pct[2]), ("X21e", "top-3 sites", zp20.cum_sites[2]), ("X21f", "top-3 lift", zp20.lift[2]),
        ("X21g", "top-7 zones area %", zp20.cum_area_pct[6]), ("X21h", "top-7 sites", zp20.cum_sites[6]), ("X21i", "top-7 lift", zp20.lift[6]),
        ("X21j", "zones needed for 100 % of sites", k100), ("X21k", "area of those zones %", zp20.cum_area_pct[k100 - 1]),
        ("X27a", "H4m attractors within 1 km of the 8 priority lines", int((d8 <= 1000).sum())),
        ("X27b", "their distance to the 8 lines, min (m)", float(d8[d8 <= 1000].min())),
        ("X27c", "their distance to the 8 lines, max (m)", float(d8[d8 <= 1000].max())),
        ("X27d", "attractors that disappear when the 8 lines are excluded (B5)", int((ch.change == "only in H4m").sum())),
        ("X27e", "attractors that appear when the 8 lines are excluded (B5)", int((ch.change == "only in H4m-excl8").sum())),
    ]
    if sweep is not None:
        num += [("X7a", "sweep: min of configuration-mean AUC", sweep.auc_media.min()), ("X7b", "sweep: max of configuration-mean AUC", sweep.auc_media.max())]
    out["cifras_texto"] = pd.DataFrame(num, columns=["id", "descripcion", "valor"])

    # ---------------- write + index
    idx = {}
    for name, df in out.items():
        p = EV / f"{name}.csv"; df.to_csv(p, index=False); idx[p.name] = fc.sha256(p)
    sig, parts = fc.firma_calculo()
    json.dump(dict(firma_calculo=sig, archivos=idx, procedimientos_sha256=fc.sha256(fc.OUT / "PROCEDIMIENTOS.md"),
                   evaluate_sha256=fc.sha256("code/final_evaluate.py")), open(EV / "indice.json", "w"), indent=1)
    print(f"{len(idx)} tables written to {EV}")


def sweep_summary():
    """S6: validate the 1,210 JSON with the sweep script's own rules; None if not yet available."""
    res = S6_DIR / "resultados_barrido_H4m"
    if not res.exists():
        return None
    spec = importlib.util.spec_from_file_location("sweep", "out/b6/ryzen/barrido_H4m/final_sweep_eps_kappa.py")
    sw = importlib.util.module_from_spec(spec); spec.loader.exec_module(sw)
    sig = json.load(open(S6_DIR / "prueba_OK.json"))["firma_paquete"]
    jobs = sw.all_jobs(res)
    rows = []
    for e, k, s, p in jobs:
        ok, why = sw.validate(p, e, k, s, sig)
        if not ok:
            raise SystemExit(f"S6: invalid {Path(p).name}: {why}")
        rows.append(json.loads(Path(p).read_text()))
    D = pd.DataFrame(rows)
    assert len(D) == 1210 and not D.duplicated(["eps", "kappa", "seed"]).any()
    extra = [f.name for f in res.glob("*.json") if f.name not in {Path(j[3]).name for j in jobs}]
    assert not extra, f"S6: unexpected files {extra[:3]}"
    return D.groupby(["eps", "kappa"]).agg(n=("auc", "size"), auc_media=("auc", "mean"), auc_de=("auc", "std"),
                                           rho20_media=("rho20_P2", "mean"), rho20_de=("rho20_P2", "std"),
                                           lcc_media=("lcc_size", "mean"), m2_media=("m2", "mean")).reset_index()


if __name__ == "__main__":
    main()
