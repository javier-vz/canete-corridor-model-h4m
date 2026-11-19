"""Pre-tests of the final recompute (validation point C). Stops at the first failure.
  python code/test_final_pipeline.py            # everything except the v2 control runs
  python code/test_final_pipeline.py --ctrl     # compare the CTRL runs (final_runs.py --bloques CTRL) with the
                                                # published Tables 3-5 rows
Reference values are the published v8_11 numbers or frozen Phase B/B4/B5 outputs. When a published
number comes from v1 (provenance), the test says so; the gate is that the CODE reproduces it."""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import argparse, json, pickle, shutil, subprocess, tempfile, time
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

import ca_kernel_v2 as ca
import final_common as fc
import final_runs as fr
import final_evaluate as fe
import final_export_ppm as fx

LOG = []


def check(name, cond, detail=""):
    LOG.append(dict(test=name, ok=bool(cond), detail=str(detail)))
    print(f"[{'OK ' if cond else 'FAIL'}] {name} {detail}", flush=True)
    if not cond:
        Path("out/final/previa").mkdir(parents=True, exist_ok=True)
        pd.DataFrame(LOG).to_csv("out/final/previa/pruebas_previas.csv", index=False)
        raise SystemExit(f"STOP: {name}")


def note(name, detail):
    LOG.append(dict(test=name, ok="nota", detail=str(detail)))
    print(f"[NOTA] {name}: {detail}", flush=True)


def r(x, d):
    return round(float(x), d)


def main_tests():
    # ---------------------------------------------------------------- S0
    p = subprocess.run([_sys.executable, "code/test_final_h4m_state.py"], capture_output=True, text=True)
    check("S0 identity and designs (test_final_h4m_state.py)", p.returncode == 0 and "ALL OK" in p.stdout, p.stdout.strip().splitlines()[-1])
    st = fc.state_for(fr.cfg())
    st2 = pickle.load(open("out/state_v2_repro.pkl", "rb"))
    V = st.valid_mask
    # ---------------------------------------------------------------- S9 helper
    sv = fc.slope_variant_state(st, "S1_minmax")
    check("S9 helper: min-max variant = S0 (attractors, pairs)", sv.attractors == st.attractors and sv.pair_candidates == st.pair_candidates)
    check("S9 helper: E', c' of min-max variant = S0", max(np.abs(sv.E_grid - st.E_grid)[V].max(), np.abs(sv.c_grid - st.c_grid)[V].max()) < 1e-12)
    # ---------------------------------------------------------------- runner integrity (temp folder, results discarded)
    tmp = Path(tempfile.mkdtemp(prefix="jcaa80_test_runs_"))      # portable temporary folder (post-freeze change 1)
    sig, _ = fc.firma_calculo()
    fr._init(sig)
    J = {(j["block"], j["id"]): j for j in fr.jobs()}
    jf = J[("S7", "fete_c_n50_s2000")]
    fr.run_job(jf, tmp)
    check("runner: FETE job valid after writing", fr.validate(jf, sig, tmp)[0] is True)
    pj = fr.paths_for(jf, tmp)[0]; good = pj.read_text()
    d = json.loads(good); d["metrics"]["auc"] = float("nan"); pj.write_text(json.dumps(d))
    check("runner: NaN metric detected", fr.validate(jf, sig, tmp)[0] is False)
    d = json.loads(good); d["spec"]["seed"] = 2001; pj.write_text(json.dumps(d))
    check("runner: changed specification detected", fr.validate(jf, sig, tmp)[0] is False)
    pj.write_text(good)
    check("runner: other code/data signature detected", fr.validate(jf, "x" * 64, tmp)[0] is False)
    js = J[("S1", "canon_s42")]
    t0 = time.time(); o = fr.run_job(js, tmp)
    F, L = np.load(fr.paths_for(js, tmp)[1])["F"], np.load(fr.paths_for(js, tmp)[1])["L"]
    b4 = np.load("out/b4/runs/H4m_s42.npz")
    check("runner: S1 canonical = B4 H4m seed 42 (F, L, AUC)", np.abs(F - b4["F"]).max() == 0 and np.abs(L - b4["L"]).max() == 0 and o["metrics"]["auc"] == float(b4["auc"]),
          f"AUC {o['metrics']['auc']:.6f} ({time.time() - t0:.0f}s)")
    check("runner: observer F rebuilt = kernel F; 90 pairs, 360 paths", o["F_rebuild_maxdiff"] < 1e-12 and o["n_path_records"] == 90 and o["n_path_alternatives"] == 360,
          f"maxdiff {o['F_rebuild_maxdiff']:.1e}, records {o['n_path_records']}, paths {o['n_path_alternatives']}")
    check("runner: S1 job valid (npz and paths hashes)", fr.validate(js, sig, tmp)[0] is True)
    pz = fr.paths_for(js, tmp)[1]; raw = pz.read_bytes(); pz.write_bytes(raw[:-10] + b"0123456789")
    check("runner: altered npz detected", fr.validate(js, sig, tmp)[0] is False)
    shutil.rmtree(tmp)
    # ---------------------------------------------------------------- evaluation on v1/v2 fields
    L1, L2, F2 = np.load("out/L_v1_s42.npy"), np.load("out/L_v2_s42.npy"), np.load("out/F_v2_s42.npy")
    fc._PART.clear()
    t10 = {"1-h": (1 - st2.h_grid, 0.8496, 0.6063), "T": (st2.T_grid, 0.8177, 0.6533), "E": (st2.E_grid, 0.7563, 0.6871)}
    for k, (f, a, rr) in t10.items():
        check(f"Table 10 row {k} (P2) reproduced", r(fe.auc(f, st2), 4) == a and r(fc.rho20(f, st2), 4) == rr, f"{fe.auc(f, st2):.4f} / {fc.rho20(f, st2):.4f}")
    check("Table 10 row L v2: AUC 0.8078, rho20 P2 0.6376 (Phase A)", r(fe.auc(L2, st2), 4) == 0.8078 and r(fc.rho20(L2, st2), 4) == 0.6376)
    # Sec. 4.4 of v8_11 comes from the v1 field (provenance); the gate is that the zonal code reproduces it.
    z = {n: fc.zonal_stats(L1, fc.partition(st2, n)) for n in (20, 30, 40)}
    zp = fc.zone_prefix(L1, fc.partition(st2, 20)); k100 = int((zp.cum_sites >= 84).idxmax()) + 1
    got = dict(rho20=(r(z[20]["rho"], 3), r(z[20]["p"], 3)), rmax=r(z[20]["rho_max"], 3), zws=z[20]["zones_with_sites"],
               t1=(r(zp.cum_area_pct[0], 1), int(zp.cum_sites[0]), r(zp.lift[0], 1)), t3=(r(zp.cum_area_pct[2], 1), int(zp.cum_sites[2]), r(zp.lift[2], 1)),
               t7=(r(zp.cum_area_pct[6], 1), int(zp.cum_sites[6]), r(zp.lift[6], 1)), r30=(r(z[30]["rho"], 2), r(z[30]["p"], 3)),
               r40=(r(z[40]["rho"], 2), r(z[40]["p"], 3)), loo20=(r(z[20]["rho_loo"], 2), r(z[20]["p_loo"], 3)),
               loo30=(r(z[30]["rho_loo"], 2), r(z[30]["p_loo"], 3)))
    pub = dict(rho20=(0.638, 0.003), rmax=0.622, zws=9, t1=(8.0, 64, 9.5), t3=(23.7, 67, 3.4), t7=(47.1, 79, 2.0),
               r30=(0.45, 0.013), r40=(0.46, 0.003), loo20=(0.57, 0.011), loo30=(0.37, 0.047))
    def same(g, q):   # rho/area/sites/lift by rounding; p within 0.0006 (paper: 0.00249 -> 0.003, double rounding)
        if isinstance(q, tuple) and len(q) == 2 and isinstance(q[1], float) and q[1] < 0.1:
            return g[0] == q[0] and abs(g[1] - q[1]) <= 0.0006   # double rounding 0.00249 -> 0.0025 -> 0.003
        return g == q
    raw = dict(rho20=(r(z[20]["rho"], 3), z[20]["p"]), r30=(r(z[30]["rho"], 2), z[30]["p"]), r40=(r(z[40]["rho"], 2), z[40]["p"]),
               loo20=(r(z[20]["rho_loo"], 2), z[20]["p_loo"]), loo30=(r(z[30]["rho_loo"], 2), z[30]["p_loo"]))
    got.update(raw)
    bad = {k: (got[k], pub[k]) for k in pub if not same(got[k], pub[k])}
    check("Sec. 4.4 numbers reproduced by the zonal code (v1 L = provenance of v8_11)", not bad, str(bad) if bad else "")
    note("Sec. 4.4 '100% capture with the top 15 zones (79.2%)' not reproducible from the stored fields",
         f"v1 and v2 give {k100} zones ({zp.cum_area_pct[k100 - 1]:.1f}%); zones ranked 14 and 15 differ by "
         f"{abs(zp.L_mean[13] - zp.L_mean[14]):.1e} in mean L (near-tie): the published count depends on that order")
    zl = {nm: fc.zonal_stats(Lx, fc.partition(st2, 40)) for nm, Lx in (("v1", L1), ("v2", L2))}
    note("Sec. 4.4 'rho_loo = 0.40 (p = 0.011) at n = 40' not reproducible from the stored fields",
         f"v1 {zl['v1']['rho_loo']:.3f} (p {zl['v1']['p_loo']:.4f}); v2 {zl['v2']['rho_loo']:.3f} (p {zl['v2']['p_loo']:.4f}); same rule reproduces n = 20 and 30")
    REG = pd.read_csv("out/e3/regionalization_comparison.csv")
    Z = np.load("out/e3/partitions.npz"); cap2 = fe.Capture(st2)
    worst = 0
    for key in [k for k in Z.files if k not in ("iy", "ix")]:
        alg, n_ = key.rsplit("_n", 1)
        lab = np.empty(cap2.N, int); lab[cap2.perm_e3] = Z[key]
        part = pd.DataFrame({"iy": cap2.iy, "ix": cap2.ix, "zone": lab, "presencia_sitio_nuevo": cap2.sv})
        zs = fc.zonal_stats(L2, part); ref = REG[(REG.algorithm == alg) & (REG.n == int(n_))].iloc[0]
        worst = max(worst, abs(zs["rho"] - ref.rho_count), abs(zs["rho_loo"] - ref.rho_loo), abs(zs["p"] - ref.p_count))
    check("Table 12 code reproduces E3 metrics of the 12 partitions (v2 L)", worst < 1e-12, f"max diff {worst:.1e}")
    KS = pd.read_csv("out/e3/kmeans_100seeds_summary.csv").set_index("n")
    iyix, srm, km = np.load("out/rowmajor_iy_ix.npy"), np.load("out/rowmajor_sites.npy"), np.load("out/kmeans_labels_100.npy")
    for n in (20, 30, 40):
        labs = km if n == 20 else fe.kmeans_labels(st2, n)
        med = np.median(fe.zonal_rho_many(L2, labs, iyix, srm)[:, 0])
        check(f"100 K-means n={n}: median rho reproduces E3", abs(med - KS.loc[n, "rho_med"]) < 1e-12, f"{med:.4f}")
    CP = pd.read_csv("out/b4/t4/capture_per_seed.csv"); cap = fe.Capture(st)
    Lh = np.load("out/b4/runs/H4m_s42.npz")["L"]; got = cap.cells(Lh)
    ref = CP[(CP.predictor == "L") & (CP.scenario == "H4m") & (CP.seed.astype(str) == "42")].set_index("budget_pct")
    check("capture by cell reproduces B4 (H4m seed 42, 4 budgets)", all((got[b]["sites"], got[b]["tie_min"], got[b]["tie_max"]) ==
          tuple(ref.loc[int(b * 100), ["sites", "tie_min", "tie_max"]].astype(int)) for b in fe.BUD))
    CZ = pd.read_csv("out/b4/t4/capture_complete_zones.csv")
    got = cap2.zones(Z["KMeans_P2_rs42_n20"], L2)
    ref = CZ[(CZ.partition == "K-means P2") & (CZ.ranking_field == "L v2 seed 42")].set_index("budget_pct")
    check("capture by complete zones reproduces B4 (v2 seed 42)", all((got[b]["zones_used"], got[b]["sites"]) ==
          tuple(ref.loc[int(b * 100), ["zones_used", "sites"]].astype(int)) and abs(got[b]["area_consumed_pct"] - ref.loc[int(b * 100), "area_consumed_pct"]) < 1e-9 for b in fe.BUD))
    DM = pd.read_csv("out/b5/t5/decision_metrics.csv").set_index("scenario")
    for scen, folder in [("H4m", "out/b4/runs"), ("H4m-excl8", "out/b5/runs"), ("H4m-main", "out/b4/runs"), ("H1+H2m", "out/b4/runs")]:
        zz = [np.load(f"{folder}/{scen}_s{s}.npz") for s in (42, 43, 44)]
        a = np.mean([float(z["auc"]) for z in zz]); rr = np.mean([fc.rho20(z["L"], st) for z in zz])
        c10 = np.mean([100 * cap.cells(z["L"])[0.10]["sites"] / cap.NS for z in zz])
        check(f"Table 11 {scen}: AUC, rho20, capture10 = B5", max(abs(a - DM.loc[scen, "AUC_mean"]), abs(rr - DM.loc[scen, "rho20_P2_mean"]),
              abs(c10 - DM.loc[scen, "capture10_cell_pct_mean"]) / 100) < 1e-9)
    ov = {nm: [r(fe.overlap(Lx, st2, K, 1), 4) for K in (42, 168, 1260)] for nm, Lx in (("v1", L1), ("v2", L2))}
    check("buffered overlap code reproduces published 0.0714/0.0952/0.0968", [0.0714, 0.0952, 0.0968] in ov.values(), str(ov))
    lcc = {nm: fe.lcc_top(Lx, st2) for nm, Lx in (("v1", L1), ("v2", L2))}
    check("LCC code = kernel run_one for v2 seed 42 (90, cleanroom)", lcc["v2"] == 90, f"v1 {lcc['v1']} (published 69), v2 {lcc['v2']}")
    pmin, pmax, npairs = fe.cheapest_prob(pickle.load(open("out/b4/runs/v2_s42_paths.pkl", "rb")))
    check("cheapest-path probability: one value per routed pair, within [1/4, 1] (weighting = kernel, see F rebuild test)",
          npairs == 91 and 0.25 <= pmin <= pmax <= 1, f"v2 seed 42: {pmin:.4f}-{pmax:.4f}")
    note("Sec. 3.3 'probability between 0.2500 and 0.2566 across the 91 pairs' not reproducible from the v2 records",
         f"v2 perturbed costs (the kernel's weighting) give {pmin:.4f}-{pmax:.4f}; base costs give 0.2508-0.3190; "
         "v1 paths were not recorded, so the published range is probably from v1")
    res2 = st2.site_mask & ~fe.dominant_zone_sites(st2)
    k = fe.ks(st2.E_grid, st2, res2)
    check("KS function: exact method, residual set 20 sites (v2 E)", k["method_used"] == "exact" and k["n_sites"] == 20 and k["n_background"] == 12742 and r(k["D"], 6) == 0.433637,
          f"D {k['D']:.6f}, p {k['p']:.4e}")
    NC = pd.read_csv("out/normalization_L_comparison.csv").set_index("surface")
    R2 = np.load("out/river_cell_mask.npy")
    n5, _ = fe.normalization(F2, L2, st2, km, iyix, srm, R2)
    n5 = n5.set_index("surface")
    check("E5 normalisation code reproduces Phase B (v2)", max(np.abs(n5[c] - NC[c]).max() for c in ("AUC", "rho20_P2", "rho20_median_100", "spearman_vs_N1",
          "sites_in_top_2n_cells", "L_on_river_cells")) < 1e-12)
    pubcsv = pd.read_csv("data/ref/fields_canonical_v2_for_ppm.csv")   # published v2 PPM input, copied unchanged (post-freeze change 2)
    with tempfile.TemporaryDirectory(prefix="jcaa80_ppm_") as td:
        mine = fx.export(F2, L2, st2.E_grid, st2, Path(td) / "ppm_export_test.csv")
    check("PPM export reproduces the published v2 PPM input", max(np.abs(mine[c] - pubcsv[c]).max() for c in ("X", "Y", "F", "L", "E", "is_site")) < 1e-15)
    pp = pd.read_csv("out/final/previa/ppm_v2_reproduccion_R_aqui.csv")
    best = pp.loc[pp.groupby("covar").AIC.idxmin()].set_index("covar")
    check("R/spatstat here reproduces Table 7 (null 2867.8; F 157.9@3500; L 136.4@2750; E 110.0@2750)",
          r(pp.null_AIC.iloc[0], 1) == 2867.8 and (r(best.loc["F", "delta_AIC_vs_null"], 1), int(best.loc["F", "radius_m"])) == (157.9, 3500)
          and (r(best.loc["L", "delta_AIC_vs_null"], 1), int(best.loc["L", "radius_m"])) == (136.4, 2750)
          and (r(best.loc["E", "delta_AIC_vs_null"], 1), int(best.loc["E", "radius_m"])) == (110.0, 2750))
    from routing import build_graph, route
    G, idx = build_graph(st.c_grid, V, 250.0, False)
    dif = 0
    for (a, b) in st.pair_candidates[:5]:
        _, c1 = route(G, idx, tuple(st.attractors[a]), tuple(st.attractors[b]))
        c2 = nx.shortest_path_length(st.G_template, tuple(st.attractors[a]), tuple(st.attractors[b]), weight="weight_base")
        dif = max(dif, abs(c1 / 250.0 - c2) / c2)          # routing.py weights edges in metres, the kernel in cell units
    check("E2-A A1 graph = kernel graph (least-cost path cost / 250 m, 5 pairs)", dif < 1e-9, f"rel diff {dif:.1e}")


def ctrl_tests():
    sig, _ = fc.firma_calculo()
    J = [j for j in fr.jobs() if j["block"] == "CTRL"]
    for j in J:
        check(f"CTRL job valid: {j['id']}", fr.validate(j, sig, fc.OUT / "runs")[0] is True)
    pub = {"v2_a_w_m20": (0.8177, 0.0015, None), "v2_r0.80": (0.8118, 0.0028, 24), "v2_q0.975": (0.8200, 0.0014, 22)}
    for key, (m, s, na) in pub.items():
        rs = [json.loads(fr.paths_for(j, fc.OUT / "runs")[0].read_text()) for j in J if j["id"].startswith(key)]
        a = np.array([x["metrics"]["auc"] for x in rs])
        ok = r(a.mean(), 4) == m and r(a.std(ddof=1), 4) == s and (na is None or all(x["n_attractors"] == na for x in rs))
        check(f"CTRL {key}: published mean ± sd reproduced", ok, f"{a.mean():.4f} ± {a.std(ddof=1):.4f} (published {m} ± {s})")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--ctrl", action="store_true"); a = ap.parse_args()
    (ctrl_tests if a.ctrl else main_tests)()
    Path("out/final/previa").mkdir(parents=True, exist_ok=True)
    out = Path("out/final/previa") / ("pruebas_controles_v2.csv" if a.ctrl else "pruebas_previas.csv")
    pd.DataFrame(LOG).to_csv(out, index=False)
    print(f"ALL {sum(x['ok'] is True for x in LOG)} TESTS PASSED ({sum(x['ok'] == 'nota' for x in LOG)} notes) -> {out}")
