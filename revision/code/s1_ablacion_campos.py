"""JCAA-80, audit after Stage 6, P1-A: diagnostic ablation of the fixation field L on the STORED fields (no model runs,
no change to H4m, to any frozen file or to the manuscript).

Reads (read-only, from the repository root): the stored fields F and L of the canonical runs
(out/final/runs/S1/canon_s42.npz and out/final/runs/S2/canon_s2000.npz ... canon_s2029.npz), the frozen H4m state
(final_common.state_for(final_runs.cfg()): T, h, valid cells, site cells), the reference K-means partition P2
(final_common.partition) and the stored contiguity-constrained partitions (out/e3/partitions.npz).

Scenarios (all min-max normalisations over the valid cells Omega, as in ca_kernel_v2.build_fixation_field):
  L_canon            the frozen expression, computed with ca_kernel_v2.build_fixation_field (control: = stored L)
  T                  terrace suitability alone (no traffic)
  F                  traffic alone (sum-normalised; its min-max version has identical ranks, so the same metrics)
  L_sin_F_directo    minmax(0.20 T + 0.15 Cn + 0.10 Dn): no direct F term, traffic still enters through C and D
  L_TF_nativo        minmax(0.20 T + 0.55 F): no C, D; F keeps its sum-normalised scale
  L_TF_escalado      minmax(0.20 T + 0.55 minmax(F)): exploratory, F rescaled to [0, 1]; nominal weights, not tuned
  L_FCD_sinT         minmax(0.55 F + 0.15 Cn + 0.10 Dn): every traffic term, no T
  baja_elevacion     1 - h (reference comparator of Table 10)
Metrics: AUC of registered site cells against the other valid cells; cell-level Spearman against L_canon and T;
Spearman between zone mean and site count (rho) for P2 and the stored Ward/SKATER partitions with 20 zones; Spearman
between zone rankings and those of L_canon and T; top-1 and top-3 zones; cell capture at 5/10/20/30 % with the frozen
tie rule (final_evaluate.Capture); overlap of the top-K cells (K = 168 and 10 % of Omega) with those of L_canon and T.
No p-values are reported (spatial autocorrelation). Outputs in revision/out/S3_ablation/ (P1A_*.csv);
this is a copy of audit_post_etapa6/code/audit_ablation_H4m.py with only the output folder changed."""
import hashlib
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AUD = ROOT / "revision/out/S3_ablation"; AUD.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "code"))
os.chdir(ROOT)
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")
sys.dont_write_bytecode = True

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, uniform_filter
from scipy.stats import spearmanr

import ca_kernel_v2 as ca
import final_common as fc
import final_runs as fr
import final_evaluate as fe

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
st = fc.state_for(fr.cfg())
V, S, T = st.valid_mask, st.site_mask, st.T_grid
low = np.where(V, 1.0 - st.h_grid, 0.0)
N, NS = int(V.sum()), int(S.sum())
assert (N, NS) == (12826, 84)
cap = fe.Capture(st)
P2 = fc.partition(st, 20)
Z = np.load("out/e3/partitions.npz")
PARTS = {"P2": (P2.iy.values, P2.ix.values, P2.zone.values)}
for key in ("Ward_rook_n20", "Ward_queen_n20", "SKATER_rook_n20"):
    PARTS[key.removesuffix("_n20")] = (Z["iy"], Z["ix"], Z[key])


def norm(a):
    out = np.zeros_like(a, dtype=float)
    v = a[V]
    out[V] = (v - v.min()) / (v.max() - v.min() + 1e-12)
    return out


def components(F):
    F0 = np.where(V, F, 0.0)
    return norm(gaussian_filter(F0, sigma=1.0, mode="nearest")), norm(uniform_filter(F0, size=5, mode="nearest"))


def scenarios(F):
    Cn, Dn = components(F)
    L_canon, _, _ = ca.build_fixation_field(F, T, V)
    return {
        "L_canon": L_canon,
        "T": T,
        "F": np.where(V, F, 0.0),
        "L_sin_F_directo": norm(0.20 * T + 0.15 * Cn + 0.10 * Dn),
        "L_TF_nativo": norm(0.20 * T + 0.55 * F),
        "L_TF_escalado": norm(0.20 * T + 0.55 * norm(F)),
        "L_FCD_sinT": norm(0.55 * F + 0.15 * Cn + 0.10 * Dn),
        "baja_elevacion": low,
    }, (Cn, Dn)


def zone_means(field, part):
    iy, ix, lab = part
    g = pd.DataFrame({"z": lab, "v": field[iy, ix], "s": S[iy, ix].astype(int)}).groupby("z")
    return g.v.mean(), g.s.sum()


def topk(field, k):
    v = field[cap.iy, cap.ix]
    return set(np.lexsort((cap.TIE, -v))[:k])


def metrics(name, field, ref):
    row = {"scenario": name, "auc": float(ca.auc_score(field, S, V)),
           "spearman_cell_vs_L_canon": spearmanr(field[V], ref["L_canon"][V])[0],
           "spearman_cell_vs_T": spearmanr(field[V], T[V])[0]}
    for pname, part in PARTS.items():
        m, s = zone_means(field, part)
        mL, _ = zone_means(ref["L_canon"], part); mT, _ = zone_means(T, part)
        row[f"rho_{pname}"] = spearmanr(m, s)[0]
        dom = s.idxmax()                                   # zone with most sites (the same for every scenario)
        row[f"rho_loo_{pname}"] = spearmanr(m.drop(dom), s.drop(dom))[0]
        row[f"zone_rank_vs_L_canon_{pname}"] = spearmanr(m, mL)[0]
        row[f"zone_rank_vs_T_{pname}"] = spearmanr(m, mT)[0]
        order = m.sort_values(ascending=False).index
        row[f"top1_zone_{pname}"] = int(order[0])
        row[f"top3_zones_{pname}"] = "-".join(map(str, sorted(order[:3])))
        row[f"top3_sites_{pname}"] = int(s.loc[order[:3]].sum())
    for b, c in cap.cells(field).items():
        row[f"capture_{int(b * 100)}_pct"] = 100 * c["sites"] / NS
        row[f"capture_{int(b * 100)}_tie_min_pct"] = 100 * c["tie_min"] / NS
        row[f"capture_{int(b * 100)}_tie_max_pct"] = 100 * c["tie_max"] / NS
    for k, lab in ((2 * NS, "K168"), (int(round(0.10 * N)), "top10pct")):
        mine = topk(field, k)
        row[f"overlap_{lab}_with_L_canon"] = len(mine & topk(ref["L_canon"], k)) / k
        row[f"overlap_{lab}_with_T"] = len(mine & topk(T, k)) / k
    return row


runs = [("S1", "canon_s42", 42)] + [("S2", f"canon_s{s}", s) for s in range(2000, 2030)]
rows, ctrl, scale, inputs = [], [], [], {}
for block, jid, seed in runs:
    pz = Path(f"out/final/runs/{block}/{jid}.npz"); pj = pz.with_suffix(".json")
    inputs[str(pz)] = sha(pz)
    rec = json.loads(pj.read_text())
    assert rec["npz_sha256"] == inputs[str(pz)], jid          # the field file is the one recorded by the run
    z = np.load(pz); F, Ls = z["F"], z["L"]
    sc, (Cn, Dn) = scenarios(F)
    d = float(np.abs(sc["L_canon"] - Ls)[V].max())
    ctrl.append(dict(run=f"{block}/{jid}", seed=seed, F_sum_over_Omega=float(F[V].sum()), max_abs_diff_L_canon_vs_stored=d,
                     auc_recomputed=float(ca.auc_score(Ls, S, V)), auc_in_run_json=rec["metrics"]["auc"],
                     n_attractors=rec["n_attractors"], n_pairs=rec["n_pairs"], state_hash_matches=rec["state_hash"] == fc.state_key_hash(st)))
    if d > 1e-12:
        raise SystemExit(f"STOP: L_canon does not reproduce the stored L of {jid} (max difference {d:.3e})")
    for term, arr in (("0.55*F", 0.55 * F), ("0.20*T", 0.20 * T), ("0.15*Cn", 0.15 * Cn), ("0.10*Dn", 0.10 * Dn)):
        scale.append(dict(seed=seed, term=term, mean=float(arr[V].mean()), max=float(arr[V].max()), sd=float(arr[V].std())))
    for name, field in sc.items():
        if name in ("T", "baja_elevacion") and block != "S1":
            continue                                               # seed-independent fields: computed once
        rows.append(dict(block=block, seed=seed if name not in ("T", "baja_elevacion") else "-", **metrics(name, field, sc)))
    if block == "S1":
        fz = F[V]
        scale.append(dict(seed=seed, term="F: share of cells with F = 0", mean=float((fz == 0).mean()), max=np.nan, sd=np.nan))
        top5 = np.sort(fz)[::-1][: int(round(0.05 * N))]
        scale.append(dict(seed=seed, term="F: share of traffic in the top 5 % of cells", mean=float(top5.sum() / fz.sum()), max=np.nan, sd=np.nan))

R = pd.DataFrame(rows)
C = pd.DataFrame(ctrl)
SC = pd.DataFrame(scale)

# summary: seed 42 apart; mean, SD (ddof = 1), min and max over the 30 S2 seeds
num = [c for c in R.columns if c not in ("block", "seed", "scenario") and pd.api.types.is_numeric_dtype(R[c])]
summ = []
for name, g in R.groupby("scenario", sort=False):
    s42 = g[g.block == "S1"].iloc[0]
    s2 = g[g.block == "S2"]
    for col in num:
        e = dict(scenario=name, metric=col, seed42=s42[col])
        if len(s2):
            e.update(S2_mean=s2[col].mean(), S2_sd=s2[col].std(ddof=1), S2_min=s2[col].min(), S2_max=s2[col].max(), S2_n=len(s2))
        summ.append(e)
SUM = pd.DataFrame(summ)

# cross-checks against the frozen evaluation tables
t10 = pd.read_csv("out/final/eval/tabla10_predictores.csv").set_index("predictor")
zp2 = pd.read_csv("out/final/eval/zonal_P2.csv").set_index("n")
cc = pd.read_csv("out/final/eval/captura_celda.csv")
s42 = R[R.block == "S1"].set_index("scenario")
checks = [
    ("AUC L_canon seed 42 = tabla10 L", s42.loc["L_canon", "auc"], t10.loc["L", "auc"]),
    ("rho P2 L_canon seed 42 = tabla10 L rho20", s42.loc["L_canon", "rho_P2"], t10.loc["L", "rho20"]),
    ("rho P2 L_canon seed 42 = zonal_P2 n=20", s42.loc["L_canon", "rho_P2"], zp2.loc[20, "rho"]),
    ("AUC T = tabla10 T", s42.loc["T", "auc"], t10.loc["T", "auc"]),
    ("rho P2 T = tabla10 T rho20", s42.loc["T", "rho_P2"], t10.loc["T", "rho20"]),
    ("AUC 1-h = tabla10 1 - h", s42.loc["baja_elevacion", "auc"], t10.loc["1 - h", "auc"]),
    ("rho P2 1-h = tabla10 1 - h rho20", s42.loc["baja_elevacion", "rho_P2"], t10.loc["1 - h", "rho20"]),
    ("AUC F seed 42 = tabla10 F", s42.loc["F", "auc"], t10.loc["F", "auc"]),
    ("rho P2 F seed 42 = tabla10 F rho20", s42.loc["F", "rho_P2"], t10.loc["F", "rho20"]),
]
for b in (5, 10, 20, 30):
    for lab, fld in (("L", "L_canon"), ("T", "T"), ("1 - h", "baja_elevacion")):
        sel = cc[(cc.predictor == lab) & (cc.budget_pct == b)]
        if lab == "L":
            sel = sel[sel.field.str.contains("42")]
        if len(sel) == 1:
            checks.append((f"capture {b} % {lab} = captura_celda", s42.loc[fld, f"capture_{b}_pct"], float(sel.sites_pct.iloc[0])))
s2L = R[(R.block == "S2") & (R.scenario == "L_canon")].set_index("seed")
for s in range(2000, 2030):
    for b in (5, 10, 20, 30):
        sel = cc[(cc.predictor == "L") & (cc.field == f"seed {s}") & (cc.budget_pct == b)]
        checks.append((f"capture {b} % L seed {s} = captura_celda", s2L.loc[s, f"capture_{b}_pct"], float(sel.sites_pct.iloc[0])))
sens = pd.read_csv("out/final/eval/sensibilidades_por_semilla.csv")
CK = pd.DataFrame([dict(check=n, audit=a, frozen=f, abs_diff=abs(a - f), ok=abs(a - f) <= 1e-12) for n, a, f in checks])
if not CK.ok.all() or not C.state_hash_matches.all() or not (C.auc_recomputed - C.auc_in_run_json).abs().max() <= 1e-12:
    print(CK.to_string()); raise SystemExit("STOP: a control does not reproduce the frozen values")

R.to_csv(AUD / "P1A_ABLACION.csv", index=False)
SUM.to_csv(AUD / "P1A_ABLACION_RESUMEN.csv", index=False)
C.to_csv(AUD / "P1A_CONTROLES_RUNS.csv", index=False)
CK.to_csv(AUD / "P1A_CONTROLES_TABLAS.csv", index=False)
SC.to_csv(AUD / "P1A_ESCALAS.csv", index=False)
(AUD / "P1A_ENTRADAS.sha256").write_text("".join(f"{h}  {p}\n" for p, h in inputs.items()), encoding="utf-8")
print(f"{len(runs)} runs, {len(R)} rows; controls: {len(CK)} table checks OK, max |L_canon - stored L| = {C.max_abs_diff_L_canon_vs_stored.max():.1e}")
show = SUM[SUM.metric.isin(["auc", "rho_P2", "spearman_cell_vs_T", "capture_10_pct", "overlap_K168_with_T"])]
print(show.pivot(index="scenario", columns="metric", values="seed42").round(4).to_string())
print(show.pivot(index="scenario", columns="metric", values="S2_mean").round(4).to_string())
