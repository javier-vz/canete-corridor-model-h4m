"""Reformulated Table 10: five fields (1 - h, T, E', F, L) by cell-scale and zone-scale metrics, from stored results only.

Columns: cell AUC (site cells vs background), share of site cells in the top 10 % of cells, Spearman rho between zone
mean and site count for the reference K-means partition P2 and for Ward with rook adjacency (Ward-4), 20 zones each.
Each row says which realisation it describes: deterministic (does not depend on a seed), seed 42 (canonical illustrative
run) or the 30 canonical realisations (seeds 2000-2029: mean and range). Values from different realisations are never mixed
in one cell. Sources: out/final/eval/tabla10_predictores.csv, out/final/eval/captura_celda.csv,
revision/out/S3_ablation/P1A_ABLACION.csv and revision/out/S5_counts_densities/RHO_CONTEOS_DENSIDADES.csv.
Also writes the FETE comparison (out/final/eval/tabla8_fete.csv) as a separate small table for the text.
Outputs (revision/out/table3/): TABLA10_REFORMULADA.csv, tabla10_filas.tex, FETE_TEXTO.csv."""
from pathlib import Path
import pandas as pd, numpy as np
W = Path(__file__).resolve().parents[2]; OUT = W / "revision/out/table3"; OUT.mkdir(parents=True, exist_ok=True)
t10 = pd.read_csv(W / "out/final/eval/tabla10_predictores.csv").set_index("predictor")
cap = pd.read_csv(W / "out/final/eval/captura_celda.csv"); cap = cap[cap.budget_pct == 10]
ab = pd.read_csv(W / "revision/out/S3_ablation/P1A_ABLACION.csv")
rz = pd.read_csv(W / "revision/out/S5_counts_densities/RHO_CONTEOS_DENSIDADES.csv")
def rho(part, field, real):
    return rz[(rz.partition == part) & (rz.field == field) & (rz.realisation == real)].rho_count.iloc[0]
rows = []
for name, key, capkey, fkey in [("Low elevation, 1 - h", "1 - h", "1 - h", "1 - h"), ("Topographic suitability, T", "T", "T", "T"),
                                ("Corrected affinity, E'", "E' (final affinity)", "E'", "E'")]:
    c = cap[cap.field == capkey].iloc[0]
    rows.append(dict(field=name, realisation="deterministic", auc=t10.loc[key, "auc"], capture10=c.sites_pct,
                     capture10_ties=f"{c.tie_min_pct:.1f}-{c.tie_max_pct:.1f}" if c.tie_max_pct > c.tie_min_pct else "",
                     rho_kmeans=t10.loc[key, "rho20"], rho_ward4=rho("Ward-4 (rook)", fkey, "deterministic")))
for name, sc, fkey in [("Simulated traffic, F", "F", "F"), ("Composite field, L", "L_canon", "L")]:
    s42 = ab[(ab.block == "S1") & (ab.scenario == sc)].iloc[0]; s2 = ab[(ab.block == "S2") & (ab.scenario == sc)]
    rows.append(dict(field=name, realisation="seed 42", auc=s42.auc, capture10=s42.capture_10_pct,
                     capture10_ties=f"{s42.capture_10_tie_min_pct:.1f}-{s42.capture_10_tie_max_pct:.1f}" if s42.capture_10_tie_max_pct > s42.capture_10_tie_min_pct else "",
                     rho_kmeans=s42.rho_P2, rho_ward4=s42.rho_Ward_rook))
    rows.append(dict(field=name, realisation="30 realisations: mean (range)", auc=s2.auc.mean(), auc_min=s2.auc.min(), auc_max=s2.auc.max(),
                     capture10=s2.capture_10_pct.mean(), capture10_min=s2.capture_10_pct.min(), capture10_max=s2.capture_10_pct.max(),
                     rho_kmeans=s2.rho_P2.mean(), rho_kmeans_min=s2.rho_P2.min(), rho_kmeans_max=s2.rho_P2.max(),
                     rho_ward4=s2.rho_Ward_rook.mean(), rho_ward4_min=s2.rho_Ward_rook.min(), rho_ward4_max=s2.rho_Ward_rook.max()))
T = pd.DataFrame(rows)
# controls: L seed 42 AUC and rho = Table 10; L capture seed 42 = captura_celda; F/L Ward-4 seed 42 = zonal diagnostic
assert abs(T[(T.field.str.startswith("Composite")) & (T.realisation == "seed 42")].auc.iloc[0] - t10.loc["L", "auc"]) < 1e-12
assert abs(T[(T.field.str.startswith("Simulated")) & (T.realisation == "seed 42")].rho_kmeans.iloc[0] - t10.loc["F", "rho20"]) < 1e-12
assert abs(T[(T.field.str.startswith("Composite")) & (T.realisation == "seed 42")].capture10.iloc[0] - cap[cap.field == "L seed 42"].sites_pct.iloc[0]) < 1e-9
for nm, f in (("Simulated", "F"), ("Composite", "L")):
    assert abs(T[(T.field.str.startswith(nm)) & (T.realisation == "seed 42")].rho_ward4.iloc[0] - rho("Ward-4 (rook)", f, "seed 42")) < 1e-9
T.to_csv(OUT / "TABLA10_REFORMULADA.csv", index=False)

def rng(v, lo, hi, d):
    if pd.isna(lo):
        return f"{v:.{d}f}"
    if abs(hi - lo) < 10 ** (-d) / 2:
        return f"{v:.{d}f} (constant)"
    return f"{v:.{d}f} ({lo:.{d}f}--{hi:.{d}f})"
tex = []
for r in T.itertuples():
    name = {"Low elevation, 1 - h": "Low elevation, $1-\\tilde h$", "Topographic suitability, T": "Topographic suitability, $T$",
            "Corrected affinity, E'": "Corrected affinity, $E'$", "Simulated traffic, F": "Simulated traffic, $F$",
            "Composite field, L": "Composite field, $L$"}[r.field]
    real = {"deterministic": "--", "seed 42": "seed 42", "30 realisations: mean (range)": "30 realisations"}[r.realisation]
    cap10 = rng(r.capture10, getattr(r, "capture10_min", np.nan), getattr(r, "capture10_max", np.nan), 1)
    if isinstance(r.capture10_ties, str) and r.capture10_ties:
        cap10 = f"{r.capture10:.1f}$^{{a}}$"
    tex.append(" & ".join([name if r.realisation != "30 realisations: mean (range)" else "", real,
                           rng(r.auc, getattr(r, "auc_min", np.nan), getattr(r, "auc_max", np.nan), 4), cap10,
                           rng(r.rho_kmeans, getattr(r, "rho_kmeans_min", np.nan), getattr(r, "rho_kmeans_max", np.nan), 3),
                           rng(r.rho_ward4, getattr(r, "rho_ward4_min", np.nan), getattr(r, "rho_ward4_max", np.nan), 3)]) + r" \\")
(OUT / "tabla10_filas.tex").write_text("\n".join(tex) + "\n", encoding="utf-8")

f8 = pd.read_csv(W / "out/final/eval/tabla8_fete.csv"); f8.to_csv(OUT / "FETE_TEXTO.csv", index=False)
pd.set_option("display.width", 250); print(T.round(4).to_string(index=False)); print("\n".join(tex)); print(f8.round(4).to_string(index=False))
