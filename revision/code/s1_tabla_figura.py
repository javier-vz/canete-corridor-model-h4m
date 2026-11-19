"""S1 (ablation of L on the stored fields): table S1 (CSV and LaTeX rows) and figure S1 with untruncated axes.
Reads 02_S1_ABLACION/P1A_ABLACION.csv and P1A_ABLACION_RESUMEN.csv (written by s1_ablacion_campos.py)."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

D = Path(__file__).resolve().parents[2] / "revision/out/S3_ablation"
A = pd.read_csv(D / "P1A_ABLACION.csv")
S = pd.read_csv(D / "P1A_ABLACION_RESUMEN.csv")
V = [("L_canon", "L as published", "F, T, C, D", "F sum-normalised (sum 1); T, C, D min-max; final min-max"),
     ("L_sin_F_directo", "L without the direct F term", "T, C, D", "as published, F term removed; final min-max"),
     ("L_TF_nativo", "T and F, without C and D", "T, F", "F sum-normalised; final min-max"),
     ("L_TF_escalado", "T and F rescaled (exploratory)", "T, F", "F min-max before weighting; final min-max"),
     ("L_FCD_sinT", "Traffic terms without T", "F, C, D", "as published, T term removed; final min-max"),
     ("F", "F alone", "F", "sum-normalised (ranks identical to min-max)"),
     ("T", "T alone (no traffic)", "T", "min-max"),
     ("baja_elevacion", "Low elevation, 1 - h", "h", "min-max")]
rows = []
for key, lab, terms, scale in V:
    g = A[A.scenario == key]; s42 = g[g.block == "S1"].iloc[0]; s2 = g[g.block == "S2"]
    rows.append({"variant": lab, "scenario_id": key, "terms": terms, "scale": scale,
                 "auc_seed42": s42.auc, "auc_30_mean": s2.auc.mean() if len(s2) else None, "auc_30_sd": s2.auc.std(ddof=1) if len(s2) else None,
                 "auc_30_min": s2.auc.min() if len(s2) else None, "auc_30_max": s2.auc.max() if len(s2) else None,
                 "rho20_P2_seed42": s42.rho_P2, "rho20_Ward_rook_seed42": s42.rho_Ward_rook, "rho20_Ward_queen_seed42": s42.rho_Ward_queen,
                 "rho20_SKATER_seed42": s42.rho_SKATER_rook, "top168_shared_with_L_pct": 100 * s42.overlap_K168_with_L_canon,
                 "top168_shared_with_T_pct": 100 * s42.overlap_K168_with_T, "capture_top10_pct_seed42": s42.capture_10_pct,
                 "spearman_cells_with_L": s42.spearman_cell_vs_L_canon, "spearman_cells_with_T": s42.spearman_cell_vs_T,
                 "site_cells": 84, "background_cells": 12742})
T1 = pd.DataFrame(rows)
T1.to_csv(D / "S1_ABLACION_TABLA.csv", index=False)

f = lambda x, d=4: "--" if pd.isna(x) else f"{x:.{d}f}"
tex = []
for r in T1.itertuples():
    tex.append(" & ".join([r.variant.replace("1 - h", "$1-\\tilde h$"), r.terms, f(r.auc_seed42),
                           (f"{r.auc_30_mean:.4f} $\\pm$ {r.auc_30_sd:.4f}" if pd.notna(r.auc_30_mean) else "--"),
                           f(r.rho20_P2_seed42, 3), f"{r.top168_shared_with_L_pct:.0f}", f"{r.top168_shared_with_T_pct:.0f}",
                           f"{r.capture_top10_pct_seed42:.1f}", f(r.spearman_cells_with_L, 3), f(r.spearman_cells_with_T, 3)]) + r" \\")
(D / "S1_tabla_filas.tex").write_text("\n".join(tex) + "\n", encoding="utf-8")

# figure: untruncated axes (AUC from 0.5, the no-discrimination value; percentages from 0 to 100)
BLUE, ORANGE, TXT, MUTED, GRID, SURF = "#2a78d6", "#eb6834", "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"
P = [("auc", "AUC (site cells vs other cells)", (0.5, 1.0)), ("capture_10_pct", "Site cells in the top 10% of cells (%)", (0, 100)),
     ("overlap_K168_with_L_canon", "Top 168 cells shared with L (%)", (0, 100)), ("overlap_K168_with_T", "Top 168 cells shared with T (%)", (0, 100))]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": MUTED, "text.color": TXT, "xtick.color": TXT, "ytick.color": TXT})
fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.9), sharey=True, facecolor=SURF)
y = {k: i for i, (k, *_ ) in enumerate(V[::-1])}
for ax, (col, title, lim) in zip(axes, P):
    ax.set_facecolor(SURF)
    mult = 100 if col.startswith("overlap") else 1
    for key, *_ in V:
        g = A[A.scenario == key]; s42 = g[g.block == "S1"][col].iloc[0] * mult; s2 = g[g.block == "S2"][col] * mult
        if len(s2):
            ax.plot([s2.min(), s2.max()], [y[key], y[key]], color=BLUE, lw=2, solid_capstyle="round", zorder=2)
            ax.plot(s2.mean(), y[key], "o", ms=6, color=BLUE, mec=SURF, mew=1.5, zorder=3)
        ax.plot(s42, y[key] + (0.22 if len(s2) else 0), "D", ms=5, color=ORANGE, mec=SURF, mew=1.2, zorder=4)
    ax.set_xlim(*lim); ax.set_title(title, fontsize=8.5, loc="left")
    ax.grid(axis="x", color=GRID, lw=0.6); ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
axes[0].set_yticks(list(y.values()), [lab for _, lab, *_ in V[::-1]])
h = [plt.Line2D([], [], marker="D", ls="none", color=ORANGE, mec=SURF, ms=5, label="Seed 42"),
     plt.Line2D([], [], marker="o", color=BLUE, mec=SURF, ms=6, lw=2, label="30 canonical realisations (seeds 2000–2029): mean and range")]
fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.6, -0.01))
fig.tight_layout(rect=(0, 0.06, 1, 1))
fig.savefig(D / "S1_ABLACION_FIGURA.pdf"); fig.savefig(D / "S1_ABLACION_FIGURA.png", dpi=200)
print(T1[["variant", "auc_seed42", "auc_30_mean", "auc_30_sd", "rho20_P2_seed42", "top168_shared_with_L_pct", "top168_shared_with_T_pct", "capture_top10_pct_seed42"]].round(4).to_string(index=False))
