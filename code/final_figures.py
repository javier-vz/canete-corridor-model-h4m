"""JCAA-80 final recompute: Figures 8 and 9, drawn only from the evaluation CSVs (single source).
Fig 8: cumulative share of sites against cumulative area for complete zones ranked by mean L, for
       the n = 20, 30 and 40 partitions (P2 procedure). One point per complete zone, steps between
       them, nothing interpolated; the diagonal is the random expectation.
Fig 9: zonal Spearman rho (all zones) and rho_loo (without the zone with most sites) by granularity.
Outputs: out/final/fig/fig8_captura_zonas.(png|pdf), fig9_correlaciones.(png|pdf)"""
import os as _os, sys as _sys
from pathlib import Path as _Path
_ROOT = _Path(_os.environ.get("JCAA_ROOT", _Path(__file__).resolve().parents[1]))
_os.chdir(_ROOT); _sys.path.insert(0, str(_ROOT / "code"))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EV = _Path("out/final/eval"); FIG = _Path("out/final/fig")
COL = {20: "#1b4f72", 30: "#c0392b", 40: "#7d8c2b"}
MK = {20: "o", 30: "s", 40: "^"}


def fig8(prefix, zonal, path):
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    ax.plot([0, 100], [0, 100], color="0.6", lw=1, ls="--", label="Random expectation")
    for n in (20, 30, 40):
        z = prefix[prefix.n == n].sort_values("rank")
        x = np.r_[0, z.cum_area_pct.values]; y = np.r_[0, z.cum_sites_pct.values]
        ax.step(x, y, where="post", color=COL[n], lw=1.4)
        ax.plot(x[1:], y[1:], MK[n], color=COL[n], ms=4.5,
                label=f"{n} zones (ρ = {zonal.loc[n, 'rho']:.2f})")
    ax.set_xlim(0, 100); ax.set_ylim(0, 102)
    ax.set_xlabel("Cumulative area of complete zones ranked by mean L (%)")
    ax.set_ylabel("Cumulative share of recorded sites (%)")
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix("." + ext), dpi=300)
    plt.close(fig)


def fig9(zonal, path):
    ns = [20, 30, 40]; x = np.arange(3); w = 0.36
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    b1 = ax.bar(x - w / 2, zonal.loc[ns, "rho"], w, color="#1b4f72", label="All zones")
    b2 = ax.bar(x + w / 2, zonal.loc[ns, "rho_loo"], w, color="#85a9c8", label="Without the zone with most sites")
    for bars, col in ((b1, "p"), (b2, "p_loo")):
        for bar, n in zip(bars, ns):
            ax.annotate(f"{bar.get_height():.2f}\np = {zonal.loc[n, col]:.3f}", (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                        ha="center", va="bottom", fontsize=7.5, xytext=(0, 2), textcoords="offset points")
    ax.set_xticks(x, [f"{n} zones" for n in ns]); ax.set_ylim(0, 1)
    ax.set_ylabel("Spearman ρ (mean L vs site count)")
    ax.axhline(0, color="0.3", lw=0.8)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(path.with_suffix("." + ext), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    zonal = pd.read_csv(EV / "zonal_P2.csv").set_index("n")
    prefix = pd.read_csv(EV / "zonas_prefijo.csv")
    fig8(prefix, zonal, FIG / "fig8_captura_zonas")
    fig9(zonal, FIG / "fig9_correlaciones")
    print("figures written to", FIG)
