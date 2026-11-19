"""Routing budget and run time of the corridor model and the two FETE baselines (Table S14). No run is repeated.
(FETE part of s4_fete_s6_mincul.py of the authors' working copy.)

S4 reads the run JSON files of S2 (canonical H4m, 30 seeds), S7 (FETE on c', 90 pairs) and S8 (FETE on 1 - E', 90 pairs) and
the path record of S1. Run time is the wall-clock time of each job, measured by final_runs.py with time.time() in its
worker process (main machine, two worker processes in parallel); it is indicative, not a benchmark.
The MINCUL summary only uses the analysis grid and the stored partitions: counts, never coordinates."""
import glob, json, os, pickle, sys
from pathlib import Path
W = Path(__file__).resolve().parents[2]; os.chdir(W); sys.path.insert(0, str(W / "code")); sys.dont_write_bytecode = True
import numpy as np, pandas as pd
import final_common as fc, final_runs as fr

E = W / "revision/out"
rows = []
for label, pat, endpoints, routes in [
        ("Corridor model H4m (canonical, eps 0.15, kappa 0.03)", "out/final/runs/S2/canon_s*.json", "90 pairs of the 25 attractors (6 nearest neighbours)", 4),
        ("FETE on c' (deterministic least-cost path)", "out/final/runs/S7/fete_c_n90_s*.json", "90 random pairs of valid cells per seed", 1),
        ("FETE on 1 - E' (deterministic least-cost path)", "out/final/runs/S8/fete_1mE_n90_s*.json", "90 random pairs of valid cells per seed", 1)]:
    J = [json.load(open(p)) for p in sorted(glob.glob(pat))]
    rt = np.array([j["runtime_s"] for j in J]); auc = np.array([j["metrics"]["auc"] for j in J]); m2 = np.array([j["metrics"]["m2"] for j in J])
    seeds = sorted(j["spec"]["seed"] for j in J)
    rows.append({"method": label, "runs": len(J), "seeds": f"{seeds[0]}-{seeds[-1]}", "endpoints": endpoints, "pairs_per_run": 90,
                 "routes_per_pair": routes, "routes_per_run": 90 * routes, "noise_and_congestion": "yes" if routes == 4 else "no",
                 "auc_mean": auc.mean(), "auc_sd": auc.std(ddof=1), "m2_mean": m2.mean(), "m2_sd": m2.std(ddof=1),
                 "runtime_s_median": float(np.median(rt)), "runtime_s_p25": float(np.percentile(rt, 25)), "runtime_s_p75": float(np.percentile(rt, 75)),
                 "runtime_status": "indicative: wall-clock per job, same machine, 2 parallel workers; not a benchmark"})
paths = pickle.load(open("out/final/runs/S1/canon_s42_paths.pkl", "rb"))
assert len(paths) == 90 and all(len(p["paths"]) == 4 for p in paths)          # 360 routes in the canonical run of seed 42
S4 = pd.DataFrame(rows)
(E / "S6_fete_budget").mkdir(parents=True, exist_ok=True)
S4.to_csv(E / "S6_fete_budget/S4_FETE_PRESUPUESTO.csv", index=False)
t8 = pd.read_csv("out/final/eval/tabla8_fete.csv")
chk = [(S4.auc_mean[0], t8[t8.method.str.contains("corridor")].auc_mean.iloc[0]), (S4.auc_mean[1], t8[t8.method == "FETE, cost c'"].auc_mean.iloc[0]),
       (S4.auc_mean[2], t8[t8.method == "FETE, cost 1 - E'"].auc_mean.iloc[0])]
assert all(abs(a - b) < 1e-12 for a, b in chk), "S4 does not reproduce Table 8"
print(S4.to_string(index=False))
