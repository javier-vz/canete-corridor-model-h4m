# The model as implemented

This page states every rule and nominal parameter of the final configuration as they are coded in
[`code/ca_kernel_v2.py`](../code/ca_kernel_v2.py) (the kernel of v1.0.1 with an optional attractor-exclusion mask; not
modified for the final configuration),
[`code/final_h4m_state.py`](../code/final_h4m_state.py) (hydrological rules) and
[`code/final_runs.py`](../code/final_runs.py) (configurations and seeds). Symbols follow the article. The internal name
of the final configuration is **H4m**.

## Grid and inputs

- Domain Ω: the 12,826 valid cells of a 93 × 144 lattice of 250 m cells (EPSG:32718). In the kernel lattice row 0 is the
  southern row.
- $h$, $p$ and $W$ are the columns `h_norm`, `p_norm` and `W_corridor` of
  `data/grilla_con_W_hidrica.csv` (restricted; the public copy without the site columns is
  [`data/grilla_con_W_hidrica_publica.csv`](../data/grilla_con_W_hidrica_publica.csv); [data/README.md](../data/README.md)).
- Min–max normalisation over Ω: $\tilde x = (x - \min x) / (\max x - \min x + 10^{-12})$.

## Terrace suitability *T*

`build_terrace_field`, with $r = 0.65$, $s = 0.20$ and a window of $w = 5$ cells:

$$T_0 = r\,(1-p) + (1-r)\,U_w(1-p), \qquad T_1 = T_0\,\big((1-s) + s\,(1-h)\big), \qquad T = \widetilde{T_1},$$

where $U_w$ is the mean over a $w \times w$ window (edge cells padded with the nearest value).

## Affinity *E* and cost *c*

`compute_environment_and_cost`, with the nominal weights $a = (0.10, 0.20, 0.50, 0.20)$ and
$b = (0.18, 0.36, 0.26, 0.20)$:

$$E = 0.10\,(1-h) + 0.20\,(1-p) + 0.50\,W + 0.20\,T, \qquad c = 0.18\,h + 0.36\,p + 0.26\,(1-W) + 0.20\,(1-T).$$

## Hydrological rules of the final configuration

Let $d$ be `dist_river_m` (distance from the cell centre to the nearest HydroRIVERS line) and $R$ the 1,020 river
cells, whose 250 m footprint intersects one of the 81 HydroRIVERS lines of the model window
([`data/celdas_fluviales_1020.csv`](../data/celdas_fluviales_1020.csv)).

- Affinity: $E' = E\,\big(1 - 0.25\,e^{-d/250\,\mathrm{m}}\big)$. $E'$ replaces $E$ everywhere in the model.
- Attractor exclusion: cells in $R$ and cells with $d \le 125$ m (the 829 cells with $d \le 125$ m are all in $R$).
  Excluded cells remain part of the routing graph.
- Cost: $c' = c + 0.10$ on the cells of $R$.

## Attractors

`extract_local_maxima_adaptive`, quantile $q = 0.92$, blocking radius 5 cells, 4 bands, at least 4 attractors per band:

1. Candidates are cells outside the exclusion set with $E' \ge \tau = Q_q(E' \mid \Omega)$ that are the maximum of $E'$
   in their 11 × 11 window.
2. Candidates are accepted in decreasing order of $E'$; each accepted cell blocks the 11 × 11 square around it.
3. The lattice is split into 4 bands of equal width along its columns. A band with fewer than 4 attractors receives
   further candidates above its own quantile $\min(0.98, q + 0.03)$ of $E'$, with the same blocking rule.

Result: 25 attractors.

## Pairs

`build_local_pairs`: each attractor is paired with its 6 nearest attractors (Euclidean distance on the lattice). The
union of these unordered pairs gives 90 pairs.

## Routing graph

`build_graph`: 8-neighbour graph on Ω. The base weight of an edge is
$w_{uv} = \tfrac{1}{2}(c'_u + c'_v)\,\ell_{uv}$, with $\ell = 1$ for orthogonal and $\sqrt 2$ for diagonal moves
(lattice units).

## One realisation

`run_one(state, eps, kappa, seed)`, with $\beta = 1.5$, $k = 4$ paths, congestion cap 2 and a single
`numpy.random.default_rng(seed)` stream:

1. **Pair order.** For each pair, $C_{ab}$ is the least-cost distance on the base weights and
   $w_{ab} = e^{2E'_a}\,e^{2E'_b}\,e^{-3C_{ab}}$. Pairs are processed in decreasing $w_{ab}$.
2. **Noise.** Before routing a pair, every edge weight is set to $w_{uv}\,g_{uv}\,\max(0.1,\,1 + \varepsilon\,\xi_{uv})$,
   with $\xi_{uv} \sim N(0,1)$ drawn anew for every edge and pair and $g_{uv}$ the congestion multiplier (initially 1).
   With $\varepsilon = 0$ the factor is 1.
3. **Paths.** The 4 shortest simple paths under these weights are found (`networkx.shortest_simple_paths`). With
   $C_1 \le \dots \le C_4$ their costs, path $j$ has probability $\pi_j \propto e^{-\beta (C_j - C_1)}$, and every cell of
   path $j$ receives $e^{2E'_a}\,e^{2E'_b}\,e^{-3C_1}\,\pi_j$ of traffic.
4. **Congestion.** Each edge used by $n_{uv}$ of the pair's paths is updated to $g_{uv} \leftarrow \min(g_{uv} +
   \kappa\,n_{uv},\,2)$.
5. **Traffic.** $F$ is the accumulated traffic divided by its total, so that $\sum_\Omega F = 1$.
6. **Composite field (historically termed “fixation field”).**
   $L = \widetilde{\big(0.55\,F + 0.20\,T + 0.15\,\tilde C + 0.10\,\tilde D\big)}$, where $C$ is a Gaussian filter of $F$
   ($\sigma = 1$ cell) and $D$ the 5 × 5 mean of $F$. $F$ enters $L$ with $\sum F = 1$, not min–max normalised, so its
   direct term is small; `out/final/eval/descomposicion_L.csv` gives the contribution of each term.
7. **Run metrics.** AUC is the Mann–Whitney statistic $U / (n_1 n_0)$ of $L$ in the 84 site-labelled cells against the other
   12,742 valid cells (background cells) (ties take average ranks). $m_1$, $m_2$ and $m_3$ describe the $K = 2 \times 84 = 168$ cells with
   the highest $L$ (8-connected components): share in the largest component, share of $F$ in it, and $e^{-S}$ with $S$
   the Shannon entropy of $F$.

The canonical setting is $(\varepsilon, \kappa) = (0.15, 0.03)$; the deterministic baseline is $(0, 0)$.

## Evaluation

The evaluation procedures were written in [`out/final/PROCEDIMIENTOS.md`](../out/final/PROCEDIMIENTOS.md) (in Spanish)
on 5 November 2026, before any final value was computed, and are implemented by
[`code/final_common.py`](../code/final_common.py) and [`code/final_evaluate.py`](../code/final_evaluate.py). In short:

- **Partitions.** Reference K-means partition with 20 zones (P2): features (column, row, $h$, $p$) standardised, cells in
  row-major order, `KMeans(random_state=42, n_init=10)`; the same with 30 and 40 zones. Ward (4 and 8 neighbours),
  SKATER and 100 K-means initialisations are read from partitions stored in an earlier phase
  ([out/README.md](../out/README.md)).
- **Zonal statistics.** $\rho_n$: Spearman correlation between the mean of $L$ in each zone and its number of site-labelled cells
  (two-sided p-value of `scipy.stats.spearmanr`); $\rho_{\mathrm{loo}}$: the same without the zone with most site-labelled cells;
  cumulative capture over complete zones ranked by mean $L$, without interpolation.
- **Cell capture.** Share of the 84 site-labelled cells among the top $\mathrm{round}(b \times 12{,}826)$ cells for
  $b = 5, 10, 20, 30\,\%$; ties broken by the fixed order `numpy.random.default_rng(0).permutation` of the valid cells.
- **Kolmogorov–Smirnov.** `scipy.stats.ks_2samp(..., alternative="two-sided", method="exact")` of $E'$ and $L$ (seed 42)
  at the site-labelled cells against the 12,742 background cells, for all 84 site-labelled cells and for the 20 outside
  the zone with most site-labelled cells. p-values do not account for spatial dependence.
- **Summaries.** Mean and sample standard deviation (ddof = 1) over the seeds of each block.
- **Point-process models.** [`code/final_ppm.R`](../code/final_ppm.R) fits multi-scale models of the centres of the 84 site-labelled cells with
  $F$, $L$ and $E'$ of seed 42 as covariates, smoothed at radii of 250–5,000 m (the column `E` of its input holds $E'$).

## Configurations and seeds

| Block | Content | Seeds | Runs |
|---|---|---|---|
| S1 | canonical realisation and deterministic baseline, with fields and paths | 42 | 2 |
| S2 | canonical realisations | 2000–2029 | 30 |
| S3 | affinity and cost weights: nominal and 16 one-at-a-time perturbations of ±20 % | 3000–3002 | 51 |
| S4 | parameters $(r, s, w)$ of $T$: nominal and 6 variants | 5000–5002 | 21 |
| S5 | attractor quantile $q$ = 0.90, 0.92, 0.95, 0.975 | 4000–4002 | 12 |
| S6 | $(\varepsilon, \kappa)$ sweep: $\varepsilon$ = 0–0.40 in steps of 0.04, $\kappa$ = 0–0.10 in steps of 0.01 | 1000–1009 | 1,210 |
| S7 | FETE baseline on $c'$ with 50, 90, 200, 500 and 1,000 random pairs | 2000–2029 | 150 |
| S8 | FETE baseline on $1 - E'$ with 90 random pairs | 2000–2029 | 30 |
| S9 | slope transformations (clipping, $\log(1+p)$, FAO classes, cell average) | 42 | 4 |
| CTRL | three configurations of the previous model version (v2: 125 m river buffer, without the rules of $E'$ and $c'$), recomputed to check that its published values are reproduced | 3000–3002, 4000–4002, 5000–5002 | 9 |

The hydrographic sensitivity of Table 11 (final, eight unverified reaches removed, main channel only, no cost penalty)
uses runs with seeds 42–44 stored in `out/b4/runs/` and `out/b5/runs/` ([out/README.md](../out/README.md)).
