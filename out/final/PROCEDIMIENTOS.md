# JCAA-80 — Procedimientos congelados del recálculo final (H4m)

Fecha: 5 de noviembre de 2026, antes de calcular cualquier valor definitivo de H4m.
Este archivo entra en el manifiesto congelado. Lo implementan `code/final_common.py` y `code/final_evaluate.py`.
**Ninguna regla se elige ni se cambia según el resultado que produce.**

## 1. Campos y realizaciones

| Uso | Campo y semillas |
|---|---|
| Mapas, Tabla 10, KS, PPM, §4.2–4.4 | realización canónica S1, semilla 42 (ε = 0.15, κ = 0.03) |
| Fila del modelo en la Tabla 8; captura por semilla; mapa medio | conjunto S2, semillas 2000–2029 |
| Tablas 3, 4 y 5 | S3 (3000–3002), S4 (5000–5002), S5 (4000–4002) |
| Rango del barrido ε/κ | S6 (1000–1009 en 121 puntos) |
| FETE | S7 (c′; 50, 90, 200, 500 y 1000 pares) y S8 (1 − E′; 90 pares), semillas 2000–2029 |

**E′ es el campo de afinidad del modelo final.** Donde el artículo dice E como campo del modelo (Figuras 4 y 5, KS, PPM, Tabla 10, FETE 1 − E), se usa E′. En la Tabla 10, E (sin H2m) puede ir como fila adicional.

## 2. Particiones

- **P2 (n = 20):** `canonical_partition(state, n_zones=20, random_state=42, n_init=10)`, con variables (columna, fila, h̃, p̃) estandarizadas y celdas en orden fila-columna.
- **n = 30 y 40:** el mismo procedimiento, con `n_zones = 30` o `40`.
- **Ward con 4 y 8 vecinos, SKATER y las 100 inicializaciones de K-means:** las particiones guardadas en la Fase B (`out/e3/partitions.npz`, `out/kmeans_labels_100.npy`). No se recalculan.

## 3. Métricas zonales

- **ρ_n:** Spearman entre L̄_z (media de L en la zona) y s_z (número de sitios). El valor p es el de `scipy.stats.spearmanr`: bilateral, aproximación t.
- **ρ con L_max:** lo mismo, con el máximo de L en cada zona.
- **ρ_loo:** se quita la zona con más sitios y se recalcula ρ con su valor p.
- **Prefijo de zonas completas:**
  - las zonas se ordenan por L̄ descendente;
  - en cada zona completa se informan el área acumulada (celdas / 12,826), los sitios acumulados (/ 84) y el lift (fracción de sitios / fracción de área);
  - **nada se interpola.**
  - Con presupuesto, se suman zonas en ese orden mientras el área acumulada no pase el presupuesto. Se informa el área realmente consumida.

## 4. Captura por celda

- **Presupuesto exacto:** n_b = round(b × 12,826) celdas para b = 5, 10, 20 y 30 %.
- **Orden y empates:**
  - las celdas se ordenan por valor descendente;
  - los empates se rompen con un orden fijo, `numpy.random.default_rng(0).permutation`, sobre las celdas válidas en orden fila-columna;
  - se informan además tie_min y tie_max.
- **Resumen:**
  - por semilla (S2): media, DE, mínimo y máximo;
  - mapa medio: la media celda a celda de los 30 L, evaluada una sola vez.
- **Comparadores:** 1 − h̃, T y E′.

## 5. Estadísticos de celda

- **AUC:** `auc_score` del kernel (U de Mann–Whitney / (n_pos · n_neg)) sobre Ω. Positivos: las 84 celdas con sitio.
- **|C₁|, m₁, m₂, m₃:** los que devuelve `run_one` del kernel. Usa las 168 celdas top de L (K = 2 × 84) y componentes 8-conexas.
- **Solapamiento con búfer Ω^K_b** (ecuación del artículo):
  - T_K = {celdas con L ≥ L_(K)}, donde L_(K) es el K-ésimo mayor valor;
  - B_b = dilatación de Chebyshev de las celdas con sitio, con radio b celdas;
  - Ω^K_b = |T_K ∩ B_b| / |T_K|.
  - Se informa para b ∈ {0, 1, 2} y K ∈ {0.5, 1, 2, 5, 10, 15} × 84.
- **Probabilidad del camino más barato (§3.3):**
  - para cada par de la semilla 42: P = exp(−β (c_π − c*)) normalizada sobre los 4 caminos, con β = 1.5 y los costos perturbados que registra el observador;
  - se informan el mínimo y el máximo entre pares.

## 6. Kolmogorov–Smirnov (Tabla 6) — procedimiento fijado

**Prueba.**
- Dos muestras, bilateral: `scipy.stats.ks_2samp(x_sitios, x_fondo, alternative="two-sided", method="exact")`.
- Se registra el método que SciPy usa de verdad. Si no puede calcular el exacto, cambia al asintótico con un aviso, y eso se anota.

**Grupos.**
- **Fondo:** las 12,742 celdas válidas sin sitio, el mismo en todas las comparaciones.
- **Conjunto residual:** los sitios fuera de la zona de P2 (n = 20) con más sitios. Esa zona tiene 64 sitios, así que quedan 20 residuales. Los sitios excluidos no pasan al fondo.

**Campos:** E′ y L de la semilla 42.

**Justificación del método exacto.**
- Las muestras de sitios son pequeñas (84 y 20).
- La distribución asintótica es una aproximación para n efectivo ≈ n₁.
- El modo `auto` de SciPy usaría la asintótica solo porque el fondo supera 10,000 celdas, que es un umbral de cálculo y no un criterio estadístico.

**Procedencia.**
- La implementación histórica se reproduce solo para documentar de dónde salieron los valores publicados (`out/final/previa/ks_procedencia.csv`).
- No interviene en la elección del método.

**Alcance.** Los valores p no consideran la dependencia espacial. Se presentan como diagnósticos.

## 7. Tablas de resumen

- **Tablas 3–5 y 8–9:** media ± DE muestral (ddof = 1) sobre las semillas de cada bloque.
- **Barrido:** media por (ε, κ) sobre 10 semillas. El texto informa el mínimo y el máximo de esas medias.
- **Única fuente numérica:**
  - toda cifra del texto y de las tablas sale de un CSV de `out/final/eval/`;
  - el redondeo se hace solo al escribir el LaTeX.

## 8. PPM

- `code/final_ppm.R` es `multi_scale_ppm_v2.R` con solo los nombres de archivo cambiados.
- Covariables: F, L y E′ de la semilla 42; los mismos radios.
- Antes de usarlo, el entorno de R debe reproducir la Tabla 7 publicada desde `fields_canonical_v2_for_ppm.csv`.
