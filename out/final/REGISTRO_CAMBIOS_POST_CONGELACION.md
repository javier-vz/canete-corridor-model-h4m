# JCAA-80 — Registro de cambios posteriores a la congelación

**Manifiesto original:** `out/final/MANIFEST_congelado.sha256` (54 archivos, congelado el 06/11 a las 09:41, sha256 2ffc3e90…).
**Firma de cálculo de las corridas:** 9bd72e1a…. Cubre `final_common.py`, `final_runs.py`, `final_h4m_state.py`, `ca_kernel_v2.py`, `canonical_partition.py`, la grilla, las celdas fluviales y las variantes de pendiente.

## Verificación antes de cualquier cambio

06/11, 10:57, al terminar S1–S9:
- `final_runs.py --verificar`: 309 de 309 trabajos válidos (S1–S9 y CTRL).
- `final_freeze.py --verificar`: **manifest OK**, los 54 archivos sin cambios.

## Cambio 1 — portabilidad de las pruebas (punto 5 de la evaluación B7)

- **Cuándo:** 06/11, 11:10, después de terminar S1–S9 y sin tocar ningún archivo de la firma de cálculo.
- **Archivo:** `code/test_final_pipeline.py`. Antes 54a4cd5c…, después bffff320….
  - El directorio temporal de las pruebas del ejecutor pasa de una ruta absoluta del entorno a `tempfile.mkdtemp(prefix="jcaa80_test_runs_")`.
  - La entrada PPM v2 publicada se lee desde `data/ref/fields_canonical_v2_for_ppm.csv` (archivo nuevo). Es una copia sin cambios del original de `pre_paper_v3`, con el mismo sha256: 7d03e51f….
  - La exportación de prueba se escribe en un `TemporaryDirectory`.
- **Después del cambio:**
  - se repitieron las pruebas: 34 OK y 3 notas, idénticas a las anteriores salvo los tiempos (`pruebas_previas_post_cambio_2026-11-06.csv`);
  - se repitieron los controles v2: 12 OK;
  - se restauró en `out/final/previa/pruebas_previas.csv` el registro congelado original.
- **Estado del manifiesto original:** `final_freeze.py --verificar` informa como cambiado solo `code/test_final_pipeline.py`. Es el cambio documentado aquí.

## Archivos nuevos posteriores a la congelación (capa de informe)

Solo leen salidas congeladas y no cambian ninguna:
- `code/b8_diferencias.py`: matriz de diferencias con el manuscrito v8_11.
- `code/b8_concordancia.py`: concordancia entre tablas, figuras y estadísticas.
- `code/b8_informe.py`: informe de recálculo.
- `code/b8_manifest.py`: manifiesto del paquete reproducible.
- `data/ref/fields_canonical_v2_for_ppm.csv`: referencia del cambio 1.

## Sin cambios

- Kernel, estado H4m, ejecutor, evaluación, figuras, exportaciones y PPM.
- El paquete del Ryzen (`JCAA80_barrido_H4m_Ryzen.zip`, 341c0840…).

## Integración de S6 (06/11, 17:55)

**Recepción y validación del zip del Ryzen.**
- Archivo recibido: `para_enviar_barrido_H4m.zip` (sha256 bc1ee832…).
- Su `MANIFEST.sha256` es idéntico al del paquete enviado (0731fd33…) y su `prueba_OK.json` es idéntico al recibido a las 09:14 (firma f7a44f52…).
- Validación con las reglas del propio script: 1,210 de 1,210 válidas, combinaciones únicas e iguales a las esperadas, sin archivos sobrantes.

**Verificación cruzada.**
- Se repitieron aquí 2 trabajos elegidos al azar con `default_rng(20261106)`: `e0.000_k0.100_s1004` y `e0.080_k0.050_s1006`.
- Se usó el mismo paquete, con prueba local de la misma firma.
- Diferencia máxima con el Ryzen: 4.3e−18 (`out/final/S6/verificacion_cruzada_2_trabajos.csv`).

**Re-ejecución de `final_evaluate.py` (sin cambios de código).**
- Agregó `barrido_S6.csv` y las cifras X7.
- Los 24 CSV ya existentes se reprodujeron **byte a byte**. Solo cambiaron `cifras_texto.csv` (filas X7 nuevas) e `indice.json`.
- El resumen calculado aquí es idéntico al `resumen_barrido_H4m.csv` del Ryzen.

**Figuras 8 y 9.** Se regeneraron desde los mismos CSV y son idénticas píxel a píxel a las entregadas. Los bytes del PNG difieren solo por la compresión.

**Archivos nuevos:** `out/final/S6/` (resultados, prueba, manifiesto, resumen del Ryzen, zip recibido y verificación cruzada).

## Fase B9: capa de informe para la reescritura (06/11, 18:50)

**Qué se agregó.** Solo archivos nuevos. Leen estados y salidas congeladas; no ejecutan corridas ni modifican archivos existentes.
- `code/b9_diag_E.py` → `out/final/informe/diagnostico_E_atractores_sitios.csv` y `rangos_campos_figura4.csv`.
- `code/b9_diag_atractores.py` → `atractores_sensibilidades.csv` y `conteo_atractores_8_lineas.csv`. Los 28 estados de las Tablas 3–5 se reconstruyeron con el código congelado y coinciden con el `state_hash` registrado en su corrida.
- `code/b9_matriz_cambios.py` → `MATRIZ_CAMBIOS_MANUSCRITO.csv`, `MATRIZ_CAMBIOS_MANUSCRITO.md` y `cifras_fuera_de_matriz.csv`.
- `code/b9_nota_conteo.py` → `NOTA_CORRECCION_CONTEO_ATRACTORES.md`.
- Manifiesto de estos archivos: `out/final/informe/MANIFEST_B9.sha256`.

**Verificación antes de agregar.**
- `final_freeze.py --verificar`: solo cambia `code/test_final_pipeline.py` (cambio 1).
- `MANIFEST_paquete_reproducible.sha256`: 1,766 de 1,766 archivos sin cambios.

**Este registro.** Se amplió con esta sección. El manifiesto del paquete B8 conserva el sha256 de la versión anterior (`ae8b5002…`); la versión nueva figura en `MANIFEST_B9.sha256`.

**Regeneración sobre la versión nueva del manuscrito (06/11, 19:30).** Los autores enviaron `main_v8 (13).tex` («13» es el contador de descarga), guardado como `data/main_v8_13.tex` (sha256 `9eca8e68…`). `b9_matriz_cambios.py` se actualizó para usar esa versión: conserva su nombre, agrega la alineación de líneas con v8_11 y adapta las entradas de los pasajes reescritos. La matriz, la lista de control y el manifiesto basados en v8_11, junto con la versión anterior del script, quedan en `out/final/informe/historial/` con el sufijo `_base_v8_11`.

## Etapa 3 de la reescritura: Figura 8 con la diagonal como referencia de densidad uniforme (06/11, 23:21)

**Qué se agregó.** Solo archivos nuevos, en `out/final/fig_etapa3/`, generados por `code/b9_etapa3_fig8.py`:
- `fig8_captura_zonas_ref_uniforme.png` y `.pdf`;
- `fig8_comprobacion.json`;
- `fig8_leyenda_antes_despues.png`.

El script usa la función congelada `final_figures.fig8()` y los mismos CSV. Solo cambia el texto de la diagonal en la leyenda.

**Sin cambios:** `code/final_figures.py` y `out/final/fig/` (sha256 comprobados).

**Verificación.** `final_freeze.py --verificar`: solo cambia `code/test_final_pipeline.py`, el cambio 1 ya documentado.

**Paquete para Dina.** Está fuera de `out/final`, en `out/b9/dina/`, y lee las capas de `out/final/qgis` sin modificarlas.
