# Glossary of Spanish names

File names, column names, command-line options and messages of the frozen computation are partly in Spanish. They
were kept because the checksums recorded before the final runs, and inside every run, depend on the exact files.

## Prefixes of script names

| Prefix | Meaning |
|---|---|
| `final_` | final computation (frozen before the final runs) |
| `b6_`, `b9_` | phases B6 and B9 of the revision: report and diagnostic scripts that only read frozen outputs |
| `b9_etapa3_` | stage 3 of the rewrite of the article |
| `test_` | tests |

## Words in file and folder names

| Spanish | English |
|---|---|
| `grilla_con_W_hidrica` | grid with the hydrological *W* |
| `celdas_fluviales` | river cells |
| `tabla3_pesos`, `tabla4_T`, `tabla5_cuantil` | Table 3, weights; Table 4, *T*; Table 5, quantile |
| `tabla9_fete_presupuestos`, `tabla10_predictores` | Table 9, FETE budgets; Table 10, predictors |
| `tabla11_hidrologia`, `tabla12_particiones` | Table 11, hydrology; Table 12, partitions |
| `cifras_texto` | numbers of the text |
| `captura_celda`, `captura_zonas`, `resumen` | capture by cells, capture by zones, summary |
| `zonas_prefijo` | zones ranked by mean *L* with cumulative values |
| `solapamiento` | overlap |
| `descomposicion_L`, `normalizacion` | decomposition of *L*, normalisation |
| `pendiente` | slope |
| `sensibilidades_por_semilla`, `hidrologia_por_semilla` | sensitivities by seed, hydrology by seed |
| `resolucion` | resolution |
| `pares_rutas` | pairs and routes |
| `indice` | index |
| `barrido`, `resultados_barrido` | sweep, sweep results |
| `para_enviar` | to be sent |
| `verificacion_cruzada_2_trabajos` | cross-check of two jobs |
| `prueba`, `prueba_OK` | test, test passed |
| `informe` | report |
| `diagnostico_E_atractores_sitios` | diagnostic of *E* at attractors and site-labelled cells |
| `rangos_campos_figura4` | ranges of the fields of Figure 4 |
| `atractores_sensibilidades` | attractors of the sensitivity configurations |
| `conteo_atractores_8_lineas` | count of attractors near the eight lines |
| `diagnostico_L_vs_T` | diagnostic of *L* against *T* |
| `previa`, `pruebas_previas`, `controles` | preliminary, pre-tests, controls |
| `ks_procedencia` | provenance of the KS values |
| `post_cambio`, `sha_antes` | after the change, sha256 before |
| `PROCEDIMIENTOS` | procedures |
| `MANIFEST_congelado` | frozen manifest |
| `REGISTRO_CAMBIOS_POST_CONGELACION` | log of changes after the freeze |
| `fig8_captura_zonas_ref_uniforme`, `fig9_correlaciones` | Figure 8, capture by zones with the uniform-density reference; Figure 9, correlations |
| `comprobacion`, `leyenda_antes_despues` | check, legend before and after |
| `rutas_s42` (layers `rutas_canonica`, `rutas_determinista`) | routes of seed 42 (canonical, deterministic) |
| `zonas_Lmedio_P2`, `zonas_id_P2`, `zonas_centroides_P2` | mean *L*, identifier and centroid of the P2 zones |
| `atractores`, `sitios` | attractors, sites (in file names: the site-labelled cells) |
| `datos`, `referencia`, `huella` | data, reference, fingerprint |
| `LEEME` | read me |

## Words in columns and JSON keys

| Spanish | English |
|---|---|
| `altitud`, `pendiente` | elevation, slope |
| `presencia_sitio_nuevo` | label of the cell: 1 when its centre lies within about 150 m of a survey record of the reference set (site-labelled cell) |
| `distancia_sitios_nuevos`, `distancia_mincul` | distance from the cell centre to the nearest survey record of the reference set; distance to the nearest site of a layer of the national heritage registry (MINCUL), not used in the evaluation |
| `descripcion`, `valor` | description, value |
| `bloque` | block |
| `media`, `de` | mean, standard deviation |
| `segundos` | seconds |
| `trabajo`, `metrica`, `aqui`, `dif_abs` | job, metric, here (main machine), absolute difference |
| `firma`, `firma_calculo`, `firma_paquete` | signature of code and data, computation signature, package signature |
| `fecha` | date |
| `archivos` | files |
| `procedimientos_sha256` | sha256 of the procedures |

## Command-line options and messages

| Spanish | English |
|---|---|
| `--bloques` | blocks to process |
| `--verificar` | only check the stored runs (or the freeze manifest) |
| `--salida` | output folder |
| `--solo N` | only the first N pending jobs |
| `--prueba`, `--lanzar`, `--resumen` | environment test, launch the runs, summary (sweep script) |
| `válidos`, `inválidos`, `faltan` | valid, invalid, missing |
| `nada pendiente` | nothing pending |
| `firma de código/datos distinta` | different code/data signature |
