# S2 — Los ocho tramos HydroRIVERS sin verificar: criterio y fuentes

No se recalculó H4m. Todo sale de archivos guardados.

## Fuentes

| Qué | Fuente |
|---|---|
| Atributos oficiales (`HYRIV_ID`, orden de Strahler, longitud) | shapefile oficial HydroRIVERS v1.0 de Sudamérica (Lehner y Grill 2013), con su .dbf |
| Identificación de los ocho tramos | los números de registro (FID) que usa `final_evaluate.py`; cada línea del modelo coincide con su registro oficial (distancia de Hausdorff 0 m) |
| Revisión de imagen | formulario de Dina, hoja «8 prioritarias», devuelto el 05/11/2026; imagen satelital de Google Earth; **la fecha de la imagen no se registró** |
| Tráfico atribuido a cada tramo | campo F canónico de la semilla 42, sumado sobre las celdas de 250 m que el tramo cruza |
| Atractores con y sin los ocho tramos | `out/final/informe/conteo_atractores_8_lineas.csv` |
| Efecto en las métricas | Tabla 11 congelada (`out/final/eval/tabla11_hidrologia.csv`), semillas 42–44 |

Los números de línea del formulario (23, 25, 26, 31, 35, 43, 44, 47) son los del paquete que recibió Dina. El script comprueba que corresponden exactamente a los ocho tramos excluidos.

## Criterio de estado (fijo, no depende de las métricas)

- **no verificable:** Dina no pudo decir qué es la línea («No se puede determinar») o la columna de interpretación repite lo que se ve. Son 5 tramos (N° 1, 2, 3, 6 y 8).
- **interpretación tentativa:** hay una interpretación, con confianza media o baja. Son 3 tramos:
  - N° 4: una acequia en un tramo (confianza media);
  - N° 5: posible drenaje antiguo, no confirmado (baja);
  - N° 7: posible quebrada seca en el primer tramo, hoy con una carretera (media).
- **verificado:** solo con confianza alta. Ninguno.
- Nada de esto prueba un uso prehispánico ni que un tramo fuera transitable en el pasado. Es revisión remota, sin campo.

## Decisión heurística

- Los ocho tramos **se mantienen** en la configuración final. Es la regla para tramos no verificables del registro interno de decisiones, fijada antes de revisar las imágenes.
- Solo se retiran en la sensibilidad H4m-excl8 (Tabla 11).
- No se cambia nada del modelo con esta auditoría.

## Lo que muestran los datos guardados

- **Tráfico en las celdas de los ocho tramos:** 2,83 % en la semilla 42 (287 celdas); 2,7 % de media en las semillas 42–44 (Tabla 11).
- **Al retirarlos (semillas 42–44):**
  - AUC 0,8291 → 0,8170;
  - ρ20 0,653 → 0,622;
  - captura en el 10 %: 44,0 → 38,9 %;
  - atractores y pares: 25 y 90 → 20 y 70.
- **Atractores:** 11 de los 25 están a menos de 1 km de esos tramos. Al retirarlos desaparecen 13 (10 de esos 11) y aparecen 8 nuevos.
- **Los % del formulario no son los de la configuración final.** Suman 40,1 %, el valor de la Tabla 11 para la penalización solo en el cauce principal. En la configuración final es 2,7 %.

## Texto propuesto para §4.8 (borrador hecho con IA: revisar y reescribir)

En el manuscrito candidato va en azul, dentro de `\propuesta{}`:

> Supplementary Material S2 lists the eight reaches by HydroRIVERS identifier, with the satellite-image review: along five of them the nature of the line could not be determined, and for the other three the interpretation is tentative (part of an irrigation ditch, a possible old drainage and a possible dry ravine now followed by a road). The imagery date was not recorded and none of the reaches was checked in the field. Over seeds 42--44, removing them lowers mean AUC from 0.8291 to 0.8170 and $\rho_{20}$ from 0.653 to 0.622, while 13 of the 25 attractors disappear and 8 new ones appear (Table~\ref{tab:hydro-sensitivity}).

## Pendiente de Dina

- Fecha de las imágenes de Google Earth que miró.
- Tramo N° 8 (línea 43): la columna «Qué es, en tu opinión» repite la observación. ¿Corresponde «no se puede determinar» o «la línea no corresponde a ningún cauce»?

## Archivos

- `S2_OCHO_TRAMOS.csv` — tabla S2 (una fila por tramo).
- `S2_ocho_tramos.gpkg` — capas: tramos con la revisión, **capa de incertidumbre separada**, atractores H4m y atractores sin los ocho tramos. Sin sitios arqueológicos.
- `S2_MAPA.pdf/png` — mismo marco, escala y colores con y sin los ocho tramos.
- `S2_log.txt` — comprobaciones. Script: `../code/s2_ocho_tramos.py`.
