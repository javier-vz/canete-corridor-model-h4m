# Barrido ε/κ en el Ryzen (JCAA-80, recálculo final)

Este paquete corre las **1,210 corridas** del barrido de ε y κ con la configuración aprobada (H4m). Es la misma malla y son las mismas semillas que el barrido original. Usa el kernel v2 sin cambios.

**Hay dos momentos:**
- **Hoy, la prueba de entorno** (unos 5 minutos). Solo verifica que tu computadora da exactamente los mismos números que aquí. No produce resultados del artículo.
- **Después, el barrido** (4–6 horas). Solo cuando te avise que ChatGPT validó la matriz de recálculo.

## 1. Preparar (una vez)

1. Descomprime `JCAA80_barrido_H4m_Ryzen.zip` dentro de `Documents\pre_paper_v3\`. Se crea la carpeta `barrido_H4m`; no toca nada de lo que ya tienes.
2. Abre **Anaconda Prompt** y entra a la carpeta:
   ```
   cd Documents\pre_paper_v3\barrido_H4m
   ```
3. Comprueba que están las librerías:
   ```
   python -c "import numpy, scipy, pandas, networkx; print('OK')"
   ```
   Si falta alguna: `pip install numpy scipy pandas networkx`.

## 2. Prueba de entorno (hoy)

```
python final_sweep_eps_kappa.py --prueba
```

**Qué revisa:**
- que los archivos sean los del paquete;
- que el estado del modelo sea idéntico: 25 atractores y 90 pares;
- que la corrida de la semilla 42 dé exactamente el mismo F, L, AUC y ρ₂₀ que aquí.

**Al final dice una de dos cosas:**
- **PRUEBA OK:** se crea `prueba_OK.json`. Mándamelo; trae las versiones y el tiempo por corrida.
- **PRUEBA FALLÓ:** no sigas. Mándame una foto o una copia de toda la pantalla.

## 3. Barrido (solo cuando te avise)

```
python final_sweep_eps_kappa.py --lanzar
```

- Usa todos los núcleos menos uno. Para usar menos procesos, por ejemplo si la PC se calienta: `--lanzar --workers 8`.
- **Antes de lanzar:**
  - deja la laptop **enchufada**;
  - en Configuración → Sistema → Energía, pon "Suspender: nunca" mientras corre;
  - cierra programas pesados.
- **Solo lanza si la prueba de este mismo paquete salió OK.** Si cambias algún archivo de la carpeta, el script se niega a lanzar y te pide repetir la prueba.
- **Si se interrumpe** (Ctrl+C, apagón, reinicio): escribe el mismo comando y sigue donde quedó.
  - Cada corrida se guarda por separado y se revisa su contenido.
  - Solo se vuelven a correr las que faltan o están dañadas.
  - Las dañadas se mueven a `resultados_barrido_H4m\invalidos\`; no se borran.

Para ver cuántas van bien en cualquier momento:

```
python final_sweep_eps_kappa.py --verificar
```

## 4. Al terminar

```
python final_sweep_eps_kappa.py --resumen
```

Crea `para_enviar_barrido_H4m.zip`. **Ese es el archivo que me mandas.**
- Solo lo crea si están las 1,210 corridas exactas, todas válidas y sin repetidas.
- Si falta alguna, te lo dice: vuelve a `--lanzar` y se completan solas.

## Qué hay en la carpeta

| Archivo | Qué es |
|---|---|
| `final_sweep_eps_kappa.py` | el script del barrido (prueba, lanzar, resumen) |
| `ca_kernel_v2.py` | el kernel del modelo, sin cambios (el mismo de `pre_paper_v3`) |
| `final_h4m_state.py` | arma la configuración H4m sobre el kernel |
| `datos/` | grilla, celdas de río y zonas P2 |
| `referencia/` | valores de la semilla 42 calculados aquí, para la prueba |
| `MANIFEST.sha256` | huellas de los archivos; la prueba las verifica |
