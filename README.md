# loxofed · análisis federado de la DBN de multimorbilidad LOXO

Paquete R + contenedor Docker que reproduce, en cada comunidad autónoma, el núcleo del análisis de Navarra
(red dinámica de dos cortes **t0→t1** + modelo de mortalidad) y exporta **solo agregados**.
Un nodo central construye la **estructura consenso** y **meta-analiza** los coeficientes.

Entrada por sitio: las dos tablas del CDM (`pacientes`, `enfermedades`). Sin red, sin microdatos fuera.

**Versión 0.7.6** (0.2.0: mortalidad con fechas exactas; 0.2.1: tope de personas por sitio; 0.3.0: pesos de muestreo e informe local; 0.3.1: informe con etiquetas cortas, bosque y red en capas; 0.4.0: meta-análisis REML + Hartung-Knapp e informe central; 0.5.0: informe central en Quarto; 0.6.0: calibración observadas/esperadas y comparador de ejecuciones; 0.7.0: efecto calendario en la mortalidad; 0.7.1: registro diagnóstico durante la pandemia). Novedad respecto a 0.1.0: la mortalidad (capa 2) usa por defecto las **fechas exactas** de muerte y de
diagnóstico (tiempo-persona), sin cortes anuales y sin dejar ninguna muerte fuera. La convención anterior sigue disponible
con `model.mortality: annual_next` para replicar el informe de Navarra.

## 1. Flujo (2 rondas)

```
 SITIO (×7)                               CENTRAL (p. ej. Aragón)                    SITIO (×7)
 ─────────                                ───────────────────────                    ─────────
 qc ──► round1 ──► r1_*.csv ─────────────► central-consensus ──► consensus/ ───────► round2 ──► r2_*.csv ──► central-meta
        · cohorte                           · nodos: prevalencia agrupada ≥ 1 %       · CPDs de incidencia          · efectos aleatorios
        · boot.strength                     · arcos: fuerza agrupada ≥ 0,60             sobre la red consenso          (DerSimonian-Laird)
        · prevalencia (60 cond.)            · dirección: precedencia por fechas       · mortalidad (tiempo-persona) · HR, I², τ²
        · precedencia por pares               agrupada (≥ 55 %)                       · personas-año y muertes
```

**Por qué dos rondas.** La estructura (qué arcos existen) solo se aprende con microdatos, y los coeficientes solo son
meta-analizables si se ajustan **sobre la misma estructura**. Ronda 1: cada sitio comparte la *fuerza* de cada arco (una
frecuencia de bootstrap, no datos). El centro fija el consenso. Ronda 2: cada sitio reajusta los modelos sobre esa red común y
comparte β + EE.

## 2. Entradas

Dos CSV delimitados por `|`, con cabecera. Los nombres de columna se configuran por sitio (`cdm.cols`).

| Tabla | Columna (por defecto) | Contenido | Granularidad |
|---|---|---|---|
| pacientes | `Id_pac` | identificador único | una fila por persona |
| | `F_Nac` | año de nacimiento (o fecha; se usa solo el año) | **año** |
| | `Sexo` | `H` (hombre) / `M` (mujer) según el modelo común; recodificado a Hombre/Mujer con `cdm.sex_map` (estricto: otro valor → ERROR en `qc`) | |
| | `F_Exitus` | fecha de muerte; vacía si sigue viva | **mes como mínimo, día ideal** (modo `exact`); el año basta solo con `annual_next` |
| | `F_Baja` *(opcional)* | baja administrativa | día o mes |
| enfermedades | `Id_pac`, `EnfC` | persona y condición (una de las 60 SNACK, `LOXO_DIS60`) | una fila por diagnóstico |
| | `F_EnfC` | fecha de diagnóstico | **día ideal, mes mínimo recomendable** |

Las fechas de muerte y de diagnóstico se usan a resolución completa en la mortalidad y en la validación de dirección
(qué enfermedad precede a cuál). `qc` avisa (`WARN`) si más del 5 % de las muertes (modo `exact`) o más del 20 % de los
diagnósticos son anuales. Las fechas de resolución parcial se imputan (año → 1 de julio; mes → día 15), pero con muerte
anual imputada un diagnóstico del mismo año posterior a esa fecha queda incoherente: pide mes o día.

## 3. Mortalidad (capa 2)

`model.mortality: exact` (por defecto)
- El seguimiento de cada persona va del 1 de enero de `year_first` hasta la **primera** de: muerte, baja administrativa o
  1 de enero de `year_last + 1`.
- Se parte en **cada diagnóstico incidente** y cada `time_split_months` meses (12 por defecto; solo sirve para que la edad avance).
- El estado de cada tramo es el vigente al inicio del tramo; el desenlace es la muerte al final. Modelo de **Poisson con offset**
  (log del tiempo expuesto): `exp(β)` = HR de tasa instantánea. La edad entra como spline de **nudos fijos** (edad en el punto medio del tramo).
- Se cuentan **todas** las muertes del seguimiento, incluida una muerte el propio día basal. Hay tests que lo verifican.
- `dx_lag_days` (0 por defecto): un diagnóstico cuenta como exposición solo `dx_lag_days` días después. Es una **sensibilidad**
  para diagnósticos codificados junto a la muerte (efecto, no causa). Debe fijarse **antes** de ver resultados y ser igual en
  todos los sitios (entra en la huella). Los diagnósticos prevalentes cuentan siempre.
- Coste: construir el tiempo-persona de ~11.000 personas tarda ~0,3 s y el ajuste ~1,3 s. Para ~150.000 personas se esperan
  del orden de 1,5-2 millones de filas; está por verificar con datos reales.

`model.mortality: annual_next`: replica el informe de Navarra (foto a 1 de enero de *t* → muerte en *t+1* entre quienes
sobreviven *t*, cloglog anual). Las muertes del propio año *t* no entran.

La **capa 1 (incidencia y estructura) sigue siendo anual** en ambos modos: `bnlearn` necesita una tabla discreta.

## 4. Qué sale de cada sitio

Solo estas tablas (más `manifest.json` con versión, huella de configuración y MD5 de cada fichero):

| Fichero | Contenido | Control de divulgación |
|---|---|---|
| `qc.csv` | resultados del control de calidad | conteos pequeños → `<10` |
| `r1_prevalence.csv` | prevalencia basal de las 60 condiciones | `n_pos` 0<n<10 → NA |
| `r1_nodes.csv`, `r1_edges.csv` | nodos locales y fuerza de arco A_t0→B_t1 | frecuencias de bootstrap, no hay celdas |
| `r1_direction.csv` | por par: A antes que B, B antes que A, empates | celdas pequeñas → NA y total anulado |
| `r1_meta.csv` | tamaños y versiones | |
| `r2_coefs.csv` | β, EE, n, eventos por coeficiente (incidencia y mortalidad) | modelos ajustados, sin filas |
| `r2_observed.csv` | `exact`: personas-año y muertes; `annual_next`: n en riesgo y muertes, por sexo × edad × año | `n_death` (y `n_risk`) 0<n<10 → NA |

`.lox_write()` **bloquea** la escritura de cualquier tabla con una columna de identificador o con más de 2 M de filas.
Umbral por defecto `disclosure.min_cell = 10`; el paquete avisa si se baja de 5.

## 5. Uso con Docker

Requisitos: Docker Desktop (o Docker Engine) y este directorio. La imagen es multiarquitectura (Mac Intel/Apple Silicon, Linux).

```bash
# 1) Construir la imagen (una vez; tarda unos minutos: compila los paquetes de R)
docker build -f docker/Dockerfile -t loxofed:0.7.0 .

# 2) Comprobar la imagen: tests del paquete dentro del contenedor
docker run --rm --entrypoint Rscript loxofed:0.7.0 -e 'testthat::test_package("loxofed")'

# 3) DEMO completa con datos sintéticos: 3 sitios + nodo central, sin red
./docker/demo.sh ./loxofed_demo      # si da "permission denied": chmod +x docker/*.sh
#    -> loxofed_demo/central/consensus/consensus_structure.csv   (red consenso)
#    -> loxofed_demo/central/meta/meta_hazard_HR.csv             (HR meta-analizados)
#    -> loxofed_demo/A/out/round1 y round2                       (lo único que sale de un sitio)

# 4) Con datos reales de UN sitio (p. ej. Navarra)
export LOXO_DATA=/ruta/a/csv LOXO_CONFIG=/ruta/a/config.yml LOXO_OUT=/ruta/a/salida
./docker/run_site.sh qc              # revisa qc.csv antes de seguir
./docker/run_site.sh round1          # -> $LOXO_OUT/round1/   (se envía al nodo central)
# central:
./docker/run_central.sh consensus config.yml central/ ARA/round1 NAV/round1 AND/round1 ...
# cada sitio:
./docker/run_site.sh round2 /ruta/a/central/consensus
# central:
./docker/run_central.sh meta config.yml central/ ARA/round2 NAV/round2 AND/round2 ...
```

Los contenedores corren con `--network none`, datos en solo lectura y como tu usuario (los ficheros de salida son tuyos).
`LOXO_MEM` (24g por defecto) limita la memoria. Configuración de ejemplo en `inst/config/config.example.yml`.
Sin Docker: `R CMD INSTALL .` y `Rscript -e 'loxofed::lox_cli()' round1 --config config.yml`.

Si ves un aviso de codificación ("cannot be translated from 'ANSI_X3.4-1968'"), es una configuración regional no UTF-8
fuera de Docker: `export LANG=en_US.UTF-8` (la imagen ya la trae).

### Qué debe ser idéntico en los 7 sitios
`cohort`, `model`, `structure` (salvo `seed`) y `consensus` entran en una **huella** (`lox_config_fingerprint`), que no depende
del orden de las claves. El central y la ronda 2 se **niegan** a combinar sitios con huellas o versiones del paquete distintas.
Lo propio de cada sitio (rutas, `site_id`, semilla, tope de personas) no entra.

**Spline de edad con nudos fijos** (`age_knots = 70, 80, 90`, frontera 60–105). Es imprescindible: con `ns(df = 4)` cada sitio
calcularía nudos propios y los coeficientes `edad_ns1..4` no serían comparables.

## 5b. Sitios grandes: tope de personas (`sample_n`)

Si un sitio es muy grande, se puede trabajar con una **muestra aleatoria de personas** (con todas sus enfermedades):

```yaml
cohort:
  sample_n: 50000      # ~ (null) = sin tope
  sample_seed: 2026
```

- El muestreo se hace sobre la cohorte **ya elegible** (>= `min_chronic` enfermedades prevalentes), así que 50.000 son
  50.000 personas analizables. Es reproducible con la semilla.
- `sample_n` y `sample_seed` son decisión del sitio: **no entran en la huella** (cada sitio puede usar un tope distinto).
  El tamaño real y el nº de elegibles quedan en `r1_meta.csv` y en el manifiesto de la ronda 2.
- El meta-análisis pondera por varianza: un sitio muestreado pesa menos, sin sesgo (muestra aleatoria simple de personas).
  Con 50.000 personas los eventos raros (p. ej. nodos con `min_events_node` < 100) pueden quedar fuera de ese sitio.
- **Pesos de muestreo.** El peso es `N elegibles / n muestreados` (se guarda como `sample_weight` en `r1_meta.csv` y en el
  manifiesto de la ronda 2). Solo se usa para **reescalar recuentos descriptivos** (`n_pos_est`, `n_est` en la prevalencia;
  `person_years_est`, `n_death_est` en la mortalidad observada). Se aplica **después** de suprimir celdas pequeñas, así que
  una celda suprimida sigue suprimida. Los **coeficientes no se ponderan**: con pesos constantes no cambian y los errores
  estándar saldrían falsos.
- Referencia (simulado, 207.790 elegibles; ronda 1 + 2; 2 núcleos): sin tope ~82 s y ~3,8 GB de RAM pico;
  con `sample_n: 50000` ~37 s y ~0,95 GB. La memoria crece ~linealmente con las personas, por eso el tope.

## 5c. Informe local

Al terminar la ronda 2 se escribe `<salida>/informe_local.html`: HTML autocontenido (R puro, sin Quarto ni red) con QC,
cohorte (y nota de muestreo), prevalencias, estructura (consenso si se dispone), HR de mortalidad, OR de incidencia y
mortalidad observada. Solo usa las tablas agregadas ya suprimidas, por lo que es tan compartible como `round1/` y `round2/`.
Se desactiva con `report: {auto: false}`. Para regenerarlo: `./docker/run_site.sh report [/ruta/consensus]`
(o `lox_report_local(output_dir, consensus_dir)` desde R). Muestra mortalidad **observada**, sin esperada ni
estandarización por edad.

## 5d. Meta-análisis (nodo central)

`lox_meta()` combina, por coeficiente, las estimaciones (β, EE) de los sitios en escala logarítmica:

- **Principal:** efectos aleatorios con τ² por **REML** e **IC95 % de Hartung-Knapp** (t con k−1 gl, truncado para no ser más estrecho que el normal),
  más **intervalo de predicción** (t con k−2 gl; k ≥ 3): dónde cabría el efecto en un sitio nuevo.
- **Sensibilidad:** DerSimonian-Laird (`beta_dl`, `tau2_dl`) y rango **dejando un sitio fuera** (`loo_min`, `loo_max`).
- **k mínimo:** `meta.min_sites` (3 por defecto). Con menos sitios, `status = "insuficiente (k<3)"` y no se da efecto combinado
  (las estimaciones por sitio quedan en `meta_site_coefs.csv`).
- **Heterogeneidad:** `tau2`, `I2` (Higgins), `Q` y `p_Q`.

Validación (`dev/validate_meta.R`, 3.000 réplicas, k = 7): cobertura del IC95 % del efecto medio con τ = 0,1 / 0,3 / 0,5 →
DL normal 0,90 / 0,88 / 0,89; REML normal 0,89 / 0,89 / 0,90; **REML + Hartung-Knapp 0,95 / 0,95 / 0,96**. El IC normal es demasiado
estrecho con pocos sitios, por eso HK es el principal. El intervalo de predicción cubre 0,85–0,95 (queda corto cuando τ es pequeño, porque
τ² es muy incierto con k = 7). Federación simulada de 7 sitios con heterogeneidad conocida (τ = 0,25; `dev/simulate_federation.R`, un
sitio con muestra aleatoria y dos con fechas mensuales): correlación entre el efecto combinado y la media verdadera 0,996, error absoluto medio 0,02,
τ estimada próxima a la verdadera, cobertura de HK 5/5.

Salidas (`meta/`): `meta_coefs.csv`, `meta_hazard_HR.csv`, `meta_site_coefs.csv`, `meta_sites.csv`, y **`informe_central.html`** (sitios, red consenso,
HR de mortalidad combinados, OR de incidencia con IC/predicción/I², arcos frágiles y diagramas de bosque por arco). Se regenera con
`central-report --out DIR`.

### Informe Quarto del meta-análisis (nodo central)

Además del HTML en R puro (`lox_report_central`), el nodo central puede generar un informe **Quarto** más completo (estilo informe de Navarra: métodos, forest por arco, mapas de calor de heterogeneidad, residuos sitio × arco, sensibilidad REML-HK vs DL, leave-one-out, red interactiva):

```r
loxofed::lox_render_central("<central_dir>")   # -> <central_dir>/informe_central_quarto.html
```

Requiere Quarto (>= 1.4) y, **solo en el nodo central**, los paquetes R: data.table, ggplot2, scales, igraph, knitr, ggraph, tidygraph, patchwork, ggrepel (opcionales: visNetwork, DT, htmltools, jsonlite). Los sitios no necesitan Quarto. El informe usa MathJax por CDN para las fórmulas: sin conexión se ven en crudo. `meta/meta_observed.csv` y `descriptive/` los genera `central-meta`/`central-consensus`.

### Calibración observadas / esperadas (0.6.0)

La ronda 2 (mortalidad con fechas exactas) devuelve en `r2_observed.csv` las muertes **esperadas** por el modelo ajustado (`n_exp`, y `n_exp_est` reescalada por el peso de muestreo) junto a las observadas, por sexo × grupo de edad × año. Las celdas suprimidas por control de revelación también suprimen la esperada. El informe local (sección 7) muestra O/E por grupo de edad, por año natural y por sexo, con IC95 % de Poisson exacto; el informe central Quarto lo compara entre sitios (sección «Calibración del modelo de mortalidad por sitio»). Por construcción el O/E global es 1; lo informativo es su forma por edad y año (el modelo no lleva tendencia temporal).

### Comparar dos ejecuciones (0.6.0)

```r
loxofed::lox_compare_runs("navarra_out_050", "navarra_out_full", labels = c("muestra 50k", "completa"))
# o: loxofed compare --a navarra_out_050 --b navarra_out_full --labels "muestra 50k,completa"
```

Genera `comparacion.html` en la carpeta B: cohorte, arcos que entran o salen del consenso (con acuerdo de Jaccard), HR de mortalidad con cambios de puesto, OR de incidencia y O/E. Útil para documentar el efecto de un filtro de datos, de un tope de muestreo o de un cambio de versión. Cada carpeta debe tener `round1/`, `round2/` y un `consensus/` (se busca en `consensus/`, `central_tmp/consensus/` y `central/consensus/`, o se pasa con `consensus_a` / `consensus_b`). Si las dos ejecuciones comparten datos de base, no son independientes: las diferencias son descriptivas.

### Efecto calendario en la mortalidad (0.7.0)

El modelo de mortalidad (`mortality = "exact"`) incluye por defecto, además de edad (spline), sexo y enfermedades, un **spline natural del año natural** con nudos fijos y comunes (`model$calendar$knots`, por defecto 2015 y 2018, con límites `year_first` y `year_last + 1`) y un **indicador de los años pandémicos** (`model$calendar$pandemic_years`, por defecto 2020 y 2021). Se desactiva con `model$calendar$enabled: false`. Los parámetros entran en la huella de configuración: todos los sitios deben usar los mismos.

Qué se devuelve y se muestra (transparencia):
- Coeficientes `anio_ns1..K` y `pandemia` en `r2_coefs.csv` (nodo `muerte`), meta-analizables como el resto.
- `r2_period_effect.csv`: tasa relativa de mortalidad por año (referencia = primer año) con IC95 % por método delta, separando la **tendencia** de la **tendencia + pandemia**.
- `r2_hr_sin_calendario.csv` y `n_exp0` en `r2_observed.csv`: el mismo modelo **sin** término calendario, para ver cuánto cambian HR y O/E al incluirlo.
- Informe local (secciones 7 y 8) e informe central Quarto («Efecto calendario en la mortalidad»): curvas por año, HR pandémico por sitio y combinado, y el efecto sobre los HR de las enfermedades.

Limitación declarada: en una cohorte cerrada el año natural y el tiempo desde el inicio del seguimiento van juntos, así que la **tendencia** suave mezcla cambio secular, selección de supervivientes y completitud del registro; se trata como ajuste. El indicador pandémico, al ser una desviación brusca respecto a esa tendencia, está mejor identificado. La capa de incidencia no lleva término calendario.

### Registro diagnóstico durante la pandemia (0.7.1)

`r2_dx_by_year.csv` devuelve, por enfermedad y año natural, las personas en riesgo al 1 de enero (vivas, en seguimiento y sin la enfermedad) y los **diagnósticos incidentes** del año (`TOTAL` suma todas las enfermedades; celdas < 10 suprimidas; `pandemia` marca los años de `model$calendar$pandemic_years`). `lox_dx_pandemic()` compara lo observado en los años pandémicos con lo esperado por la **tendencia de los años previos** (Poisson con tendencia lineal del año, extrapolado) y da el cociente observados/esperados con IC95 %. El informe local (sección 9) y el central Quarto («Registro diagnóstico durante la pandemia») lo muestran en conjunto y por enfermedad y sitio, y lo relacionan con el HR pandémico de la mortalidad. Un cociente < 1 sugiere infrarregistro (confinamiento, saturación asistencial, menor acceso); un déficit común a todos los sitios apunta a la pandemia, y uno solo de un sitio, a un problema local de codificación.

## 6. Qué se ha validado y qué no

- **Hecho (R 4.3.3, bnlearn 4.9.4):** 16 bloques de test, 96 comprobaciones, 0 fallos (`tests/testthat`). Incluyen un flujo
  federado completo de 2 sitios sintéticos, ambos modos de mortalidad, y el tiempo-persona caso por caso (muerte el día basal,
  baja, muerte tras la baja, muerte en el último día y fuera del seguimiento, diagnóstico incidente, contigüidad sin huecos,
  estados absorbentes). `R CMD check` sin errores; avisos pendientes: generar la ayuda `.Rd` con roxygen2.
- **Hecho:** la red consenso recupera las progresiones verdaderas con su dirección (HTA→IC, DM→ERC, depresión→demencia, FA→ictus).
- **Hecho (lógica de los scripts):** `demo.sh`, `run_site.sh` y `run_central.sh` se ejecutaron de punta a punta con un
  sustituto de `docker` que traduce los montajes a rutas locales (3 sitios + central, exit 0).
- **No hecho:** la **imagen Docker no se ha construido ni ejecutado** (sin daemon de Docker en el entorno de desarrollo).
  El `Dockerfile` usa `rocker/r-ver:4.4.2`. Primer paso: `docker build` y el paso 2-3 de arriba.
- **No hecho:** ejecución con los datos reales de Navarra (~150 k personas, ~56 nodos): tiempos y memoria por verificar.
- **No incluido (a propósito):** XGBoost, familias de variables DBN, Fase 2 estructural, sub-análisis Cox y el Shiny.

### Las tres formas de alinear la mortalidad, sobre una simulación con verdad conocida
`dev/check_death_alignment.R`, 40.000 personas simuladas en pasos mensuales. `OR_sim` es exp(β) en la escala logit anual del
simulador; en escala de tasa instantánea el HR verdadero es algo menor cuando la mortalidad basal es alta, así que la última
columna no tiene por qué llegar exactamente a `OR_sim`.

| Estado | OR_sim | anual *t→t+1* | mismo año | **exacto** |
|---|---|---|---|---|
| Fibrilación auricular | 1,22 | 1,20 | 1,21 | 1,22 |
| Enf. cerebrovascular | 1,35 | 1,25 | 1,29 | 1,31 |
| Enf. renal crónica | 1,49 | 1,34 | 1,39 | 1,43 |
| Demencia | 2,23 | 1,81 | 1,99 | 2,15 |
| Insuficiencia cardiaca | 2,01 | 1,68 | 1,83 | 1,89 |
| **Muertes contadas** (de 18.528) | | 17.371 | 18.528 | 18.528 |

La convención anual *t→t+1* atenúa los HR y pierde el 6 % de las muertes; el modo exacto cuenta todas y es el menos atenuado.
Es un resultado en datos **simulados** con tasa continua; en datos reales el tamaño del efecto dependerá de la cinética real.

## 7. Otras convenciones

1. **Bajas administrativas.** Si se informa `cdm.cols.baja`, se excluye a quien tiene baja en el año basal o antes
   (`exclude_baja_baseline`) y se censura al llegar la baja (en `exact` por fecha; en `annual_next` por año). En Navarra no hay
   baja, así que este camino solo está cubierto por casos construidos a mano, no por un dato real.
2. **Nodos por sitio.** Cada sitio aprende la red con los nodos que superan el 1 % localmente; el consenso solo usa pares
   donde ambos nodos eran nodos locales. Un nodo raro en un sitio no puede aportar ni vetar arcos allí.
3. **Dirección.** Los empates de fecha no cuentan como orden. Se aplica el mismo criterio que en el informe (`dir_min = 0,55`)
   sobre los conteos agrupados de todos los sitios.
4. **Edad.** La cohorte filtra por `año basal − año de nacimiento` (igual que el informe). El modelo `exact` usa la edad en el
   punto medio de cada tramo suponiendo nacimiento a mitad de año; la capa 1 anual usa la convención del informe.

## 0.7.3 — informe local para lectores no especializados

El informe local incluye un bloque «Resumen para lectores no especializados» (cifras dinámicas), un recuadro
«Cómo leerlo» tras cada sección, un recuadro inicial sobre el proceso federado y un glosario. El contenido técnico no cambia. Se exporta `lox_calendar_basis`.
Nota de interpretación: con el efecto calendario, la razón O/E por año se aleja de 1 en el primer y último año
(se comprobó que las ecuaciones de score del ajuste se cumplen y que añadir nudos no lo corrige); es una limitación
de la cohorte cerrada y se advierte en el informe.

## 0.7.4 — QC de baja frente a muerte

Regla común a todos los sitios (decidida): la fecha de baja **censura** el seguimiento y una muerte solo cuenta si ocurre en la fecha de baja o antes.
`qc` ahora informa, si el sitio tiene `baja`: personas con ambas fechas, muertes con baja anterior/igual/posterior, **muertes perdidas por la regla**
(n y % de las muertes del seguimiento, y cuántas caen en los 30 días tras la baja) y diagnósticos posteriores a la baja. Se avisa (WARN) si se pierde más del 1% de las muertes.
El informe local lo describe en la sección 1.

## 0.7.5 — ronda 2 ampliada y análisis del nodo central

- **Baja:** `F_Baja` es ahora la columna por defecto (`baja: ~` la desactiva).
- **Estructura provisional:** si el consenso se calculó con un solo sitio, los informes lo dicen.
- **Errores estándar robustos por persona** (sandwich agrupado) en la capa de incidencia: con varias filas por persona, los EE del GLM eran optimistas y sobreponderaban a los sitios grandes en el metaanálisis.
- **Contraste de edad** (`edad_70a80`, log-OR/log-HR de pasar de 70 a 80 años, con su EE por método delta) en incidencia y mortalidad: resume la spline en un número meta-analizable.
- **HR por sexo** (`r2_hr_por_sexo.csv`): el modelo de mortalidad ajustado por separado en hombres y mujeres; el centro metaanaliza cada sexo y el contraste hombre − mujer (`meta_sex.csv`).
- **Tablas 2×2 por arco** (`r2_arc_tables.csv`) e **incidencia por año × sexo × grupo de edad** (`r2_incidence_strata.csv`), con celdas < 10 suprimidas.
- **Nodo central:** `lox_network_metrics()` (grado, fuerza, centralidad armónica, intermediación: disparadores y receptores), `lox_communities()` (Louvain/Leiden; requiere `igraph`) y el índice de Rand ajustado entre sitios. El informe Quarto añade centralidad, comunidades, sexo, edad y sexo sobre la incidencia y OR crudos con Mantel-Haenszel.
- Con solo 2 sitios (banco de pruebas) hay que poner `meta: {min_sites: 2}` en la configuración del centro; los resultados se leen entonces como ensayo del procedimiento.

## 0.7.6 — enfermedades de un solo sexo
En la capa 1 (incidencia), si casi todos los casos de un nodo son de un sexo (el otro sexo tiene ≤ 5 casos o < 1 %), el modelo se ajusta solo en ese sexo y **no se estima el efecto del sexo** (evita la separación perfecta, β≈18, p. ej. en enfermedades de próstata). `r2_coefs.csv` lleva la columna `restriccion` (`solo_Hombre`/`solo_Mujer`) y el informe local lo avisa. La mortalidad no cambia. Todos los sitios deben usar la misma versión.
