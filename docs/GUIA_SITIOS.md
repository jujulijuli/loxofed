# loxofed · Guía para los sitios participantes

Análisis federado de multimorbilidad (red bayesiana dinámica + mortalidad). **Los datos individuales no salen del sitio**:
el contenedor se ejecuta **sin red** y solo escribe tablas agregadas, con las celdas pequeñas suprimidas.

## 0. Resumen del recorrido

| Paso | Quién | Qué | Resultado |
|---|---|---|---|
| 1 | sitio | instalar la imagen Docker | `loxofed:<versión>` |
| 2 | sitio | preparar los 2 CSV y el `config.yml` | |
| 3 | sitio | `qc` (control de calidad) | `qc.csv`; corregir los ERROR |
| 4 | sitio | `round1` | carpeta `round1/` → **se envía** |
| 5 | central | consenso con las `round1/` de todos | carpeta `consensus/` → se reparte |
| 6 | sitio | `round2 <consensus>` | carpeta `round2/` → **se envía**; `informe_local.html` |
| 7 | central | meta-análisis de las `round2/` | resultados finales |

## 1. Requisitos
- Docker (Linux, macOS, Windows con WSL2). Para los scripts `.sh` en Windows: WSL o Git Bash.
- RAM para Docker (en Docker Desktop: Settings → Resources → Memory, y «Apply & restart»). Referencia medida en Navarra
  (≈140.000 personas, sin muestreo): la ronda 2 fallaba con 7,7 GB y funcionó con 15,4 GB. Para sitios de ~300.000 personas se recomiendan **32 GB**;
  con `sample_n: 50000` bastan unos 8 GB. En Linux sin Docker Desktop, `run_site.sh` limita el contenedor a `LOXO_MEM` (por defecto 24g).
- Tiempo: en Navarra, ronda 1 ≈ 14 min y ronda 2 ≈ 2 min; en un sitio del doble de tamaño, cuenta aproximadamente el doble.
- Disco: ~3 GB para la imagen y unos MB para las salidas.
- **No hace falta internet** para ejecutar (solo para descargar la imagen si no se recibe como fichero).

## 2. Instalar la imagen
Se recibe `loxofed_<versión>_image.tar.gz` y su `.md5`. Comprobar y cargar:
```bash
md5sum loxofed_<versión>_image.tar.gz        # (macOS: md5 ...)  debe coincidir con el .md5
gunzip -c loxofed_<versión>_image.tar.gz | docker load
docker run --rm loxofed:<versión> version     # imprime la versión
```
Todos los sitios deben usar **la misma versión**: el nodo central rechaza combinar versiones distintas.
(Alternativa: construirla con `docker build -f docker/Dockerfile -t loxofed:<versión> .` en la raíz del paquete.)

## 3. Preparar los datos
Dos CSV delimitados por `|`, con cabecera, en una carpeta **local** (no sincronizada con OneDrive/Dropbox ni
dentro de una carpeta compartida con terceros). El contenedor los monta en **solo lectura**.

| Tabla | Columna (por defecto) | Contenido | Granularidad mínima |
|---|---|---|---|
| pacientes | `Id_pac` | identificador único | una fila por persona |
| | `F_Nac` | nacimiento | **año** |
| | `Sexo` | sexo | `H` / `M` (modelo común) |
| | `F_Exitus` | fecha de muerte (vacía si vive) | **mes**; **día ideal** |
| | `F_Baja` *(opcional)* | baja administrativa | mes o día |
| enfermedades | `Id_pac`, `EnfC` | persona y condición (código del catálogo de 60) | una fila por diagnóstico |
| | `F_EnfC` | fecha de diagnóstico | **mes**; **día ideal** |

Formatos de fecha admitidos: `YYYY-MM-DD`, `YYYY-MM` o `YYYY`. Cuanto más fina la fecha, mejor la validación de la
dirección de los arcos y la mortalidad. `qc` avisa si más del 5 % de las muertes o el 20 % de los diagnósticos son anuales.

## 4. Configurar (`config.yml`)
Copia `inst/config/config.example.yml` y edita **solo la sección «propio del sitio»**:
`site_id` y, si hiciera falta, las rutas. Con el modelo de datos común (mismos nombres de columna, `Sexo` = H/M) no hay que cambiar nada más. La sección «común» debe quedar **idéntica**
en todos los sitios (el central compara una huella y se niega a combinar configuraciones distintas).

- **`Sexo`:** el modelo común lo codifica `H` (hombre) / `M` (mujer). El mapeo por defecto es estricto: si aparece
  cualquier otro valor, `qc` da ERROR en lugar de asignarlo en silencio.
- **`cohort.sample_n`** (opcional, propio del sitio): tope de personas. Muestra aleatoria de la cohorte elegible, con todas
  sus enfermedades. Por defecto `~` = sin tope (recomendado: pruébalo primero así). Si no cabe en memoria, usa un valor como `150000` (o `50000`);
  la muestra es reproducible (`sample_seed`), los modelos se ajustan sin ponderar y los recuentos descriptivos se reescalan al tamaño real.
  El tamaño real queda en `r1_meta.csv`. Avisa al nodo central si usas muestreo.
- **Fecha de baja (`cdm.cols.baja`):** por defecto es `F_Baja`, columna del modelo común (puede ir **vacía**). Si tu tabla no la tiene, pon `baja: ~`; si se llama distinto, indica su nombre.
  La regla del análisis es común a todos los sitios: la baja **censura** el seguimiento y una muerte solo cuenta si ocurre en la fecha de baja o antes. `qc` informa de cuántas muertes se pierden
  por tener la baja anterior (y cuántas dentro de los 30 días siguientes); si el porcentaje es alto, avisa antes de continuar.

## 5. Ejecutar
```bash
cd loxofed                                   # carpeta del paquete
export LOXO_DATA=/ruta/a/tus/csv             # carpeta con los 2 CSV
export LOXO_CONFIG=/ruta/a/config.yml
export LOXO_OUT=/ruta/a/salida               # donde se escribe todo
export LOXO_IMAGE=loxofed:<versión>

./docker/run_site.sh qc                      # 1) control de calidad: no debe haber ERROR
./docker/run_site.sh round1                  # 2) estructura + agregados  -> $LOXO_OUT/round1
#   ... enviar $LOXO_OUT/round1 al nodo central y esperar la carpeta consensus/ ...
./docker/run_site.sh round2 /ruta/a/consensus   # 3) ajuste -> $LOXO_OUT/round2 + informe_local.html
```
Las variables valen para esa ventana de terminal. Memoria del contenedor: `LOXO_MEM` (por defecto 24g; baja a la RAM real,
p. ej. `export LOXO_MEM=8g`).

**Tiempos de referencia** (Navarra, `sample_n: 50000`, MacBook Pro): `qc` segundos · `round1` ~14 min · `round2` ~1–2 min.
La ronda 1 domina por el bootstrap de la estructura.

## 6. Qué enviar y qué no
Enviar **solo** las carpetas `round1/` y `round2/` completas (con su `manifest.json`), por el canal acordado. Contienen
tablas agregadas con celdas < 10 suprimidas; el paquete se niega a escribir columnas de identificador o tablas > 2 M filas.
Antes de enviar:
- [ ] `round1/` y `round2/` contienen solo `.csv` y `manifest.json`.
- [ ] `qc.csv` sin ERROR.
- [ ] Revisión rápida de `informe_local.html` (se puede compartir; usa las mismas tablas agregadas).

**No enviar:** los CSV de entrada, el contenedor con datos montados, ni ningún fichero fuera de esas carpetas.

## 7. Problemas frecuentes
| Síntoma | Causa y solución |
|---|---|
| `sexo_sin_mapear` > 0 (ERROR) | hay valores de `Sexo` distintos de H/M: corrige la extracción (o, si es un caso conocido, avisa para ampliar `sex_map`) |
| `docker build ... requires 1 argument` | falta el punto final: `docker build -f docker/Dockerfile -t loxofed:<v> .` (con espacio antes del punto) |
| «huella distinta» | la sección común del config no coincide con la del central; usa la plantilla sin tocarla |
| versiones distintas de paquete | todos los sitios con la misma imagen (`loxofed:<versión>`) |
| no encuentra los CSV | carpeta en OneDrive/nube: cópialos a una carpeta local; revisa que `LOXO_DATA` es la carpeta, no el fichero |
| código de salida 137 / «killed» / sin memoria | Docker se quedó sin RAM: sube la memoria de Docker Desktop (y «Apply & restart»), ajusta `LOXO_MEM` a la RAM real o usa `sample_n` |
| `exitus_perdidos_por_baja_anterior` (WARN) | más del 1% de las muertes del seguimiento tienen la baja anterior y no se cuentan (regla de censura); revisa el origen de ambas fechas y avisa |
| `dx_posterior_al_exitus` (WARN) | diagnósticos fechados tras la muerte: se ignoran; avisa si son muchos |
| acentos raros en el informe | la imagen ya fija UTF-8; no la modifiques |

## 8. Para el nodo central
```bash
./docker/run_central.sh consensus config.yml salida/ siteA/round1 siteB/round1 ...   # -> salida/consensus/
./docker/run_central.sh meta      config.yml salida/ siteA/round2 siteB/round2 ...   # -> salida/meta/
```
`consensus/` se reparte a todos los sitios. `meta/` contiene `meta_coefs.csv` (efecto combinado REML + Hartung-Knapp, intervalo de
predicción, τ², I², sensibilidad), `meta_hazard_HR.csv`, `meta_site_coefs.csv` y `meta_sites.csv`; y se escribe `informe_central.html`
(`run_central.sh` lo deja en la carpeta de salida). Con menos de 3 sitios no se da resultado combinado (`meta.min_sites`).

## 9. Qué hace el análisis (resumen)
Red bayesiana dinámica de dos capas sobre 60 condiciones crónicas (catálogo SNACK), con **edad y sexo como covariables**:
capa 1 = incidencia anual de cada condición dada la estructura consenso (regresión logística), capa 2 = mortalidad con
fechas exactas (tiempo-persona, Poisson). Los sitios no comparten pacientes ni modelos individuales, sino coeficientes
(β y error estándar) que el central combina por meta-análisis de efectos aleatorios.
