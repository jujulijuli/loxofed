#!/usr/bin/env bash
# DEMO COMPLETA con datos sinteticos: 3 sitios + nodo central, todo dentro de Docker y sin red.
# Sirve para comprobar que la imagen funciona y para ver que se genera en cada paso.
#   docker build -f docker/Dockerfile -t loxofed:0.7.6 .
#   ./docker/demo.sh [carpeta_de_trabajo]        # por defecto ./loxofed_demo
set -euo pipefail
# Una ruta relativa se interpreta desde donde ejecutas el script (no desde docker/)
ARG="${1:-loxofed_demo}"; case "$ARG" in /*) DEMO="$ARG" ;; *) DEMO="$PWD/$ARG" ;; esac
mkdir -p "$DEMO"; DEMO="$(cd "$DEMO" && pwd)"; cd "$(dirname "$0")"
export LOXO_N="${LOXO_N:-6000}"

# Configuracion comun a TODOS los sitios (misma huella). Valores reducidos para que sea rapido.
cat > "$DEMO/config.yml" <<'YML'
site_id: DEMO
paths: {pacientes: /data/pacientes.csv, enfermedades: /data/enfermedades.csv, output_dir: /output}
structure: {boot_R: 20, subsample: 30000}
model: {min_events_node: 30}
YML

echo "== 1) Cada sitio: simulate, qc, round1"
for S in A B C; do
  seed=$(printf '%d' "'$S")
  export LOXO_OUT="$DEMO/$S/out" LOXO_CONFIG="$DEMO/config.yml" LOXO_SEED="$seed"
  ./run_site.sh simulate
  export LOXO_DATA="$DEMO/$S/out/sim"
  # el site_id debe ser distinto en cada sitio: se cambia en una copia local de la configuracion
  mkdir -p "$DEMO/$S/cfg"; sed "s/site_id: DEMO/site_id: $S/" "$DEMO/config.yml" > "$DEMO/$S/cfg/config.yml"
  export LOXO_CONFIG="$DEMO/$S/cfg/config.yml"
  ./run_site.sh qc
  ./run_site.sh round1
done

echo "== 2) Central: consenso"
./run_central.sh consensus "$DEMO/config.yml" "$DEMO/central" "$DEMO/A/out/round1" "$DEMO/B/out/round1" "$DEMO/C/out/round1"

echo "== 3) Cada sitio: round2 sobre el consenso"
for S in A B C; do
  export LOXO_OUT="$DEMO/$S/out" LOXO_DATA="$DEMO/$S/out/sim" LOXO_CONFIG="$DEMO/$S/cfg/config.yml"
  ./run_site.sh round2 "$DEMO/central/consensus"
done

echo "== 4) Central: meta-analisis"
./run_central.sh meta "$DEMO/config.yml" "$DEMO/central" "$DEMO/A/out/round2" "$DEMO/B/out/round2" "$DEMO/C/out/round2"

echo
echo "Listo. Mira estos ficheros:"
echo "  $DEMO/central/consensus/consensus_structure.csv   (red consenso: arcos, direccion, fuerza)"
echo "  $DEMO/central/meta/meta_hazard_HR.csv             (HR de mortalidad meta-analizados)"
echo "  $DEMO/A/out/round1/ y round2/                     (lo UNICO que sale de un sitio)"
