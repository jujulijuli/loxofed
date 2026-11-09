#!/usr/bin/env bash
# Ejecuta una orden de loxofed en un sitio, sin red y con los datos en solo lectura.
#   ./run_site.sh qc                         # control de calidad de las entradas
#   ./run_site.sh round1                     # estructura + agregados
#   ./run_site.sh round2 /ruta/a/consensus   # ajuste sobre la estructura consenso
#   ./run_site.sh report [/ruta/a/consensus] # (re)genera informe_local.html desde las salidas ya hechas
#   ./run_site.sh simulate                   # CDM sintetico en $LOXO_OUT/sim para probar la instalacion
# Variables: LOXO_DATA (carpeta con los CSV), LOXO_CONFIG (config.yml), LOXO_OUT (salida),
#            LOXO_IMAGE (por defecto loxofed:0.7.6), LOXO_MEM (por defecto 24g),
#            LOXO_N / LOXO_SEED (solo para simulate)
set -euo pipefail
IMAGE="${LOXO_IMAGE:-loxofed:0.7.6}"
cmd="${1:-help}"; shift || true
OUT="${LOXO_OUT:?define LOXO_OUT con la carpeta de salida}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"

# --user: los ficheros de salida pertenecen a quien ejecuta, no a root ni a uid 1000
base=(docker run --rm --network none --memory "${LOXO_MEM:-24g}" --tmpfs /tmp
      --user "$(id -u):$(id -g)" -v "$OUT":/output)

case "$cmd" in
  simulate)
    "${base[@]}" "$IMAGE" simulate --out /output/sim --n "${LOXO_N:-5000}" --seed "${LOXO_SEED:-1}" ;;
  qc|round1|round2)
    DATA="${LOXO_DATA:?define LOXO_DATA con la carpeta de los CSV}"
    CONFIG="${LOXO_CONFIG:?define LOXO_CONFIG con la ruta de config.yml}"
    mounts=(-v "$(cd "$DATA" && pwd)":/data:ro -v "$(cd "$(dirname "$CONFIG")" && pwd)":/config:ro)
    args=("$cmd" --config "/config/$(basename "$CONFIG")")
    if [ "$cmd" = round2 ]; then
      cons="${1:?indica la carpeta consensus/ recibida del nodo central}"
      mounts+=(-v "$(cd "$cons" && pwd)":/consensus:ro); args+=(--consensus /consensus)
    fi
    "${base[@]}" "${mounts[@]}" "$IMAGE" "${args[@]}" ;;
  report)
    CONFIG="${LOXO_CONFIG:?define LOXO_CONFIG con la ruta de config.yml}"
    mounts=(-v "$(cd "$(dirname "$CONFIG")" && pwd)":/config:ro); args=(report --config "/config/$(basename "$CONFIG")")
    if [ -n "${1:-}" ]; then mounts+=(-v "$(cd "$1" && pwd)":/consensus:ro); args+=(--consensus /consensus); fi
    "${base[@]}" "${mounts[@]}" "$IMAGE" "${args[@]}" ;;
  *) "${base[@]}" "$IMAGE" help ;;
esac
