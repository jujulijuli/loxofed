#!/usr/bin/env bash
# Pasos del nodo central. Cada sitio entrega su carpeta round1/ (o round2/).
#   ./run_central.sh consensus  config.yml  salida/  siteA/round1 siteB/round1 ...
#   ./run_central.sh meta       config.yml  salida/  siteA/round2 siteB/round2 ...
set -euo pipefail
IMAGE="${LOXO_IMAGE:-loxofed:0.7.6}"
step="${1:?consensus|meta}"; config="${2:?config.yml}"; out="${3:?carpeta de salida}"; shift 3
mkdir -p "$out"
mounts=(); sites=""; i=0
for d in "$@"; do
  i=$((i+1)); mounts+=(-v "$(cd "$d" && pwd)":/in/s$i:ro); sites="$sites,/in/s$i"
done
docker run --rm --network none --tmpfs /tmp --user "$(id -u):$(id -g)" "${mounts[@]}" \
  -v "$(cd "$(dirname "$config")" && pwd)":/config:ro -v "$(cd "$out" && pwd)":/output \
  "$IMAGE" "central-$step" --config "/config/$(basename "$config")" \
  --sites "${sites#,}" --out /output
