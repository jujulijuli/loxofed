#!/usr/bin/env bash
# Construye la imagen de loxofed y la empaqueta en UN fichero para repartirla a los sitios
# (cada sitio hace `docker load`; no necesita acceso a CRAN ni a GitHub).
#   ./docker/build_release.sh            # imagen linux/amd64 (servidores) -> loxofed_<version>_image.tar.gz
#   PLATFORM=linux/arm64 ./docker/build_release.sh   # solo para probar en un Mac con Apple Silicon
# IMPORTANTE: construir en amd64 aunque tu Mac sea arm64; una imagen arm64 en un servidor amd64 va emulada (muy lenta).
set -euo pipefail
cd "$(dirname "$0")/.."
VERSION="$(sed -n 's/^Version: *//p' DESCRIPTION)"
IMAGE="loxofed:${VERSION}"; PLATFORM="${PLATFORM:-linux/amd64}"
OUT="loxofed_${VERSION}_image.tar.gz"
echo ">> Construyendo ${IMAGE} para ${PLATFORM}"
docker buildx build --platform "$PLATFORM" --load -f docker/Dockerfile -t "$IMAGE" .
echo ">> Tests dentro de la imagen"
docker run --rm --platform "$PLATFORM" --entrypoint Rscript "$IMAGE" -e 'testthat::test_package("loxofed")'
echo ">> Guardando ${OUT}"
docker save "$IMAGE" | gzip > "$OUT"
if command -v md5sum >/dev/null; then md5sum "$OUT" | tee "${OUT}.md5"; else md5 "$OUT" | tee "${OUT}.md5"; fi
echo ">> Listo. Reparte ${OUT} y su .md5; cada sitio ejecuta:  gunzip -c ${OUT} | docker load"
