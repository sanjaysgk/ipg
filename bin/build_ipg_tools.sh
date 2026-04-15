#!/usr/bin/env bash
# Compile the kescull / Li-Purcell Lab C tools that the IPG pipeline
# expects on PATH when running under -profile pixi (no containers).
#
# The Dockerfile at containers/ipg-tools/Dockerfile builds these for
# the singularity/docker container. This script is the pixi equivalent.
# All six binaries land in bin/ and are gitignored.
#
# Run once per clone:
#   bash bin/build_ipg_tools.sh
#
# Re-run after pulling a new version of containers/ipg-tools/src/.
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SRC_DIR="${REPO_ROOT}/containers/ipg-tools/src"
OUT_DIR="${REPO_ROOT}/bin"

[[ -d "${SRC_DIR}" ]] || { echo "missing source dir: ${SRC_DIR}" >&2; exit 2; }
mkdir -p "${OUT_DIR}"

cd "${SRC_DIR}"

CC=${CC:-gcc}
CFLAGS=${CFLAGS:-"-O2 -Wall"}

echo "[build_ipg_tools] using ${CC} ${CFLAGS}"

${CC} ${CFLAGS} -o "${OUT_DIR}/curate_vcf"       curate_vcf.c
${CC} ${CFLAGS} -o "${OUT_DIR}/revert_headers"   revert_headers.c
${CC} ${CFLAGS} -o "${OUT_DIR}/alt_liftover"     alt_liftover.c
${CC} ${CFLAGS} -o "${OUT_DIR}/triple_translate" triple_translate.c
${CC} ${CFLAGS} -o "${OUT_DIR}/squish"           squish.c -pthread -lm
${CC} ${CFLAGS} -o "${OUT_DIR}/origins"          origins.c cJSON.c -lcurl -lm

echo "[build_ipg_tools] built:"
ls -la "${OUT_DIR}"/{curate_vcf,revert_headers,alt_liftover,triple_translate,squish,origins}
