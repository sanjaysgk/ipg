#!/usr/bin/env bash
#
# fetch_test_bundle.sh — download and extract the chr22 test bundle into
# ${projectDir}/.test-bundle/ (gitignored). Run once per fresh clone.
#
# After this, `nextflow run . -profile test,docker --outdir results_test`
# works on any host with Nextflow + a container engine.
#
# Env overrides:
#   IPG_TEST_BUNDLE_URL  default points at v0.1.0-test-data Release asset
#   IPG_TEST_BUNDLE_DIR  default ${REPO_ROOT}/.test-bundle

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
URL="${IPG_TEST_BUNDLE_URL:-https://github.com/sanjaysgk/ipg/releases/download/v0.1.0-test-data/ipg-test-bundle-chr22.tar.gz}"
DEST="${IPG_TEST_BUNDLE_DIR:-${REPO_ROOT}/.test-bundle}"

log()  { printf '\033[1;34m[fetch_test_bundle]\033[0m %s\n' "$*"; }
fail() { printf '\033[1;31m[fetch_test_bundle] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

# Skip if already extracted.
if [[ -f "${DEST}/samplesheet_test.csv" ]]; then
    log "Bundle already present at ${DEST} — skipping download."
    log "To force re-download: rm -rf ${DEST}"
    exit 0
fi

mkdir -p "${DEST}"
log "Downloading bundle from ${URL}"
log "Destination: ${DEST}"

# Stream tarball straight into tar — no intermediate file, no /tmp pressure.
if ! curl -fL --retry 3 --retry-delay 2 "${URL}" | tar -xzf - -C "${DEST}" --strip-components=1; then
    fail "download failed from ${URL} — check network / URL"
fi

[[ -f "${DEST}/samplesheet_test.csv" ]] || fail "extraction completed but samplesheet_test.csv not found — bundle layout unexpected"

# Substitute __BUNDLE_DIR__ placeholder in the samplesheet with the actual
# extraction dir. nf-validation requires absolute paths in samplesheets, but
# the bundle is portable so we can't bake them in at build time.
ABS_DEST="$(readlink -f "${DEST}")"
sed -i.bak "s|__BUNDLE_DIR__|${ABS_DEST}|g" "${DEST}/samplesheet_test.csv"
rm -f "${DEST}/samplesheet_test.csv.bak"
log "Patched samplesheet paths -> ${ABS_DEST}"

log "Bundle ready at ${DEST}"
log "Run: nextflow run . -profile test,docker --outdir results_test"
