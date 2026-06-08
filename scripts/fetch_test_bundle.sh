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

# Download only if not already extracted; patch step below is idempotent
# and always runs so re-runs heal half-patched / re-cloned trees.
if [[ -f "${DEST}/samplesheet_test.csv" ]]; then
    log "Bundle already present at ${DEST} — skipping download."
    log "To force re-download: rm -rf ${DEST}"
else
    mkdir -p "${DEST}"
    log "Downloading bundle from ${URL}"
    log "Destination: ${DEST}"
    if ! curl -fL --retry 3 --retry-delay 2 "${URL}" | tar -xzf - -C "${DEST}" --strip-components=1; then
        fail "download failed from ${URL} — check network / URL"
    fi
    [[ -f "${DEST}/samplesheet_test.csv" ]] || fail "extraction completed but samplesheet_test.csv not found — bundle layout unexpected"
fi

# Absolutise fastq paths in the samplesheet. nf-validation requires
# absolute paths, but the bundle is portable so paths are shipped
# relative to the bundle root (e.g. `fastq/test_R1.fastq.gz`). Rewrite
# any non-absolute fastq_1 / fastq_2 entry to point under ABS_DEST.
ABS_DEST="$(readlink -f "${DEST}")"
python3 - "${DEST}/samplesheet_test.csv" "${ABS_DEST}" <<'PY'
import csv, os, sys
path, root = sys.argv[1], sys.argv[2]
with open(path, newline='') as f:
    rows = list(csv.DictReader(f))
fields = list(rows[0].keys()) if rows else []
for r in rows:
    for col in ('fastq_1', 'fastq_2'):
        v = r.get(col, '')
        if v and not os.path.isabs(v):
            r[col] = os.path.join(root, v)
with open(path, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
PY
log "Patched samplesheet paths -> ${ABS_DEST}"

log "Bundle ready at ${DEST}"
log "Run: nextflow run . -profile test,docker --outdir results_test"
