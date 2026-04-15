#!/usr/bin/env bash
# Download nf-core/mhcquant's HepG2 test fixtures for --step ms_search
# regression testing.
#
# Data source: https://github.com/nf-core/test-datasets/tree/mhcquant
# Original PRIDE: PXD009752 (lung adenocarcinoma immunopeptidome, Q Exactive)
#
# USAGE:
#   bash bin/build_test_ms_search_bundle.sh
#
# Output lives under tests/data/ms_search/ which is .gitignored — fetch
# on demand, never commit the ~160 MB of mzMLs.
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT_DIR="${REPO_ROOT}/tests/data/ms_search"
BASE_URL="https://raw.githubusercontent.com/nf-core/test-datasets/mhcquant/testdata"

mkdir -p "${OUT_DIR}"
cd "${OUT_DIR}"

echo "[build_test_ms_search_bundle] downloading HepG2 mzMLs + Swiss-Prot FASTA into ${OUT_DIR}"

# Three HepG2 replicate mzMLs (~50 MB each)
for f in HepG2_rep1_small HepG2_rep2_small HepG2_rep3_small; do
    if [[ -s "${f}.mzML" ]]; then
        echo "[skip] ${f}.mzML already present"
    else
        wget --quiet --show-progress "${BASE_URL}/${f}.mzML"
    fi
done

# Human Swiss-Prot FASTA (~14 MB)
if [[ -s "human_swissprot.fasta" ]]; then
    echo "[skip] human_swissprot.fasta already present"
else
    wget --quiet --show-progress -O human_swissprot.fasta \
        "${BASE_URL}/2025_02_human_swissprot.fasta"
fi

# HepG2 HLA allele sheet + sample sheet
wget --quiet -O HepG2_allele_sheet.tsv "${BASE_URL}/HepG2_allele_sheet.tsv"
wget --quiet -O HepG2_sample_sheet.tsv "${BASE_URL}/HepG2_sample_sheet.tsv"

# Build an ipg-format samplesheet pointing at the three local mzMLs.
cat > samplesheet.csv <<EOF
sample,ms_file
HepG2_rep1,${OUT_DIR}/HepG2_rep1_small.mzML
HepG2_rep2,${OUT_DIR}/HepG2_rep2_small.mzML
HepG2_rep3,${OUT_DIR}/HepG2_rep3_small.mzML
EOF

echo "[build_test_ms_search_bundle] done"
ls -lh "${OUT_DIR}"
