#!/usr/bin/env bash
# Download nf-core/mhcquant's HepG2 test fixtures for --step ms_search
# regression testing, then convert to a clean mzML format that ALL
# search engines (MSFragger, Comet, Sage) can parse.
#
# Data source: https://github.com/nf-core/test-datasets/tree/mhcquant
# Original PRIDE: PXD009752 (lung adenocarcinoma immunopeptidome, Q Exactive)
#
# The ThermoRawFileParser mzMLs from nf-core/mhcquant lack instrument
# metadata that MSFragger's Batmass-IO requires.  Fix: double-convert
# mzML → mzXML → mzML through ProteoWizard msconvert, which fills in
# the missing instrumentConfiguration componentList.
#
# Requires: ProteoWizard Singularity container at $PWIZ_SIF
# (default: ~/xy86_scratch2/sanjay/singularity_cache/proteowizard.sif)
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
PWIZ_SIF="${PWIZ_SIF:-${HOME}/xy86_scratch2/sanjay/singularity_cache/proteowizard.sif}"

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

# Mini FASTA — first 5000 Swiss-Prot entries (~5 MB). Sized so:
#   - Sage loads in <16 GB RAM (vs. ~120 GB for the full 14 MB proteome)
#   - MokaPot gets enough PSMs below 10% FDR to train (500 proteins was
#     too few — mokapot errors "No PSMs found below the 'eval_fdr'")
# Override --search_fasta on the CLI for a full-proteome realistic run.
if [[ -s "human_swissprot_mini.fasta" ]]; then
    echo "[skip] human_swissprot_mini.fasta already present"
else
    awk '/^>/ {n++} n>5000 {exit} {print}' human_swissprot.fasta \
        > human_swissprot_mini.fasta
    echo "[build_test_ms_search_bundle] wrote mini FASTA with 5000 entries"
fi

# HepG2 HLA allele sheet + sample sheet
wget --quiet -O HepG2_allele_sheet.tsv "${BASE_URL}/HepG2_allele_sheet.tsv" 2>/dev/null || true
wget --quiet -O HepG2_sample_sheet.tsv "${BASE_URL}/HepG2_sample_sheet.tsv" 2>/dev/null || true

# Double-convert mzML → mzXML → mzML to produce clean files that all
# engines can read.  ThermoRawFileParser mzMLs have missing
# instrumentConfiguration metadata that breaks MSFragger (Batmass-IO NPE).
if [[ -f "${PWIZ_SIF}" ]]; then
    echo "[build_test_ms_search_bundle] double-converting mzMLs through ProteoWizard"
    for rep in 1 2 3; do
        base="HepG2_rep${rep}"
        if [[ -s "${base}_clean.mzML" ]]; then
            echo "[skip] ${base}_clean.mzML already present"
            continue
        fi
        # Step 1: mzML → mzXML
        singularity exec --bind "${OUT_DIR}:/data" "${PWIZ_SIF}" \
            wine msconvert "/data/${base}_small.mzML" --mzXML \
            --outdir /data --outfile "${base}.mzXML" 2>/dev/null
        # Step 2: mzXML → clean mzML
        singularity exec --bind "${OUT_DIR}:/data" "${PWIZ_SIF}" \
            wine msconvert "/data/${base}.mzXML" --mzML \
            --outdir /data --outfile "${base}_clean.mzML" 2>/dev/null
        echo "[build_test_ms_search_bundle] ${base}_clean.mzML ready"
    done
else
    echo "[WARN] ProteoWizard container not found at ${PWIZ_SIF}"
    echo "       MSFragger will NOT work with the raw ThermoRawFileParser mzMLs."
    echo "       Pull the container:  singularity pull proteowizard.sif docker://chambm/pwiz-skyline-i-agree-to-the-vendor-licenses:latest"
fi

# Build an ipg-format samplesheet pointing at the clean mzMLs.
cat > samplesheet.csv <<EOF
sample,ms_file
HepG2_rep1,${OUT_DIR}/HepG2_rep1_clean.mzML
HepG2_rep2,${OUT_DIR}/HepG2_rep2_clean.mzML
HepG2_rep3,${OUT_DIR}/HepG2_rep3_clean.mzML
EOF

echo "[build_test_ms_search_bundle] done"
ls -lh "${OUT_DIR}"/*.mzML "${OUT_DIR}"/*.fasta 2>/dev/null
