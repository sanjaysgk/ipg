#!/usr/bin/env bash
# Build HepG2 chr22 end-to-end test fixture.
#
# Downloads ENCODE HepG2 paired-end polyA+ RNA-seq (ENCSR985KAT rep1)
# pre-aligned BAM, extracts chr22 reads to FASTQ, and pairs with the
# existing HepG2 clean mzMLs so the FULL pipeline (db_construct →
# ms_search) can be tested on one cell line.
#
# Re-uses the chr22 reference bundle from tests/data/ipg-db-constructions/
# (GRCh38.chr22.fa, gencode.v44.chr22.gtf, STAR index, GATK VCFs).
#
# Prerequisites:
#   - tests/data/ipg-db-constructions/ populated (run bin/build_test_bundle.sh)
#   - tests/data/ms_search/HepG2_rep*_clean.mzML (run bin/build_test_ms_search_bundle.sh)
#   - samtools in PATH (or pixi env)
#
# USAGE:
#   bash bin/build_e2e_hepg2_chr22.sh
#
# Output: tests/data/e2e_hepg2_chr22/
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
OUT_DIR="${REPO_ROOT}/tests/data/e2e_hepg2_chr22"
REF_DIR="${REPO_ROOT}/tests/data/ipg-db-constructions/reference"
VC_DIR="${REPO_ROOT}/tests/data/ipg-db-constructions/variant_calling"
MS_DIR="${REPO_ROOT}/tests/data/ms_search"
WORK_DIR="/fs04/scratch2/xy86/sanjay/e2e-build-tmp"

mkdir -p "${OUT_DIR}/fastq" "${WORK_DIR}"

echo "[e2e] Building HepG2 chr22 E2E fixture"

# ─────────────────────────────────────────────────────────────
# Step 1: Download ENCODE HepG2 GRCh38-aligned BAM (rep1, ~2.8 GB)
# ENCSR985KAT / ENCFF916YZY — polyA+ RNA-seq, Graveley lab UConn
# ─────────────────────────────────────────────────────────────
BAM_URL="https://encode-public.s3.amazonaws.com/2016/03/14/0cce4404-7222-44f1-ab2a-7c0759aeb620/ENCFF916YZY.bam"
BAM="${WORK_DIR}/ENCFF916YZY.bam"

if [[ -s "${BAM}" ]]; then
    echo "[skip] BAM already downloaded"
else
    echo "[e2e] Downloading ENCODE HepG2 GRCh38 BAM (~2.8 GB)..."
    wget -c -q --show-progress -O "${BAM}" "${BAM_URL}"
fi

# ─────────────────────────────────────────────────────────────
# Step 2: Extract chr22 paired reads to FASTQ
# ─────────────────────────────────────────────────────────────
R1_CHR22="${OUT_DIR}/fastq/HepG2_chr22_R1.fastq.gz"
R2_CHR22="${OUT_DIR}/fastq/HepG2_chr22_R2.fastq.gz"

if [[ -s "${R1_CHR22}" ]] && [[ -s "${R2_CHR22}" ]]; then
    echo "[skip] chr22 FASTQs already present"
else
    echo "[e2e] Indexing BAM..."
    samtools index "${BAM}"

    echo "[e2e] Extracting chr22 properly-paired reads..."
    # -f 3: read paired + both mapped;  chr22 region
    samtools view -b -f 3 "${BAM}" chr22 | \
        samtools sort -n -@ 2 - | \
        samtools fastq -@ 2 -1 "${R1_CHR22}" -2 "${R2_CHR22}" -0 /dev/null -s /dev/null -

    R1_COUNT=$(zcat "${R1_CHR22}" | awk 'END{print NR/4}')
    echo "[e2e] Extracted ${R1_COUNT} chr22 paired reads"
fi

# ─────────────────────────────────────────────────────────────
# Step 3: Create samplesheets + symlinks
# ─────────────────────────────────────────────────────────────
echo "[e2e] Writing samplesheets..."

# db_construct samplesheet (RNA-seq input)
cat > "${OUT_DIR}/samplesheet_rnaseq.csv" <<EOF
sample,fastq_1,fastq_2,strandedness
HepG2_chr22,${OUT_DIR}/fastq/HepG2_chr22_R1.fastq.gz,${OUT_DIR}/fastq/HepG2_chr22_R2.fastq.gz,auto
EOF

# ms_search samplesheet (MS input — same cell line!)
cat > "${OUT_DIR}/samplesheet_ms.csv" <<EOF
sample,ms_file
HepG2_rep1,${MS_DIR}/HepG2_rep1_clean.mzML
EOF

echo "[e2e] Done! Fixture at ${OUT_DIR}"
echo ""
echo "=== Run the full E2E test ==="
echo ""
echo "pixi run nextflow run . -profile pixi \\"
echo "  --input ${OUT_DIR}/samplesheet_rnaseq.csv \\"
echo "  --ms_input ${OUT_DIR}/samplesheet_ms.csv \\"
echo "  --genome_fasta ${REF_DIR}/GRCh38.chr22.fa \\"
echo "  --gtf ${REF_DIR}/gencode.v44.chr22.gtf \\"
echo "  --star_index ${REF_DIR}/star_index_chr22 \\"
echo "  --dbsnp ${VC_DIR}/dbsnp138.chr22.vcf.gz \\"
echo "  --known_indels ${VC_DIR}/known_indels.chr22.vcf.gz \\"
echo "  --known_mills ${VC_DIR}/Mills.chr22.vcf.gz \\"
echo "  --germline_resource ${VC_DIR}/small_exac_common_3.chr22.vcf.gz \\"
echo "  --search_fasta ${MS_DIR}/human_swissprot_mini.fasta \\"
echo "  --ms_engines comet,sage \\"
echo "  --skip_ms2rescore true --skip_integrate_engines true \\"
echo "  --outdir test_e2e_hepg2_chr22 \\"
echo "  -work-dir /fs04/scratch2/xy86/sanjay/nf-work"
