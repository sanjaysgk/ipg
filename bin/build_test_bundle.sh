#!/usr/bin/env bash
#
# build_test_bundle.sh — build a tiny chr22-subset test data bundle for
# sanjaysgk/ipg integration tests.
#
# Derives the bundle from the full Monash M3 reference bundle and a real
# D100-liver FASTQ sample by subsetting to chromosome 22. Writes everything
# to $TEST_BUNDLE_DIR (default /fs04/scratch2/xy86/sanjay/ipg-test-data).
#
# The bundle contains:
#
#     $TEST_BUNDLE_DIR/
#     ├── samplesheet_test.csv
#     ├── fastq/
#     │   ├── test_R1.fastq.gz             subsampled from D100-liver R1
#     │   └── test_R2.fastq.gz             subsampled from D100-liver R2
#     ├── reference/
#     │   ├── GRCh38.chr22.fa              samtools faidx chr22
#     │   ├── GRCh38.chr22.fa.fai
#     │   ├── GRCh38.chr22.dict            samtools dict
#     │   ├── gencode.v44.chr22.gtf        awk '$1 == "chr22"'
#     │   ├── gencode_assembly.chr22.bed   awk '$1 == "chr22"' (for RSeQC)
#     │   └── star_index_chr22/            STAR --runMode genomeGenerate
#     └── variant_calling/
#         ├── dbsnp138.chr22.vcf.gz   + .tbi   bcftools view -r chr22
#         ├── known_indels.chr22.vcf.gz + .tbi
#         ├── Mills.chr22.vcf.gz + .tbi
#         └── small_exac_common_3.chr22.vcf.gz + .tbi
#
# Usage (from the repo root, inside the pixi env):
#     pixi run bash bin/build_test_bundle.sh
#
# Or override locations with env vars:
#     TEST_BUNDLE_DIR=/some/path \
#     REFERENCE_DIR=/some/ref \
#     FASTQ_R1=/some/r1.fq.gz \
#     FASTQ_R2=/some/r2.fq.gz \
#     pixi run bash bin/build_test_bundle.sh
#
# Re-running is idempotent: completed stages are skipped on subsequent
# invocations. Delete a sub-directory of $TEST_BUNDLE_DIR to force a re-run
# of just that stage.

set -euo pipefail

# ---- Configuration ---------------------------------------------------------
SOURCE_ROOT="${SOURCE_ROOT:-/fs04/scratch2/xy86/sanjay/ATLANTIS/RNAseq/Analysis/Cryptic}"
REFERENCE_DIR="${REFERENCE_DIR:-${SOURCE_ROOT}/GRCh38}"
VARIANT_DIR="${VARIANT_DIR:-${SOURCE_ROOT}/variant_calling_resources}"
RAW_DIR="${RAW_DIR:-/fs04/scratch2/xy86/sanjay/ATLANTIS/RNAseq/Raw_Data}"
FASTQ_R1="${FASTQ_R1:-${RAW_DIR}/D100-liver_S2_L001_R1_001.fastq.gz}"
FASTQ_R2="${FASTQ_R2:-${RAW_DIR}/D100-liver_S2_L001_R2_001.fastq.gz}"

TEST_BUNDLE_DIR="${TEST_BUNDLE_DIR:-/fs04/scratch2/xy86/sanjay/ipg-test-data}"
CHR="${CHR:-chr22}"
SUBSAMPLE_READS="${SUBSAMPLE_READS:-200000}"    # 200k paired; ~3-5k will land on chr22 — small enough that STAR finishes in ~3 min on Lustre
THREADS="${THREADS:-4}"

# ---- Helpers ---------------------------------------------------------------
log()  { printf '\033[1;34m[build_test_bundle]\033[0m %s\n' "$*"; }
fail() { printf '\033[1;31m[build_test_bundle] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

require_file() {
    [[ -s "$1" ]] || fail "missing required input: $1"
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || fail "missing required command: $1 (is the pixi env active?)"
}

# ---- Sanity checks ---------------------------------------------------------
log "Checking required inputs..."
require_cmd samtools
require_cmd bcftools
require_cmd seqtk
require_cmd STAR
require_cmd awk
require_cmd bgzip
require_cmd tabix

require_file "${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa"
require_file "${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa.fai"
require_file "${REFERENCE_DIR}/gencode.v44.primary_assembly.annotation.gtf"
require_file "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
require_file "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
require_file "${VARIANT_DIR}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
require_file "${VARIANT_DIR}/small_exac_common_3.hg38.vcf.gz"
require_file "${VARIANT_DIR}/gencode_assembly.bed"
require_file "${FASTQ_R1}"
require_file "${FASTQ_R2}"

mkdir -p "${TEST_BUNDLE_DIR}"/{fastq,reference,variant_calling}
log "Test bundle destination: ${TEST_BUNDLE_DIR}"
log "Subsetting to chromosome: ${CHR}"

# ---- Stage 1: reference FASTA + indexes ------------------------------------
FA_OUT="${TEST_BUNDLE_DIR}/reference/GRCh38.${CHR}.fa"
if [[ ! -s "${FA_OUT}" ]]; then
    log "Extracting ${CHR} from reference FASTA..."
    samtools faidx "${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa" "${CHR}" > "${FA_OUT}"
    samtools faidx "${FA_OUT}"
    samtools dict "${FA_OUT}" -o "${FA_OUT%.fa}.dict"
else
    log "[skip] reference FASTA already built at ${FA_OUT}"
fi

# ---- Stage 2: GTF and RSeQC BED --------------------------------------------
GTF_OUT="${TEST_BUNDLE_DIR}/reference/gencode.v44.${CHR}.gtf"
if [[ ! -s "${GTF_OUT}" ]]; then
    log "Extracting ${CHR} from GTF (grep -w for speed)..."
    # GTF header lines start with '#'; keep them + chr-matching rows
    awk -v c="${CHR}" '$1 == c || /^#/' \
        "${REFERENCE_DIR}/gencode.v44.primary_assembly.annotation.gtf" \
        > "${GTF_OUT}"
else
    log "[skip] GTF already built at ${GTF_OUT}"
fi

BED_OUT="${TEST_BUNDLE_DIR}/reference/gencode_assembly.${CHR}.bed"
if [[ ! -s "${BED_OUT}" ]]; then
    log "Extracting ${CHR} from RSeQC BED..."
    awk -v c="${CHR}" '$1 == c' "${VARIANT_DIR}/gencode_assembly.bed" > "${BED_OUT}"
else
    log "[skip] RSeQC BED already built at ${BED_OUT}"
fi

# ---- Stage 3: STAR index for chr22 -----------------------------------------
STAR_OUT="${TEST_BUNDLE_DIR}/reference/star_index_${CHR}"
if [[ ! -s "${STAR_OUT}/Genome" ]]; then
    log "Building STAR index for ${CHR} (may take 5-10 min)..."
    mkdir -p "${STAR_OUT}"
    # sjdbOverhang = read length - 1 = 150 - 1 = 149 (matches D100 read length)
    # genomeSAindexNbases = min(14, log2(GenomeLength)/2 - 1) = ~11 for chr22 (50 Mb)
    STAR --runMode genomeGenerate \
         --runThreadN "${THREADS}" \
         --genomeDir "${STAR_OUT}" \
         --genomeFastaFiles "${FA_OUT}" \
         --sjdbGTFfile "${GTF_OUT}" \
         --sjdbOverhang 149 \
         --genomeSAindexNbases 11 \
         --outFileNamePrefix "${STAR_OUT}/"
else
    log "[skip] STAR index already built at ${STAR_OUT}"
fi

# ---- Stage 4: variant-calling VCF subsets ----------------------------------
# dbsnp138 is a PLAIN VCF (10.9 GB uncompressed). bcftools accepts plain vcf
# but is faster on bgzipped+indexed input. Use grep for the uncompressed one.
DBSNP_OUT="${TEST_BUNDLE_DIR}/variant_calling/dbsnp138.${CHR}.vcf.gz"
if [[ ! -s "${DBSNP_OUT}.tbi" ]]; then
    log "Subsetting dbsnp138 to ${CHR} (this is the slow one — ~5 min over a 10 GB VCF)..."
    # head + grep is ~2x faster than bcftools on uncompressed input
    {
        awk '/^#/ {print; next} /^[^#]/ {exit}' "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        grep -E "^${CHR}[[:space:]]" "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
    } | bgzip --threads "${THREADS}" > "${DBSNP_OUT}"
    tabix -p vcf "${DBSNP_OUT}"
else
    log "[skip] dbsnp138 already subset at ${DBSNP_OUT}"
fi

subset_vcf_gz() {
    local in_vcf="$1" out_vcf="$2"
    if [[ ! -s "${out_vcf}.tbi" ]]; then
        log "Subsetting $(basename "${in_vcf}") to ${CHR}..."
        bcftools view -r "${CHR}" -Oz -o "${out_vcf}" "${in_vcf}"
        tabix -p vcf "${out_vcf}"
    else
        log "[skip] $(basename "${out_vcf}") already built"
    fi
}

subset_vcf_gz "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz" \
              "${TEST_BUNDLE_DIR}/variant_calling/known_indels.${CHR}.vcf.gz"

subset_vcf_gz "${VARIANT_DIR}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
              "${TEST_BUNDLE_DIR}/variant_calling/Mills.${CHR}.vcf.gz"

subset_vcf_gz "${VARIANT_DIR}/small_exac_common_3.hg38.vcf.gz" \
              "${TEST_BUNDLE_DIR}/variant_calling/small_exac_common_3.${CHR}.vcf.gz"

# ---- Stage 5: subsampled FASTQ ---------------------------------------------
# Strategy: seqtk random subsample to SUBSAMPLE_READS paired reads with a
# fixed seed for reproducibility. At ~100M reads full and 2M subsampled,
# ~40k reads will map to chr22 (chr22 ≈ 1.6% of genome) — enough for STAR
# alignment + Mutect2 to produce a non-empty VCF on the test.
FQ1_OUT="${TEST_BUNDLE_DIR}/fastq/test_R1.fastq.gz"
FQ2_OUT="${TEST_BUNDLE_DIR}/fastq/test_R2.fastq.gz"
if [[ ! -s "${FQ1_OUT}" || ! -s "${FQ2_OUT}" ]]; then
    log "Subsampling FASTQ to ${SUBSAMPLE_READS} read pairs (seed=42)..."
    seqtk sample -s 42 "${FASTQ_R1}" "${SUBSAMPLE_READS}" | bgzip --threads "${THREADS}" > "${FQ1_OUT}"
    seqtk sample -s 42 "${FASTQ_R2}" "${SUBSAMPLE_READS}" | bgzip --threads "${THREADS}" > "${FQ2_OUT}"
else
    log "[skip] FASTQ already subsampled at ${FQ1_OUT} / ${FQ2_OUT}"
fi

# ---- Stage 6: samplesheet --------------------------------------------------
SHEET_OUT="${TEST_BUNDLE_DIR}/samplesheet_test.csv"
cat > "${SHEET_OUT}" <<EOF
sample,fastq_1,fastq_2,strandedness
D100_liver_chr22,${FQ1_OUT},${FQ2_OUT},reverse
EOF
log "Wrote samplesheet: ${SHEET_OUT}"

# ---- Summary ---------------------------------------------------------------
log "Bundle build complete. Contents:"
du -sh "${TEST_BUNDLE_DIR}"/* 2>&1 | awk '{print "    " $0}'
log "Total size:"
du -sh "${TEST_BUNDLE_DIR}" | awk '{print "    " $1}'
log ""
log "Run the pipeline against this bundle with:"
log "    pixi run nextflow run . -profile test,monash,singularity \\"
log "        --outdir /fs04/scratch2/xy86/sanjay/ipg-test-results"
