#!/usr/bin/env bash
#
# build_test_bundle.sh — MAINTAINER-ONLY script. Builds the chr22 test bundle
# that ships as the v0.1.0-test-data GitHub Release asset.
#
# End users should NOT run this. They run `scripts/fetch_test_bundle.sh` instead,
# which downloads the prebuilt tarball.
#
# Required env vars (no defaults; this script fails fast if missing):
#   REFERENCE_DIR   dir containing GRCh38.primary_assembly.genome.fa + gencode.v44 GTF
#   VARIANT_DIR     dir containing GATK bundle VCFs + gencode_assembly.bed
#   FASTQ_R1        R1 FASTQ.gz from a real human RNAseq sample
#   FASTQ_R2        matching R2
#   TEST_BUNDLE_DIR output bundle dir
#
# Optional:
#   CHR=chr22        target chromosome
#   SUBSAMPLE_READS=200000  paired reads to keep (seqtk -s 42)
#   THREADS=4
#
# Public download hints (in case you don't have these on hand):
#   GENCODE GRCh38 + v44 GTF: https://www.gencodegenes.org/human/release_44.html
#   GATK bundle:              gs://gcp-public-data--broad-references/hg38/v0/
#
# Output: $TEST_BUNDLE_DIR/{fastq,reference,variant_calling,samplesheet_test.csv}
#         + $(dirname $TEST_BUNDLE_DIR)/ipg-test-bundle-chr22.tar.gz
#
# Re-runs are idempotent; delete a sub-dir to force re-build of that stage.

set -euo pipefail

# ---- Configuration ---------------------------------------------------------
# All paths are REQUIRED env vars. No Monash defaults — this script must
# fail loudly on someone else's HPC, not silently use Sanjay's paths.
: "${REFERENCE_DIR:?REFERENCE_DIR must be set — directory containing GRCh38.primary_assembly.genome.fa + gencode.v44.primary_assembly.annotation.gtf}"
: "${VARIANT_DIR:?VARIANT_DIR must be set — directory containing GATK bundle VCFs (dbsnp138, known_indels, Mills, small_exac_common_3) + gencode_assembly.bed}"
: "${FASTQ_R1:?FASTQ_R1 must be set — path to a real R1 FASTQ.gz from a human RNAseq sample}"
: "${FASTQ_R2:?FASTQ_R2 must be set — path to matching R2 FASTQ.gz}"
: "${TEST_BUNDLE_DIR:?TEST_BUNDLE_DIR must be set — output directory for the bundle}"

CHR="${CHR:-chr22}"
SUBSAMPLE_READS="${SUBSAMPLE_READS:-200000}"
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

subset_vcf_gz \
    "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz" \
    "${TEST_BUNDLE_DIR}/variant_calling/known_indels.${CHR}.vcf.gz"

subset_vcf_gz \
    "${VARIANT_DIR}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
    "${TEST_BUNDLE_DIR}/variant_calling/Mills.${CHR}.vcf.gz"

subset_vcf_gz \
    "${VARIANT_DIR}/small_exac_common_3.hg38.vcf.gz" \
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
# nf-validation requires ABSOLUTE paths in the samplesheet — relative paths
# fail with "the file or directory '<rel>' does not exist". The bundle is
# portable across machines, so we can't hardcode the absolute path here.
# Use __BUNDLE_DIR__ placeholder; scripts/fetch_test_bundle.sh substitutes the
# actual extraction dir on first run.
SHEET_OUT="${TEST_BUNDLE_DIR}/samplesheet_test.csv"
cat > "${SHEET_OUT}" <<EOF
sample,fastq_1,fastq_2,strandedness
D100_liver_chr22,__BUNDLE_DIR__/fastq/test_R1.fastq.gz,__BUNDLE_DIR__/fastq/test_R2.fastq.gz,reverse
EOF
log "Wrote samplesheet (with __BUNDLE_DIR__ placeholder): ${SHEET_OUT}"

# ---- Stage 7: tarball for GitHub Release upload ---------------------------
TARBALL="$(dirname "${TEST_BUNDLE_DIR}")/ipg-test-bundle-chr22.tar.gz"
if [[ ! -s "${TARBALL}" || "${TEST_BUNDLE_DIR}" -nt "${TARBALL}" ]]; then
    log "Creating release tarball: ${TARBALL}"
    tar -czf "${TARBALL}" -C "$(dirname "${TEST_BUNDLE_DIR}")" "$(basename "${TEST_BUNDLE_DIR}")"
    log "Tarball size: $(du -h "${TARBALL}" | cut -f1)"
else
    log "[skip] Tarball already up to date at ${TARBALL}"
fi

# ---- Summary ---------------------------------------------------------------
log "Bundle build complete. Contents:"
du -sh "${TEST_BUNDLE_DIR}"/* 2>&1 | awk '{print "    " $0}'
log "Total size:"
du -sh "${TEST_BUNDLE_DIR}" | awk '{print "    " $1}'
log ""
log "Upload to GitHub Release with:"
log "    gh release create v0.1.0-test-data ${TARBALL} \\"
log "        --title 'Test data bundle v0.1.0 (chr22)' \\"
log "        --notes 'chr22-subset reference + 200k-pair FASTQ for -profile test'"
