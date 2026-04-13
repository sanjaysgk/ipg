#!/usr/bin/env bash
#
# build_full_reference.sh — download and prepare the complete GRCh38
# reference bundle that sanjaysgk/ipg needs for a real (non-test) run.
#
# What it produces:
#
#     $REFERENCE_BASE/
#     ├── genome/
#     │   ├── GRCh38.primary_assembly.genome.fa
#     │   ├── GRCh38.primary_assembly.genome.fa.fai
#     │   ├── GRCh38.primary_assembly.genome.dict
#     │   ├── gencode.v44.primary_assembly.annotation.gtf
#     │   ├── gencode_assembly.bed         (RSeQC format)
#     │   └── star_index/
#     │       ├── Genome
#     │       ├── SAindex
#     │       └── ...
#     └── gatk_resources/
#         ├── dbsnp138.vcf.gz{,.tbi}
#         ├── known_indels.vcf.gz{,.tbi}
#         ├── Mills.vcf.gz{,.tbi}
#         └── small_exac_common_3.vcf.gz{,.tbi}
#
# ---- Public download sources ----
#
# All files are fetched from official repositories:
#
#   Genome + GTF:
#     GENCODE release 44 (GRCh38 primary assembly)
#     https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
#
#   GATK variant calling resources (Broad hg38 bundle):
#     https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0
#     (all files below accessed via direct HTTPS from storage.googleapis.com)
#
#   RSeQC BED:
#     Derived from GENCODE GTF using gtfToGenePred + genePredToBed (UCSC tools)
#     OR from SourceForge: https://sourceforge.net/projects/rseqc/files/BED/
#
# ---- Usage ----
#
# From the repo root, inside the pixi env:
#     pixi run bash bin/build_full_reference.sh
#
# Or override locations with env vars:
#     REFERENCE_BASE=/my/ref/dir THREADS=16 \
#         pixi run bash bin/build_full_reference.sh
#
# To reuse existing local files (e.g. on Monash M3 where they already
# exist), set GENOME_DIR and VARIANT_DIR to skip downloading:
#     GENOME_DIR=/fs04/.../GRCh38 \
#     VARIANT_DIR=/fs04/.../variant_calling_resources \
#         pixi run bash bin/build_full_reference.sh
#
# IMPORTANT: STAR genomeGenerate needs ~32 GB RAM. Do NOT run this on a
# login node — submit it via SLURM or use an interactive compute node:
#     smux new-session --time=2:00:00 --ntasks=1 --cpus-per-task=16 --mem=48G
#
# Re-running is idempotent: completed stages are skipped.

set -euo pipefail

# ---- Configuration ---------------------------------------------------------
REFERENCE_BASE="${REFERENCE_BASE:-$(pwd)/references/GRCh38}"
THREADS="${THREADS:-16}"
SJDB_OVERHANG="${SJDB_OVERHANG:-149}"  # read length - 1 (150bp ATLANTIS RNA-seq)

# Optional: point at existing local files to skip download
GENOME_DIR="${GENOME_DIR:-}"
VARIANT_DIR="${VARIANT_DIR:-}"

# Output directories
GENOME_OUT="${REFERENCE_BASE}/genome"
GATK_OUT="${REFERENCE_BASE}/gatk_resources"

# ---- Public URLs -----------------------------------------------------------
# GENCODE release 44 (GRCh38 primary assembly)
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44"
FASTA_URL="${GENCODE_BASE}/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="${GENCODE_BASE}/gencode.v44.primary_assembly.annotation.gtf.gz"

# Broad GATK hg38 resource bundle (Google Cloud public bucket)
BROAD_BASE="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0"
DBSNP_URL="${BROAD_BASE}/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
DBSNP_TBI_URL="${BROAD_BASE}/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi"
KNOWN_INDELS_URL="${BROAD_BASE}/Homo_sapiens_assembly38.known_indels.vcf.gz"
KNOWN_INDELS_TBI_URL="${BROAD_BASE}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
MILLS_URL="${BROAD_BASE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
MILLS_TBI_URL="${BROAD_BASE}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"

# gnomAD germline resource for Mutect2 (small exac common sites)
GNOMAD_URL="${BROAD_BASE}/small_exac_common_3.hg38.vcf.gz"
GNOMAD_TBI_URL="${BROAD_BASE}/small_exac_common_3.hg38.vcf.gz.tbi"

# RSeQC BED (hg38, SourceForge mirror of standard RSeQC BED)
RSEQC_BED_URL="https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_GENCODE_V44_Basic.bed.gz/download"

# ---- Helpers ---------------------------------------------------------------
log()  { printf '\033[1;34m[build_full_reference]\033[0m %s\n' "$*"; }
warn() { printf '\033[1;33m[build_full_reference] WARN:\033[0m %s\n' "$*"; }
fail() { printf '\033[1;31m[build_full_reference] ERROR:\033[0m %s\n' "$*" >&2; exit 1; }

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || fail "missing required command: $1 (is the pixi env active?)"
}

# Download a file if the target doesn't exist. Resumes partial downloads.
download() {
    local url="$1" dest="$2"
    if [[ -s "${dest}" ]]; then
        log "[skip] already exists: ${dest}"
        return 0
    fi
    log "Downloading: ${url}"
    log "         to: ${dest}"
    mkdir -p "$(dirname "${dest}")"
    curl -fSL --retry 3 -C - -o "${dest}.tmp" "${url}"
    mv "${dest}.tmp" "${dest}"
}

# Copy a local file if available, otherwise download from URL
fetch() {
    local local_path="$1" url="$2" dest="$3"
    if [[ -n "${local_path}" && -s "${local_path}" ]]; then
        if [[ -s "${dest}" ]]; then
            log "[skip] already exists: ${dest}"
        else
            log "Copying local: ${local_path} -> ${dest}"
            mkdir -p "$(dirname "${dest}")"
            cp "${local_path}" "${dest}"
        fi
    else
        download "${url}" "${dest}"
    fi
}

# ---- Sanity checks ---------------------------------------------------------
log "Checking required tools..."
require_cmd curl
require_cmd samtools
require_cmd bgzip
require_cmd tabix
require_cmd STAR

mkdir -p "${GENOME_OUT}" "${GATK_OUT}"
log "Reference base: ${REFERENCE_BASE}"
log "Threads: ${THREADS}"
if [[ -n "${GENOME_DIR}" ]]; then
    log "Using local genome dir: ${GENOME_DIR}"
fi
if [[ -n "${VARIANT_DIR}" ]]; then
    log "Using local variant dir: ${VARIANT_DIR}"
fi

# ---- Stage 1: Genome FASTA ------------------------------------------------
FA="${GENOME_OUT}/GRCh38.primary_assembly.genome.fa"
if [[ -s "${FA}" ]]; then
    log "[skip] Genome FASTA already at ${FA}"
else
    if [[ -n "${GENOME_DIR}" && -s "${GENOME_DIR}/GRCh38.primary_assembly.genome.fa" ]]; then
        log "Copying local genome FASTA..."
        cp "${GENOME_DIR}/GRCh38.primary_assembly.genome.fa" "${FA}"
    else
        download "${FASTA_URL}" "${FA}.gz"
        log "Decompressing genome FASTA..."
        gunzip "${FA}.gz"
    fi
fi

# ---- Stage 2: FASTA index + dict ------------------------------------------
if [[ -s "${FA}.fai" ]]; then
    log "[skip] FASTA index already at ${FA}.fai"
else
    if [[ -n "${GENOME_DIR}" && -s "${GENOME_DIR}/GRCh38.primary_assembly.genome.fa.fai" ]]; then
        cp "${GENOME_DIR}/GRCh38.primary_assembly.genome.fa.fai" "${FA}.fai"
    else
        log "Building FASTA index..."
        samtools faidx "${FA}"
    fi
fi

DICT="${FA%.fa}.dict"
if [[ -s "${DICT}" ]]; then
    log "[skip] Sequence dictionary already at ${DICT}"
else
    if [[ -n "${GENOME_DIR}" && -s "${GENOME_DIR}/GRCh38.primary_assembly.genome.dict" ]]; then
        cp "${GENOME_DIR}/GRCh38.primary_assembly.genome.dict" "${DICT}"
    else
        log "Building sequence dictionary..."
        samtools dict "${FA}" -o "${DICT}"
    fi
fi

# ---- Stage 3: GTF ----------------------------------------------------------
GTF="${GENOME_OUT}/gencode.v44.primary_assembly.annotation.gtf"
if [[ -s "${GTF}" ]]; then
    log "[skip] GTF already at ${GTF}"
else
    if [[ -n "${GENOME_DIR}" && -s "${GENOME_DIR}/gencode.v44.primary_assembly.annotation.gtf" ]]; then
        log "Copying local GTF..."
        cp "${GENOME_DIR}/gencode.v44.primary_assembly.annotation.gtf" "${GTF}"
    else
        download "${GTF_URL}" "${GTF}.gz"
        log "Decompressing GTF..."
        gunzip "${GTF}.gz"
    fi
fi

# ---- Stage 4: RSeQC BED ---------------------------------------------------
BED="${GENOME_OUT}/gencode_assembly.bed"
if [[ -s "${BED}" ]]; then
    log "[skip] RSeQC BED already at ${BED}"
else
    if [[ -n "${VARIANT_DIR}" && -s "${VARIANT_DIR}/gencode_assembly.bed" ]]; then
        cp "${VARIANT_DIR}/gencode_assembly.bed" "${BED}"
    elif [[ -n "${GENOME_DIR}" && -s "${GENOME_DIR}/gencode_assembly.bed" ]]; then
        cp "${GENOME_DIR}/gencode_assembly.bed" "${BED}"
    else
        log "Downloading RSeQC BED from SourceForge..."
        curl -fSL --retry 3 -C - -o "${BED}.gz" -L "${RSEQC_BED_URL}"
        gunzip "${BED}.gz"
    fi
fi

# ---- Stage 5: STAR index ---------------------------------------------------
STAR_DIR="${GENOME_OUT}/star_index"
if [[ -s "${STAR_DIR}/Genome" && -s "${STAR_DIR}/SAindex" ]]; then
    log "[skip] STAR index already built at ${STAR_DIR}"
else
    if [[ -n "${GENOME_DIR}" && -d "${GENOME_DIR}/star_index" && -s "${GENOME_DIR}/star_index/Genome" ]]; then
        log "Copying local STAR index..."
        mkdir -p "${STAR_DIR}"
        cp "${GENOME_DIR}/star_index"/* "${STAR_DIR}/"
    else
        log "Building full-genome STAR index (needs ~32 GB RAM, ~30-60 min)..."
        log "  genomeDir   : ${STAR_DIR}"
        log "  sjdbOverhang: ${SJDB_OVERHANG}"
        log "  threads     : ${THREADS}"
        mkdir -p "${STAR_DIR}"
        STAR --runMode genomeGenerate \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_DIR}" \
            --genomeFastaFiles "${FA}" \
            --sjdbGTFfile "${GTF}" \
            --sjdbOverhang "${SJDB_OVERHANG}" \
            --outFileNamePrefix "${STAR_DIR}/"
        log "STAR index built: $(du -sh "${STAR_DIR}" | awk '{print $1}')"
    fi
fi

# ---- Stage 6: GATK variant calling resources --------------------------------
# dbSNP138
DBSNP="${GATK_OUT}/dbsnp138.vcf.gz"
if [[ -s "${DBSNP}" && -s "${DBSNP}.tbi" ]]; then
    log "[skip] dbSNP138 already at ${DBSNP}"
else
    if [[ -n "${VARIANT_DIR}" && -s "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf" ]]; then
        # Local Monash copy is plain VCF — bgzip + tabix it
        log "Bgzipping local dbSNP138 (11 GB plain VCF, ~5 min)..."
        bgzip -@ "${THREADS}" -c \
            "${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf" \
            > "${DBSNP}.tmp"
        mv "${DBSNP}.tmp" "${DBSNP}"
        tabix -p vcf "${DBSNP}"
    else
        download "${DBSNP_URL}" "${DBSNP}"
        download "${DBSNP_TBI_URL}" "${DBSNP}.tbi"
    fi
fi

# Known indels
KI="${GATK_OUT}/known_indels.vcf.gz"
if [[ -s "${KI}" && -s "${KI}.tbi" ]]; then
    log "[skip] known_indels already at ${KI}"
else
    fetch \
        "${VARIANT_DIR:+${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz}" \
        "${KNOWN_INDELS_URL}" "${KI}"
    if [[ ! -s "${KI}.tbi" ]]; then
        fetch \
            "${VARIANT_DIR:+${VARIANT_DIR}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi}" \
            "${KNOWN_INDELS_TBI_URL}" "${KI}.tbi"
    fi
fi

# Mills + 1000G gold standard indels
MILLS="${GATK_OUT}/Mills.vcf.gz"
if [[ -s "${MILLS}" && -s "${MILLS}.tbi" ]]; then
    log "[skip] Mills already at ${MILLS}"
else
    fetch \
        "${VARIANT_DIR:+${VARIANT_DIR}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz}" \
        "${MILLS_URL}" "${MILLS}"
    if [[ ! -s "${MILLS}.tbi" ]]; then
        fetch \
            "${VARIANT_DIR:+${VARIANT_DIR}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi}" \
            "${MILLS_TBI_URL}" "${MILLS}.tbi"
    fi
fi

# gnomAD / ExAC germline resource (for Mutect2 tumour-only calling)
GNOMAD="${GATK_OUT}/small_exac_common_3.vcf.gz"
if [[ -s "${GNOMAD}" && -s "${GNOMAD}.tbi" ]]; then
    log "[skip] germline resource already at ${GNOMAD}"
else
    fetch \
        "${VARIANT_DIR:+${VARIANT_DIR}/small_exac_common_3.hg38.vcf.gz}" \
        "${GNOMAD_URL}" "${GNOMAD}"
    if [[ ! -s "${GNOMAD}.tbi" ]]; then
        fetch \
            "${VARIANT_DIR:+${VARIANT_DIR}/small_exac_common_3.hg38.vcf.gz.tbi}" \
            "${GNOMAD_TBI_URL}" "${GNOMAD}.tbi"
    fi
fi

# ---- Stage 7: Generate params YAML -----------------------------------------
PARAMS="${REFERENCE_BASE}/params_reference.yaml"
cat > "${PARAMS}" <<YAML
# Auto-generated by build_full_reference.sh on $(date -u +%Y-%m-%dT%H:%M:%SZ)
# Source: GENCODE v44 + Broad GATK hg38 bundle
#
# Usage:
#   nextflow run sanjaysgk/ipg -params-file ${PARAMS} \\
#       --input samplesheet.csv --outdir results

# Genome
fasta:                 ${FA}
fasta_fai:             ${FA}.fai
fasta_dict:            ${DICT}
gtf:                   ${GTF}
star_index:            ${STAR_DIR}
rseqc_bed:             ${BED}

# GATK variant calling resources
dbsnp:                 ${DBSNP}
dbsnp_tbi:             ${DBSNP}.tbi
known_indels:          ${KI}
known_indels_tbi:      ${KI}.tbi
mills:                 ${MILLS}
mills_tbi:             ${MILLS}.tbi
germline_resource:     ${GNOMAD}
germline_resource_tbi: ${GNOMAD}.tbi
YAML

# ---- Summary ---------------------------------------------------------------
log ""
log "========================================="
log "  Reference bundle ready"
log "========================================="
log ""
log "Location: ${REFERENCE_BASE}"
du -sh "${GENOME_OUT}" "${GATK_OUT}" 2>&1 | awk '{print "    " $0}'
log ""
log "Auto-generated params file: ${PARAMS}"
log ""
log "Run the pipeline:"
log ""
log "  nextflow run sanjaysgk/ipg \\"
log "      -params-file ${PARAMS} \\"
log "      --input samplesheet.csv \\"
log "      --outdir results"
log ""
log "Public sources used:"
log "  Genome + GTF : GENCODE v44  (ftp.ebi.ac.uk)"
log "  GATK bundle  : Broad hg38   (storage.googleapis.com)"
log "  RSeQC BED    : SourceForge   (rseqc project)"
