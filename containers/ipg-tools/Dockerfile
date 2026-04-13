# ipg-tools: reproducible build of the custom C tools used by the
# sanjaysgk/ipg cryptic peptide pipeline.
#
# The five .c sources live alongside this Dockerfile under ./src/. Three of
# them are unmodified copies from kescull/immunopeptidogenomics; two are
# locally improved versions maintained by the Li/Purcell labs:
#
#   - curate_vcf.c       — locally improved (V2): larger field buffers
#                          (chrom 200, ref/alt 4096, filter 1024), dynamic
#                          INFO field allocation, explicit memory cleanup.
#                          Originally curate_vcfV2.c in IPG_slurm/.
#   - revert_headers.c   — locally improved: optional 3rd positional
#                          argument for the output FASTA name (defaults to
#                          tmpc.fasta when omitted, for backward
#                          compatibility with the upstream 2-arg form).
#                          Originally revert_headers.c in IPG_slurm/.
#   - alt_liftover.c     — unmodified, kescull/immunopeptidogenomics.
#   - triple_translate.c — unmodified, kescull/immunopeptidogenomics.
#   - squish.c           — unmodified, kescull/immunopeptidogenomics.
#
# Bundling the sources directly in the repository (rather than git-cloning
# from kescull or sanjaysgk/immunopeptidogenomics at build time) makes the
# image build hermetic and the source-of-truth for every binary lives in
# the same git history as the Nextflow modules that consume them.
#
# Published as ghcr.io/sanjaysgk/ipg-tools:<tag> via
# .github/workflows/build-ipg-tools.yml on every push to main and on every
# git tag matching ipg-tools-v*.

FROM debian:bookworm-slim AS builder

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /src
COPY src/ ./

# Linker flag order matters under ld --as-needed (default on Debian bookworm):
# libraries MUST appear after the translation units that reference them, or
# the linker may drop the library resolution. squish.c uses pthread_create
# and pow() so requires -pthread -lm; the others are plain libc.
RUN set -eux; \
    mkdir -p /out; \
    gcc -O2 -Wall -o /out/curate_vcf       curate_vcf.c; \
    gcc -O2 -Wall -o /out/revert_headers   revert_headers.c; \
    gcc -O2 -Wall -o /out/alt_liftover     alt_liftover.c; \
    gcc -O2 -Wall -o /out/triple_translate triple_translate.c; \
    gcc -O2 -Wall -o /out/squish           squish.c -pthread -lm; \
    gcc -O2 -Wall -o /out/origins         origins.c cJSON.c -lcurl -lm; \
    strip /out/*; \
    echo "ipg-tools built from local sources at $(date -u +%Y-%m-%dT%H:%M:%SZ)" \
        > /out/ipg.build-info

# ---- Runtime image ---------------------------------------------------------
FROM debian:bookworm-slim

# Runtime deps: the binaries are dynamically linked against glibc. On Debian
# bookworm libm and libpthread are folded into libc6 (glibc 2.36+), so libc6
# alone is sufficient.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libc6 \
        libcurl4 && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /out/curate_vcf       /usr/local/bin/curate_vcf
COPY --from=builder /out/revert_headers   /usr/local/bin/revert_headers
COPY --from=builder /out/alt_liftover     /usr/local/bin/alt_liftover
COPY --from=builder /out/triple_translate /usr/local/bin/triple_translate
COPY --from=builder /out/squish           /usr/local/bin/squish
COPY --from=builder /out/origins          /usr/local/bin/origins
COPY --from=builder /out/ipg.build-info   /etc/ipg-tools.build-info

LABEL org.opencontainers.image.source="https://github.com/sanjaysgk/ipg"
LABEL org.opencontainers.image.description="Custom C tools (curate_vcf [V2], revert_headers [3-arg], alt_liftover, triple_translate, squish, origins) used by the sanjaysgk/ipg cryptic peptide pipeline. Built from sources bundled in the pipeline repo at containers/ipg-tools/src/."
LABEL org.opencontainers.image.licenses="MIT"

# Sanity check: each binary must be present, executable, and linkable.
# We do NOT run `<tool> -h` here because curate_vcf.c and revert_headers.c
# do not implement -h (curate_vcf getopt string is "d" only; revert_headers
# uses positional args). Invoking them with -h would exit nonzero and mask
# real build problems. Instead we verify existence, executability, and that
# the dynamic linker can resolve every shared library dependency.
RUN set -eu; \
    for t in curate_vcf revert_headers alt_liftover triple_translate squish origins; do \
        echo "=== $t ==="; \
        test -x "/usr/local/bin/$t" || { echo "missing or non-executable: $t"; exit 1; }; \
        ldd "/usr/local/bin/$t" | grep -E 'not found|undefined' && exit 1 || true; \
        ldd "/usr/local/bin/$t"; \
    done

CMD ["bash"]
