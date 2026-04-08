# ipg-tools: reproducible build of the custom C tools from the
# immunopeptidogenomics repository (Scull et al. 2020).
#
# Published as ghcr.io/sanjaysgk/ipg-tools:<tag> via
# .github/workflows/build-ipg-tools.yml
#
# Pin both the source repo and commit SHA so that rebuilding this image
# at any point in the future produces byte-identical binaries.

FROM debian:bookworm-slim AS builder

# ---- Pinned source ---------------------------------------------------------
# Fork of kescull/immunopeptidogenomics under sanjaysgk control, pinned to
# commit fef8e68 ("Fix bug when squishing more than 9 fasta files").
# To bump: update IPG_REPO_SHA, rebuild, verify tests.
ARG IPG_REPO_URL=https://github.com/sanjaysgk/immunopeptidogenomics.git
ARG IPG_REPO_SHA=fef8e68ab86dac49070542fe40b477ec76058aaa

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        git \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /src
RUN git clone --no-checkout "${IPG_REPO_URL}" ipg \
 && cd ipg \
 && git checkout "${IPG_REPO_SHA}" \
 && git log -1 --format='%H %s' > /src/ipg.build-info

WORKDIR /src/ipg
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
    strip /out/*

# ---- Runtime image ---------------------------------------------------------
FROM debian:bookworm-slim

# Runtime deps: the binaries are dynamically linked against glibc. On Debian
# bookworm libm and libpthread are folded into libc6 (glibc 2.36+), so libc6
# alone is sufficient. ca-certificates is not needed at runtime because the
# git clone happens in the builder stage.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        libc6 \
 && rm -rf /var/lib/apt/lists/*

COPY --from=builder /out/curate_vcf       /usr/local/bin/curate_vcf
COPY --from=builder /out/revert_headers   /usr/local/bin/revert_headers
COPY --from=builder /out/alt_liftover     /usr/local/bin/alt_liftover
COPY --from=builder /out/triple_translate /usr/local/bin/triple_translate
COPY --from=builder /out/squish           /usr/local/bin/squish
COPY --from=builder /src/ipg.build-info   /etc/ipg-tools.build-info

LABEL org.opencontainers.image.source="https://github.com/sanjaysgk/ipg"
LABEL org.opencontainers.image.description="Custom C tools from immunopeptidogenomics (curate_vcf, revert_headers, alt_liftover, triple_translate, squish) used by the sanjaysgk/ipg nf-core pipeline."
LABEL org.opencontainers.image.licenses="MIT"
LABEL org.opencontainers.image.documentation="https://github.com/sanjaysgk/ipg/blob/main/docs/ipg-tools.md"

# Sanity check: each binary must be present, executable, and linkable.
# We do NOT run `<tool> -h` here because curate_vcf.c and revert_headers.c
# do not implement -h (curate_vcf getopt string is "d" only; revert_headers
# uses positional args). Invoking them with -h would exit nonzero and mask
# real build problems. Instead we verify existence, executability, and that
# the dynamic linker can resolve every shared library dependency.
RUN set -eu; \
    for t in curate_vcf revert_headers alt_liftover triple_translate squish; do \
        echo "=== $t ==="; \
        test -x "/usr/local/bin/$t" || { echo "missing or non-executable: $t"; exit 1; }; \
        ldd "/usr/local/bin/$t" | grep -E 'not found|undefined' && exit 1 || true; \
        ldd "/usr/local/bin/$t"; \
    done

CMD ["bash"]
