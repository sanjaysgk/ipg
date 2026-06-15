process MSCONVERT {
    tag "${meta.id}"
    label 'process_low'

    // ProteoWizard msconvert with vendor readers (Bruker timsTOF .d, etc.). The
    // vendor-reader build is NOT on biocontainers (vendor licence) — this is the
    // ProteoWizard "i-agree-to-the-vendor-licenses" image, run via Wine on Linux.
    // Heavy (~GB). Used only for Bruker .d; Thermo .raw goes through THERMORAWFILEPARSER.
    // Container + invocation pattern follow lehtiolab/nf-msconvert (a proven wrapper).
    container "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:3.0.25182-6ca4f8b"

    // Singularity is read-only; msconvert's Wine needs a writable prefix.
    containerOptions { workflow.containerEngine == 'singularity' ? '-S /mywineprefix' : '' }

    input:
    tuple val(meta), path(spectra)            // a Bruker .d directory (or other vendor format)

    output:
    tuple val(meta), path("*.mzML"), emit: spectra
    path "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Under singularity the image exposes `mywine` (sets WINEPREFIX); plain `wine` elsewhere.
    def winecmd = workflow.containerEngine == 'singularity' ? 'mywine' : 'wine'
    def args    = task.ext.args ?: '--zlib --filter "peakPicking vendor msLevel=1-" --filter "combineIonMobilitySpectra"'
    def prefix  = task.ext.prefix ?: "${meta.id}"
    // NOTE: ion-mobility filters + exact flags should be validated against real
    // Bruker .d data — we have none yet. `combineIonMobilitySpectra` collapses
    // timsTOF ion mobility for standard DDA search.
    """
    # Nextflow stages a .d as a symlinked directory; the Bruker reader needs real
    # files, so dereference into place (pattern from lehtiolab/nf-msconvert).
    if [ -d "${spectra}" ]; then
        mv ${spectra} tmpdir && cp -rL tmpdir ${spectra} && rm -rf tmpdir
    fi

    ${winecmd} msconvert ${spectra} \\
        --mzML --outfile ${prefix}.mzML -o . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msconvert: "3.0.25182-6ca4f8b"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mzML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        msconvert: stub
    END_VERSIONS
    """
}
