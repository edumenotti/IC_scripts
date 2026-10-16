// Build the panel store from the public 1011 genomes data (maintainers only).
// All logic lives in mutation_hunter.build and bin/annotate.sh. Needs, on PATH: bcftools >= 1.20 and
// mutation-hunter installed with the `build` extra (pip install 'mutation-hunter[build]').
nextflow.enable.dsl = 2

params.inputs = null          // folder with the pinned inputs (mutation_hunter.build.INPUTS); VCF indexed (.tbi)
params.outdir = 'panel'
params.panel_version = null
params.min_gq = 20
params.min_dp = 5

process CHECK {
    input:
    path inputs

    output:
    stdout

    script:
    "python -m mutation_hunter.build check-inputs $inputs"
}

process REFERENCE {
    input:
    path inputs
    val checked

    output:
    path 'reference'

    script:
    "python -m mutation_hunter.build reference $inputs reference"
}

process STRAINS {
    conda 'bioconda::bcftools=1.20'

    input:
    path inputs
    val checked

    output:
    path 'strains.tsv', emit: strains
    path 'sites.tsv', emit: sites

    script:
    """
    bcftools query -l $inputs/1011Matrix.gvcf.gz > samples.txt
    python -m mutation_hunter.build strains $inputs samples.txt strains.tsv
    python -m mutation_hunter.build sites $inputs sites.tsv
    """
}

process ANNOTATE {
    tag "$chrom"
    conda 'bioconda::bcftools=1.20'

    input:
    val chrom
    path inputs
    path reference

    output:
    tuple val(chrom), path("${chrom}.bcf"), path("${chrom}.bcf.csi")

    script:
    "annotate.sh $inputs/1011Matrix.gvcf.gz $chrom $reference/ref.fa $reference/ref.gff3 ${chrom}.bcf"
}

process GENOTYPES {
    tag "$chrom"
    memory '4 GB'

    input:
    tuple val(chrom), path(bcf), path(csi)
    path strains
    path reference
    path sites

    output:
    path "${chrom}.*.{tsv.gz,npz}"

    script:
    """
    python -m mutation_hunter.build genotypes $bcf $strains $reference/ref.gff3 $chrom \\
        --min-gq $params.min_gq --min-dp $params.min_dp --sites $sites
    """
}

process ASSEMBLE {
    publishDir params.outdir, mode: 'copy'
    conda 'bioconda::bcftools=1.20'  // same bcftools as ANNOTATE, so the manifest records the one that ran
    memory '16 GB'

    input:
    path parts
    path inputs
    path reference
    path strains
    val provenance

    output:
    path "panel-${params.panel_version}"

    script:
    """
    mkdir parts && mv ${parts} parts/
    python -m mutation_hunter.build assemble parts $inputs $reference $strains panel-${params.panel_version} \\
        --panel-version ${params.panel_version} \\
        --params '{"min_gq": $params.min_gq, "min_dp": $params.min_dp, "bcftools": "'"\$(bcftools --version | head -1)"'", ${provenance}}'
    """
}

workflow {
    if (!params.inputs || !params.panel_version) error "set --inputs <folder> and --panel_version <x.y.z>"
    inputs = file(params.inputs, checkIfExists: true)
    commit = "git -C $projectDir describe --always --dirty --abbrev=40".execute().text.trim()
    provenance = "\"nextflow\": \"${nextflow.version}\", \"pipeline_commit\": \"${commit}\""
    checked = CHECK(inputs)
    reference = REFERENCE(inputs, checked)
    STRAINS(inputs, checked)
    strains = STRAINS.out.strains
    chroms = Channel.of(1..16).map { "chromosome$it" }
    ANNOTATE(chroms, inputs, reference)
    GENOTYPES(ANNOTATE.out, strains, reference, STRAINS.out.sites)
    ASSEMBLE(GENOTYPES.out.flatten().collect(), inputs, reference, strains, provenance)
}
