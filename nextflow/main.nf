// Thin wrapper: all logic lives in the mutation-hunter package. Runs one GWAS per phenotype in parallel,
// then recomputes q over every test. Requires `mutation-hunter` on PATH.
nextflow.enable.dsl = 2

params.phenotypes = 'all'      // 'all' or comma-separated condition names
params.level = 'gene'          // 'gene' or 'variant'
params.max_freq = 0.02
params.min_carriers = 3
params.outdir = 'results'

process FETCH {
    output:
    stdout

    script:
    'mutation-hunter fetch'
}

process GWAS {
    tag "$pheno"
    publishDir params.outdir, mode: 'copy'

    input:
    val pheno

    output:
    path "${pheno}.tsv", emit: tsv
    path "${pheno}.{json,html}"

    script:  // one BLAS thread per task: parallelism comes from running phenotypes side by side
    """
    export OMP_NUM_THREADS=$task.cpus OPENBLAS_NUM_THREADS=$task.cpus MKL_NUM_THREADS=$task.cpus
    mutation-hunter gwas $pheno --level $params.level --max-freq $params.max_freq \\
        --min-carriers $params.min_carriers -o $pheno
    """
}

process FDR {
    publishDir params.outdir, mode: 'copy'

    input:
    path tsvs

    output:
    path 'gwas_all.tsv'

    script:
    "mutation-hunter fdr $tsvs -o gwas_all"
}

workflow {
    wanted = params.phenotypes == 'all' ? null : params.phenotypes.tokenize(',')
    phenos = FETCH().splitText().map { it.trim() }.filter { it && (!wanted || it in wanted) }
    GWAS(phenos)
    FDR(GWAS.out.tsv.collect())
}
