#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.inputDir = "$baseDir/flat_files/" // Diretório padrão para os arquivos .flat é o diretório atual
params.scriptPath = "$baseDir/genbank_extract.py" // O script genbank.py está no mesmo diretório que o script .nf
params.blastDb = "$baseDir/blastdb/PE2" // Caminho para o banco de dados do blastp
params.outputDir = "$baseDir/mutations_results/"
params.analysisScript = "$baseDir/ult_pc_2.py"

workflow {

    flatFiles = Channel.fromFilePairs("${params.inputDir}/*.flat", size: 1, flat: true)

    extractSequences = flatFiles.map { name, flat ->
        tuple(name, file(flat), params.scriptPath)
    }

    ExtractSequences(extractSequences)
        .flatMap { name, fastaFiles -> fastaFiles.collect { file -> tuple(name, file) } }
        .set { fastaFiles }
    
    RunBlastp(fastaFiles)
        .map { name, dir -> tuple(name, dir) }
        .set { xmlDirs }

    AnalyzeMutations(xmlDirs)
}

process ExtractSequences {
    publishDir "$baseDir/flat_files", mode: 'symlink'

    input:
    tuple val(name), path(flat), val(scriptPath)
    
    output:
    tuple val(name), path("${name}_fasta/*.fasta")
    
    script:
    """
    mkdir -p ${name}_fasta
    python ${scriptPath} ${flat}
    mv *.fasta ${name}_fasta/
    """
}

process RunBlastp {
    publishDir "$baseDir/blast_results", mode: 'symlink'

    input:
    tuple val(name), path(fastaFile)
    
    output:
    tuple val(name), path("${name}_xml"), emit: xmlDir

    errorStrategy 'ignore'
        
    script:
    """
    baseName=\$(basename "${fastaFile}" .fasta)
    outputDir="${name}_xml"
    mkdir -p "\$outputDir"
    outputFile="\$outputDir/\$baseName.xml"

    if [ -s "${fastaFile}" ]; then
        blastp -query "${fastaFile}" -db ${params.blastDb} -out "\$outputFile" -outfmt 5
        if [ ! -s "\$outputFile" ]; then
            rm "\$outputFile"
        fi
    else
        echo "Warning: {fastaFile} is empty or does not exist, skipping..." >&2
    fi
    """
}

process AnalyzeMutations {

    input:
    tuple val(name), path(xmlDir)

    errorStrategy 'ignore'

    script:
    """
    echo "Caminho xmlDir: ${xmlDir}"
    outputBaseDir=\$(basename ${xmlDir} _xml)  # Correção na atribuição
    outputDir="${params.outputDir}/\${outputBaseDir}_results"  # Uso correto das variáveis do Nextflow e shell
    mkdir -p "\$outputDir"
    for xmlFile in \$(ls ${xmlDir}/*.xml); do
        xmlBaseName=\$(basename "\$xmlFile" .xml)
        if ! python ${params.analysisScript} -f "\$xmlFile" -e -od \$outputDir; then
            echo "Error: Failed to analyze mutations for \${xmlBaseName}" >&2
        fi
    done
    """
}
