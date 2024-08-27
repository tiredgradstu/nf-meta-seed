// QC and trim reads
process FASTP {
    container 'oras://community.wave.seqera.io/library/fastp:0.23.4--4ea6310369653ec7'
    publishDir "$params.outdir/fastp", mode: 'copy', pattern: '*html'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*fastq.gz'), emit: fastq
    path '*fastp_report*', emit: report

    script:
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample_id}_R1_trimmed.fastq.gz \
        -O ${sample_id}_R2_trimmed.fastq.gz \
        -h ${sample_id}_fastp_report.html \
        -j ${sample_id}_fastp_report.json \
        -w $task.cpus
    """
}

// Download and extract a pre-made Kraken DB from a URL
process KRAKEN_DL_DB {
    input:
    val db_url

    output:
    path 'kraken_db'

    script:
    """
    wget $db_url
    tar -xzvf *.tar.gz -C kraken_db
    """
}

// Download the NCBI taxonomy as a starting point for a custom Kraken DB
process KRAKEN_DL_TAX {
    container 'oras://community.wave.seqera.io/library/kraken2:2.1.3--fb44221536fbccbe'

    output:
    path 'kraken_db_tax'

    script:
    """
    kraken2-build --download-taxonomy --db kraken_db_tax
    """
}

// Download taxon-specific libraries as a starting point for a custom Kraken DB
process KRAKEN_DL_LIB {
    container 'oras://community.wave.seqera.io/library/kraken2:2.1.3--fb44221536fbccbe'

    input:
    val library

    output:
    path "kraken_db_lib_${library}"

    script:
    """
    kraken2-build --download-library $library --db kraken_db_lib_${library}
    """
}

// Add a set of custom genomes in a dir to a custom Kraken DB
process KRAKEN_ADD {
    container 'oras://community.wave.seqera.io/library/kraken2:2.1.3--fb44221536fbccbe'

    input:
    path genome_dir
    path tax_dir
    path lib_dir

    output:
    path "kraken_db_unbuilt"

    script:
    """
    mkdir -p kraken_db_unbuilt/library
    cp -r ${tax_dir}/taxonomy kraken_db_unbuilt
    cp -r ${lib_dir}/library/* kraken_db_unbuilt/library

    for fasta in ${genome_dir}/*; do
        kraken2-build --add-to-library \$fasta --db kraken_db_unbuilt
    done
    """
}

// Build the final custom Kraken DB
process KRAKEN_BUILD {
    container 'oras://community.wave.seqera.io/library/kraken2:2.1.3--fb44221536fbccbe'

    input:
    path kraken_db_unbuilt

    output:
    path 'kraken_db'

    script:
    """
    mkdir kraken_db
    cp -r $kraken_db_unbuilt/* kraken_db/
    kraken2-build --build --db kraken_db
    """
}

// Run Kraken to assign taxonomy to reads
process KRAKEN_RUN {
    container 'oras://community.wave.seqera.io/library/kraken2:2.1.3--fb44221536fbccbe'
    publishDir "$params.outdir/kraken", mode: 'copy', pattern: '*txt'

    input:
    tuple val(sample_id), path(reads)
    path(db)

    output:
    tuple val(sample_id), path('*fastq'), emit: classified_fq
    path '*main.txt', emit: main        // Classifications for each read, e.g. for Pavian
    path '*report.txt', emit: report    // Summary, e.g. for MultiQC

    script:
    """
    kraken2 \
        --threads $task.cpus \
        --db ${db} \
        --classified-out ${sample_id}#.fastq \
        --gzip-compressed \
        --paired \
        --report ${sample_id}_report.txt \
        ${reads[0]} \
        ${reads[1]} \
        > ${sample_id}_main.txt
    """
}

// Assemble reads
process ASSEMBLY {
    container 'oras://community.wave.seqera.io/library/spades:3.15.5--5ae53542733e7564'

    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "${sample_id}.fasta"
    path "${sample_id}_scaffolds.fasta"
    path "${sample_id}_spades.log"
    
    script:
    def memory_gb = MemoryUnit.of("${task.memory}").toUnit('GB')
    """
    spades.py \
        -1 ${reads[0]} \
        -2 ${reads[1]} \
        -o outdir \
        --only-assembler \
        --meta \
        --memory $memory_gb
    
    mv outdir/contigs.fasta ${sample_id}.fasta
    mv outdir/scaffolds.fasta ${sample_id}_scaffolds.fasta
    mv outdir/spades.log ${sample_id}_spades.log
    """
}

process MULTIQC {
    container 'oras://community.wave.seqera.io/library/multiqc:1.22.1--ac0a91c1ae1c160c'
    publishDir "$params.outdir/multiqc", mode: 'copy'

    input:
    path multiqc_input

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
