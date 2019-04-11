#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

def helpMessage() {
    log.info"""
    # mitoflow - postprocessing

    Filter out mitochondrial contigs from an assembly.

    ## Usage

    ```bash
    nextflow run -resume ./post.nf "*{_genomic.fasta,_mitochondrial.fasta,.bam}"
    ```

    Note that because of the way the glob works, you need something unique
    in the bracket expansion for the pattern to match properly.
    E.G. Say you had the files `genome.fasta` and `genome_mitochondrial.fasta`
    the pattern, `*{.fasta,_mitochondrial.fasta,.bam} doesn't work because
    `genome_mitochondrial.fasta` will match twice.
    It's annoying, but you might just have to rename your files. Sorry!

    Requirements:
    minimap2
    bedtools
    Python3 with Biopython installed
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.genomes = false
params.coverage = 0.95
params.percentile = 0.992


if ( params.genomes ) {
    genomes = Channel.fromFilePairs(
        params.genomes,
        checkIfExists: true,
        size: 3,
        type: "file",
        flat: true
    )
} else {
    log.info "Hey I need some assemblies and mitochondria to filter out."
    exit 1
}


genomes.into {
    genomes4Align;
    genomes4Coverage;
}


// Identity filtering is just done as part of the minimap alignment, rather
// than filtering out the results.
// I'll need to add some options here to customise, but the preset value works
// pretty well.
process align {
    label "minimap2"
    tag { name }

    input:
    set val(name), file(asm), file(mito_asm), file(bam) from genomes4Align

    output:
    set val(name), file("matches.paf"), file(asm) into alignedGenomes

    """
    minimap2 \
      -x asm20 \
      -o matches.paf \
      ${mito_asm} \
      ${asm}
    """
}


process coverage {
    label "bedtools"

    input:
    set val(name), file(asm), file(mito_asm), file(bam) from genomes4Coverage

    output:
    set val(name), file("per_base_cov.tsv") into pbCoverages

    """
    bedtools genomecov -ibam ${bam} -d > per_base_cov.tsv
    """
}


process filter {

    label "python3"
    tag { name }
    publishDir "${params.outdir}/filtered_assemblies"

    input:
    set val(name), file("matches.paf"), file("genome.fasta"),
        file("per_base_cov.tsv") from alignedGenomes
            .join(pbCoverages, by: 0)

    output:
    set val(name),
        file("${name}_matches.tsv"),
        file("${name}_genomic.fasta"),
        file("${name}_mitochondrial.fasta") into filteredGenomes

    """
    filter_mitochondria.py \
      -i matches.paf \
      -f genome.fasta \
      -p per_base_cov.tsv \
      -g ${name}_genomic.fasta \
      -m ${name}_mitochondrial.fasta \
      -t ${name}_matches.tsv \
      --coverage ${params.coverage} \
      --percentile ${params.percentile}
    """
}
