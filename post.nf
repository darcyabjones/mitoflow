#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

def helpMessage() {
    log.info"""
    # mitoflow - postprocessing

    Filter out mitochondrial contigs from an assembly.
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.genomes = false
params.coverage = 0.95


if ( params.genomes ) {
    genomes = Channel.fromFilePairs(
        params.genomes,
        checkIfExists: true,
        size: 2,
        type: "file",
        flat: true
    )
} else {
    log.info "Hey I need some assemblies and mitochondria to filter out."
    exit 1
}


process align {
    label "minimap2"
    tag { name }

    input:
    set val(name), file(asm), file(mito_asm) from genomes

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


process filter {

    label "python3"
    tag { name }
    publishDir "${params.outdir}/filtered_assemblies"

    input:
    set val(name), file("matches.paf"), file(asm) from alignedGenomes

    output:
    set val(name),
        file("${name}_matches.tsv"),
        file("${name}_genomic.fasta"),
        file("${name}_mitochondrial.fasta") into filteredGenomes

    """
    filter_mitochondria.py \
      -i matches.paf \
      -f ${asm} \
      -g ${asm.baseName}_genomic.fasta \
      -m ${asm.baseName}_mitochondrial.fasta \
      -t ${asm.baseName}_matches.tsv \
      --coverage ${params.coverage}
    """
}
