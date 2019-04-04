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
    nextflow run -resume ./post.nf "*{_genomic,_mitochondrial}.fasta"
    ```

    Note that because of the way the glob works, you need something unique
    in the bracket expansion for the pattern to match properly.
    E.G. Say you had the files `genome.fasta` and `genome_mitochondrial.fasta`
    the pattern, `*{,_mitochondrial}.fasta` doesn't work because
    `genome_mitochondrial.fasta` will match twice.
    It's annoying, but you might just have to rename your files. Sorry!

    Requirements:
    minimap2
    Python3 with Biopython installed
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

params.reference = false
params.asm_table = false

if ( params.asm_table ) {
    asmTable = Channel.fromPath(params.asm_table)
        .splitCsv(by: 1, sep: '\t', header: true)
        .map { [it.assembly, file(it.read1_file), file(it.read2_file)] }
} else {
    log.info "Hey I need some reads to assemble into a mitochondrial genome."
    exit 1
}


if ( params.reference ) {
    reference = Channel.fromPath(params.reference, checkIfExists: true, type: "file")
        .collectFile(name: 'reference.fasta', newLine: true, sort: "deep")
        .first()
} else {
    log.info "Hey I need a mitochondrial reference."
    exit 1
}


asmTable.into {
    asmTable4Recycle;
    asmTable4ReadAlignment;
}

process recycle {
    label "hced"
    tag { assembly.simpleName }
    publishDir "${params.outdir}/recycled"

    input:
    file assembly from asmTable4Recycle
        .map { a, r1, r2 -> a }
        .unique()
    file reference from reference

    output:
    set val(assembly.name), file("${assembly.baseName}_recycled.fasta") into recycledFasta

    """
    # The reference (one that is not cycled) should be second sequence in file.
    cat ${assembly} > combined.fasta
    # In case no trailing newline.
    echo "" >> combined.fasta
    cat ${reference} >> combined.fasta

    hCED \
      --input-file combined.fasta \
      --output-file recycled.fasta

    # Takes just the first seq from multifasta
    awk 'BEGIN { NSEQS=0 } /^>/ {NSEQS++} NSEQS < 2 {print $0} ' \
      recycled.fasta > ${assembly.baseName}_recycled.fasta
    """
}

recycledFasta.into {
    recycledFasta4MSA;
    recycledFasta4IndexAssembly;
}

process msa {
    label "mafft"
    publishDir "${params.outdir}/msa"

    input:
    file "*.fasta" from recycledFasta4MSA
        .map { n, f -> f}
        .collect()
    file "reference.fa" from reference

    output:
    file "aligned.fasta" into alignedAssemblies

    """
    cat reference.fa > combined.fasta
    cat *.fasta >> combined.fasta

    mafft combined.fasta > aligned.fasta
    """
}


process indexAssembly {
    label "bwa"
    tag { asm_name }

    input:
    set val(asm_name), file(asm) from recycledFasta4IndexAssembly

    output:
    set val(asm_name), file("index") into indexedAssembly

    """
    mkdir index
    cp ${asm} index/ref.fa
    samtools faidx index/ref.fa
    bwa index -p index/ref ${asm}
    """
}


joined4ReadAlignment = indexedAssembly
    .combine(asmTable4ReadAlignment.map { n, r1, r2 -> [n.name, r1, r2]}, by: 0)

process readAlignment {
    label "bwa"
    tag { asm_name }

    input:
    set val(asm_name), file("index"), file(read1_file),
        file(read2_file) from joined4ReadAlignment

    output:
    set val(asm_name), file("aligned.sam") into alignedReads

    """
    bwa mem index/ref ${read1_file} ${read2_file} > aligned.sam
    """
}


process combineAlignments {
    label "samtools"
    tag { asm_name }

    publishDir "read_alignment"

    input:
    set val(asm_name), file("*.sam") from alignedReads
        .groupTuple(by: 0)

    output:
    set val(asm_name), file("${asm_name}.bam"),
        file("${asm_name}.bam.bai") into combinedAlignedReads
    set val(asm_name), file("${asm_name}.idxstats"),
        file("${asm_name}.flagstats"), file("${asm_name}.stats") into alignedStats

    """
    cat *.sam > combined.sam

    samtools sort -O bam -o "${asm_name}.bam" combined.sam
    samtools index "${asm_name}.bam"

    samtools idxstats "${asm_name}.bam" > "${asm_name}.idxstats"
    samtools flagstat "${asm_name}.bam" > "${asm_name}.flagstat"
    samtools stats "${asm_name}.bam" > "${asm_name}.stats"

    rm combined.sam
    """
}
