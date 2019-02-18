#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    =================================
    qcflow
    =================================
    Usage:
    abaaab
    Mandatory Arguments:
      --fastq              description
      --seed
      --references         description
    Options:
    Outputs:
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

// Default parameter values

params.asm_table = false
params.filter_table = false
params.seed = false
params.reference = false
params.min_size = 12000
params.max_size = 100000
params.kmer = 39
params.read_length = 125
params.insert_size = 350


// Find common prefix between 2 strings.
// Used here to emulate basename given from the Channel.fromPairs glob
static String lcp(String r1, String r2){
    def l = [r1.toList(), r2.toList()]
        .transpose()
        .takeWhile {it[0] == it[1]}
        .collect {it[0]}
        .join("")
    return l.replaceAll(/[-\._]$/, "")
}


if ( params.asm_table ) {
    asmTable = Channel.fromPath(params.asm_table)
        .splitCsv(by: 1, sep: '\t', header: true)
        .map { [it.sample, file(it.read1_file), file(it.read2_file), it.read_length, it.insert_size] }
} else {
    log.info "Hey I need some reads to assemble into a mitochondrial genome."
    exit 1
}

if ( params.filter_table ) {
    pairTable = Channel.fromPath(params.filter_table)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.read1_file != null && it.read2_file != null }
        .map { [it.sample, file(it.read1_file), file(it.read2_file)] }
        .map { s, f, r -> [s, lcp(f.name, r.name), f, r] }

    mergedTable = Channel.fromPath(params.filter_table)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.merged_file != null }
        .map { [it.sample, file(it.merged_file)] }
}

if ( params.reference ) {
    reference = Channel.fromPath(params.reference, checkIfExists: true, type: "file")
        .collectFile(name: 'reference.fasta', newLine: true, sort: "deep")
        .first()
} else {
    log.info "Hey I need a mitochondrial reference."
    exit 1
}


if ( params.seed ) {
    seed = Channel.fromPath(params.seed, checkIfExists: true, type: "file")
        .collectFile(name: 'seed.fasta', newLine: true, sort: "deep")
        .first()
} else if ( params.reference ) {
    seed = reference
} else {
    log.info "Hey I need some fasta sequences to seed the mitochondrial assembly from."
    exit 1
}

// END OF PARAMETER VALIDATION


process assembleMito {
    label "posix"
    label "novoplasty"
    label "perl"
    label "medium_task"

    publishDir "${params.outdir}/assemblies"

    tag { name }

    input:
    file "ref.fasta" from reference
    file "seed.fasta" from seed
    set val(name), file("*R1.fastq.gz"), file("*R2.fastq.gz"),
        val(read_length), val(insert_size) from asmTable
            .groupTuple(by: 0)

    output:
    set val(name), file("${name}_mitochondrial.fasta"),
        file("${name}_log.txt") into mitoAssemblies

    script:

    if ( read_length.every { it == null } ) {
        read_length = params.read_length
    } else {
        read_length = (read_length - null)[0]
    }

    if ( insert_size.every { it == null } ) {
        insert_size = params.insert_size
    } else {
        insert_size = (insert_size - null)[0]
    }

    """
    cat *R1.fastq.gz > forward.fastq.gz
    cat *R2.fastq.gz > reverse.fastq.gz

    cat << EOF > config.txt
Project:
-----------------------
Project name          = ${name}
Type                  = mito
Genome Range          = ${params.min_size}-${params.max_size}
K-mer                 = ${params.kmer}
Extended log          = 0
Save assembled reads  = no
Seed Input            = seed.fasta
Reference sequence    = ref.fasta
Variance detection    = no

Dataset 1:
-----------------------
Read Length           = ${read_length}
Insert size           = ${insert_size}
Platform              = illumina
Single/Paired         = PE
Forward reads         = forward.fastq.gz
Reverse reads         = reverse.fastq.gz

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
Use Quality Scores    = no
EOF

    NOVOPlasty.pl -c config.txt

    # Rename for better sorting and consistency.
    if [[ -f  Circularized_assembly_1_${name}.fasta ]]; then
        mv Circularized_assembly_1_${name}.fasta ${name}_mitochondrial.fasta
    elif [[ -f  Uncircularized_assemblies_1_${name}.fasta ]]; then
        mv Uncircularized_assemblies_1_${name}.fasta ${name}_mitochondrial.fasta
    fi

    mv log_${name}.txt ${name}_log.txt

    rm forward.fastq.gz
    rm reverse.fastq.gz
    """
}

mitoAssemblies.into {
    mitoAssemblies4ReadFiltering;
    mitoAssemblies4MergedReadFiltering;
    mitoAssemblies4Alignment;
}


process alignGenomes {
    label "mummer"
    label "medium_task"

    publishDir "${params.outdir}/alignments"

    tag { name }

    input:
    file "ref.fasta" from reference
    set val(name), file(asm), file(stats) from mitoAssemblies4Alignment

    output:
    file "${name}.delta" into deltas
    file "${name}.snps" into snps
    file "${name}.1coords" into onecoords
    file "${name}.1delta" into onedeltas
    file "${name}.mcoords" into mcoords
    file "${name}.mdelta" into mdeltas
    file "${name}.qdiff" into qdiff
    file "${name}.rdiff" into rdiff
    file "${name}.report" into dnadiffReports
    file "${name}.unqry" optional true into unqrys

    """
    nucmer --prefix="${name}" "ref.fasta" "${asm}"
    dnadiff --prefix="${name}" --delta "${name}.delta"
    """
}


if ( params.filter_table ) {
    joined4ReadFiltering = pairTable
        .map { n, b, f, r -> [n, b, f.name, r.name, f.simpleName, r.simpleName, f, r] }
        .join( mitoAssemblies4ReadFiltering.map { n, a, s -> [n, a] }, by: 0)

    process readFiltering {
        label "bbmap"
        label "java"
        label "biggish_task"

        publishDir "${params.outdir}/filtered_reads"

        tag { name }

        input:
        file "ref.fasta" from reference
        set val(name), val(base_name), val(read1_name), val(read2_name),
            val(read1_simplename), val(read2_simplename),
            file("fwd"), file("rev"), file("asm.fasta") from joined4ReadFiltering

        output:
        set val(name), file(read1_name), file(read2_name) into filteredReads
        set val(name), file("${read1_simplename}_mitochondrial.fastq.gz"),
            file("${read2_simplename}_mitochondrial.fastq.gz") into mitochondrialReads
        set val(name), file("${base_name}_scafstats.txt"),
            file("${base_name}_scafstats.txt") into filteredStats

        """
        FWD="in_${read1_name}"
        ln -sf fwd "\${FWD}"
        REV="in_${read2_name}"
        ln -sf rev "\${REV}"

        bbsplit.sh \
            -Xmx${task.memory.toGiga()}g \
            t=${task.cpus} \
            ref_ref="ref.fasta" \
            ref_${name}="asm.fasta" \
            in1="\${FWD}" \
            in2="\${REV}" \
            outu1="${read1_name}" \
            outu2="${read2_name}" \
            outm1="${read1_simplename}_mitochondrial.fastq.gz" \
            outm2="${read2_simplename}_mitochondrial.fastq.gz" \
            scafstats="${base_name}_scafstats.txt" \
            refstats="${base_name}_refstats.txt"

        # Remove the index because not needed
        rm -rf -- ref
        """
    }


    joined4MergedReadFiltering = mergedTable
        .map { n, m -> [n, m.name, m.simpleName, m] }
        .join( mitoAssemblies4MergedReadFiltering.map { n, a, s -> [n, a] }, by: 0)

    process mergedReadFiltering {
        label "bbmap"
        label "java"
        label "biggish_task"

        publishDir "${params.outdir}/filtered_reads"

        tag { name }

        input:
        file "ref.fasta" from reference
        set val(name), val(reads_name), val(reads_simplename),
            file("reads"), file("asm.fasta") from joined4MergedReadFiltering

        output:
        set val(name), file(reads_name) into filteredMergedReads
        set val(name), file("${reads_simplename}_mitochondrial.fastq.gz") into mitochondrialMergedReads
        set val(name), file("${reads_simplename}_scafstats.txt"),
            file("${reads_simplename}_scafstats.txt") into filteredMergedStats

        """
        READS="in_${reads_name}"
        ln -sf reads "\${READS}"

        bbsplit.sh \
            -Xmx${task.memory.toGiga()}g \
            t=${task.cpus} \
            ref_ref="ref.fasta" \
            ref_${name}="asm.fasta" \
            in="\${READS}" \
            outu="${reads_name}" \
            outm="${reads_simplename}_mitochondrial.fastq.gz" \
            scafstats="${reads_simplename}_scafstats.txt" \
            refstats="${reads_simplename}_refstats.txt"

        # Remove the index because not needed
        rm -rf -- ref
        """
    }
}
