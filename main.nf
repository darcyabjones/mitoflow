#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # mitoflow

    A short and sweet pipeline for assembling mitochondria and filtering
    reads prior to genome assembly. Intended for use with Illumina paired
    reads.

    ## Usage

    ```bash
    cat <<EOF > reads.tsv
    sample	read1_file	read2_file
    one	pair1_R1.fastq.gz	pair1_R2.fastq.gz
    one	pair2_R1.fastq.gz	pair2_R2.fastq.gz
    two other_R1.fastq.gz	other_R2.fastq.gz
    EOF

    nextflow run -resume darcyabjones/mitoflow \
      --reference mito.fasta \
      --seed "seeds/*.fasta" \
      --asm_table reads.tsv \
      --filter_table reads.tsv \
      --read_length 150 \
      --insert_size 320
    ```

    ## Mandatory Arguments
    
    ```
    param                   | description
    ---------------------------------------------------------------------------
    `--asm_table <tsv>`     | A table mapping fastq read pairs to samples
                            | (See Tables section).

    `--reference <fasta>`   | A fasta formatted reference mitochondrial
                            | assembly from a closely related isolate.
                            | Multiple files can be provided using glob
                            | patterns.

    `--read_length <int>`   | The length of the fastq reads. This can also
                            | be provided in the tsv provided by
                            | `--asm_table` (See Tables section).

    `--insert_length <int>` | The insert size of the fastq pairs. This can
                            | also be provided in the tsv provided by
                            |`--asm_table` (See Tables section).
    ```

    ## Options

    ```
    param                  | default     | description
    ---------------------------------------------------------------------------
    `--filter_table <tsv>` | none        | A table of reads to filter with
                           |             | their corresponding samples
                           |             | (See tables section). If not
                           |             | provided will skip read filtering
                           |             | steps.

    `--seed <fasta>`       | --reference | Fasta formatted sequences to seed
                           |             | the mitochondrial assembly. This
                           |             | could be a mitochondrial gene or
                           |             | assembly from a closely related
                           |             | isolate. Multiple fasta files can
                           |             | be specified using a glob pattern.

    `--min_size <int>`     | 12000       | The minimum size (bp) of the
                           |             | mitochondrial assemblies.

    `--max_size <int>`     | 100000      | The maximum size (bp) of the
                           |             | mitochondrial assemblies

    `--kmer <int>`         | 39          | The K-mer size (bp) to use for the
                           |             | NOVOplasty assembly
    ```

    ## Tables

    Because individual samples are often sequenced in multiple runs,
    and given in multiple fastq pairs, mitoflow takes input to `--asm_table`
    and `--filter_table` as tab separated values (tsv) files. Three columns
    are mandatory for both tables: `sample`, `read1_file`, and `read2_file`.
    The column order is not important, but a header line **must** be included.
    Filenames should be either absolute paths or relative to the executing
    path. The `sample` column is the factor how the fastq pairs will be
    grouped for assemblies and for deciding which assemblies to filter against.
    It will also be the base name of the resultant assemblies.

    `--asm_table` can also use two additional columns `read_length` and
    `insert_size`, which will override the `--read_length` or `--insert_size`
    parameters for this sample. Note that if you don't provide `--read_length`
    or `--insert_size` _all_ rows must have corresponding values in these
    columns. `--filter_table` can also use the additional column `merged_file`
    which can be a single "stitched" fastq file to be filtered separately from
    the pairs.

    The `--asm_table` and `--filter_table` can be the same file. Additional
    columns will be ignored. However, be aware that NOVOplasty can't use
    merged reads, so don't use your pre-filtered reads for the `--asm_table`.
    See the `examples` folder in the github repo for examples.

    ## Outputs

    * `assemblies/*_mitochondrial.fasta`:
        The assembled mitochondria per sample. 

    * `assemblies/*_log.txt`:
        Logs from NOVOPlasty for assemblies.

    * `alignments/*.{delta,mcoords,...}`:
        MUMmer files from alignment between assemblies and `reference`.

    * `filtered_reads/<sample>`:
        Filtered fastq pairs named same as input.

    * `filtered_reads/*_mitochondrial.fastq.gz`:
        Filtered fastq pairs aligning to Mitochondria.

    * `filtered_reads/*_scafstats.txt`:
        BBSplit statistics containing number of reads aligned to different
        scaffolds.

    * `filtered_reads/*_refstats.txt`:
        BBSplit statistics containing number of reads aligned to different
        references (either `--reference` or assembly for this sample).
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
params.read_length = false
params.insert_size = false


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

    // Take the read pairs from the filter tables, split 
    // and remove any rows where all relevant filenames aren't
    // given.
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
    log.error "You've reached a point in the code that shouldn't be reached."
    log.error "please raise a bug report."
    exit 1
}


if ( !params.read_length ) {
    log.warn "You didn't provide a `read_length` parameter." 
    log.warn "Make sure _all_ rows in `asm_table` have the `read_length` column set."
    log.warn "Otherwise will raise a runtime error at the assembly step."
}


if ( !params.insert_size ) {
    log.warn "You didn't provide a `insert_size` parameter." 
    log.warn "Make sure _all_ rows in `asm_table` have the `insert_size` column set."
    log.warn "Otherwise will raise a runtime error at the assembly step."
}


// END OF PARAMETER VALIDATION


process assembleMito {
    label "novoplasty"
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

    if ( read_length.any { it != null } ) {
        read_length = (read_length - null)[0]
    } else if ( params.read_length ) {
        read_length = params.read_length
    } else {
        log.error "ERROR processing assembly for sample: ${name}."
        log.error "The `read_length` parameter is not set or provided in the `asm_table`."
        log.error "Please provide one of those and re-run :)."
        exit 1
    }

    if ( insert_size.any { it != null } ) {
        insert_size = (insert_size - null)[0]
    } else if ( params.insert_size ) {
        insert_size = params.insert_size
    } else {
        log.error "ERROR processing assembly for sample: ${name}."
        log.error "The `insert_size` parameter is not set or provided in the `asm_table`."
        log.error "Please provide one of those and re-run :)."
        exit 1
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

    // Split the read files to have value names (extensions) and simplenames (noextensions).
    // We do this so that we can have the outgoing reads named the same as the ingoing reads.
    // We also join the mitochondrial assembly for this isolate into the channel.
    joined4ReadFiltering = pairTable
        .map { n, b, f, r -> [n, b, f.name, r.name, f.simpleName, r.simpleName, f, r] }
        .join( mitoAssemblies4ReadFiltering.map { n, a, s -> [n, a] }, by: 0)

    process readFiltering {
        label "bbmap"
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
        # This is just so that the input reads have the same extension 
        # as the original file.
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


    // As with `joined4ReadFiltering`, split out the read filenames and merge the assembly.
    joined4MergedReadFiltering = mergedTable
        .map { n, m -> [n, m.name, m.simpleName, m] }
        .join( mitoAssemblies4MergedReadFiltering.map { n, a, s -> [n, a] }, by: 0)

    process mergedReadFiltering {
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
        # This is just so that the input reads have the same extension 
        # as the original file.
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
