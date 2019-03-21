# mitoflow

A short and sweet pipeline for assembling mitochondria and filtering
reads prior to genome assembly. Intended for use with Illumina paired
reads.

Note that in practice, I haven't really found a large difference between
assemblies produced with or without filtering mitochondrial reads.
That is, not with SPAdes in _Parastagonospora nodorum_.
Given the small amount of mitochondrial sequence genuinely present in the
nuclear genome (based on long read assemblies), I play it safe and assemble
with the full read set. I then use these assemblies to fish out the
mitochondrial contigs using [mashmap](https://github.com/marbl/MashMap).


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
