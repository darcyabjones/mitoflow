#!/usr/bin/env python3

import sys
import argparse
from itertools import chain
from collections import namedtuple
from collections import defaultdict
from statistics import median
from Bio import SeqIO
from Bio.Alphabet import generic_dna

__version__ = "0.0.1"


PAF = namedtuple(
    "PAF",
    [
        "query",
        "query_len",
        "query_start",
        "query_end",
        "strand",
        "target",
        "target_len",
        "target_start",
        "target_end",
        "nmatches",
        "ali_len",
        "remainder",
    ]
)

PBCov = namedtuple(
    "PBCov",
    ["seqid", "pos", "depth"]
)

Filtered = namedtuple(
    "Filtered",
    [
        "query",
        "target",
        "query_len",
        "target_len",
        "query_start",
        "query_end",
        "target_start",
        "target_end",
        "query_coverage",
        "query_missing",
        "is_filtered",
    ]
)


def parse_paf(handle):
    columns = [
        ("query", str),
        ("query_len", int),
        ("query_start", int),
        ("query_end", int),
        ("strand", str),
        ("target", str),
        ("target_len", int),
        ("target_start", int),
        ("target_end", int),
        ("nmatches", int),
        ("ali_len", int),
    ]

    for line in handle:
        sline = line.strip().split("\t")

        row = {}
        for (col, fn), val in zip(columns, sline):
            row[col] = fn(val)

        row["remainder"] = sline[len(columns):]

        yield PAF(**row)
    return


def parse_pbcov(handle):
    columns = [
        ("seqid", str),
        ("pos", int),
        ("depth", int),
    ]

    for line in handle:
        sline = line.strip().split("\t")
        row = {}
        for (col, fn), val in zip(columns, sline):
            row[col] = fn(val)

        yield PBCov(**row)
    return


def format_missing(missing):
    if len(missing) == 0:
        return ""

    missing = sorted(list(missing))
    blocks = []

    start = missing[0]
    prev = start
    for i in missing[1:]:
        if i != prev + 1:
            if prev == start:
                block = str(start)
            else:
                block = "{}-{}".format(start, prev)
            blocks.append(block)
            start = i

        prev = i

    if prev == start:
        block = str(start)
    else:
        block = "{}-{}".format(start, prev)
    blocks.append(block)
    return ",".join(blocks)


def filter_paf(paf, cov):
    out = []
    candidates = defaultdict(list)

    for aln in paf:
        candidates[aln.query].append(aln)

    for candidate, alns in candidates.items():
        covered_bases = set()
        candidate_len = alns[0].query_len
        for aln in alns:
            covered_bases.update(set(range(aln.query_start, aln.query_end)))

        this_cov = len(covered_bases) / candidate_len
        this_missing = set(range(candidate_len)).difference(covered_bases)
        is_filtered = this_cov >= cov

        for aln in alns:
            filt = Filtered(
                query=candidate,
                target=aln.target,
                query_len=candidate_len,
                target_len=aln.target_len,
                query_start=aln.query_start,
                query_end=aln.query_end,
                target_start=aln.target_start,
                target_end=aln.target_end,
                query_coverage=this_cov,
                query_missing=format_missing(this_missing),
                is_filtered=is_filtered,
            )
            out.append(filt)

    return out


def coverages(pbcovs):
    depth_dict = defaultdict(list)
    for cov in pbcovs:
        depth_dict[cov.seqid].append(cov.depth)
    return depth_dict


def coverage_medians(covs):
    out = dict()
    for seqid, depths in covs.items():
        out[seqid] = median(depths)
    return out


def percentile(vec, p):
    from math import floor
    svec = sorted(vec)
    index = floor(len(svec) * p)
    return svec[index]


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="""
        Filters mitochondrial contigs from a genome assemblyusing contig
        alignment and read coverage information.

        Given an alignment (in paf format, e.g. from minimap2) between some
        contigs and a mitochondrial genome assembly, and a per-base coverage
        file e.g. from `bedtools genomecov -d -ibam my.bam` filters out
        contigs based on alignment coverage of the contigs and the median read
        depth coverage of the aligned reads. E.g. if contigs are almost
        entirely contained within the mitochondrial genome, then that threshold
        would pass. Then if the contig in question has higher median read
        coverage than a specified percentile (ie. it is multicopy) then the
        contig is filtered out.

        Note that the contigs to filter should always be the query in the
        minimap job (and the mitochondiral genome the reference).
        Note also that things like percent identity etc are determined by the
        aligners rather than us doing post-filtering.
        """
    )

    parser.add_argument(
        "-i", "--input",
        type=argparse.FileType('r'),
        required=True,
        help="Input paf file. Use '-' for stdin (default)."
    )

    parser.add_argument(
        "-f", "--fasta",
        type=argparse.FileType('r'),
        required=True,
        help="Input fasta file."
    )

    parser.add_argument(
        "-p", "--per-base-cov",
        type=argparse.FileType('r'),
        required=True,
        help="Per base coverage."
    )

    parser.add_argument(
        "-g", "--genomic",
        type=argparse.FileType('w'),
        required=True,
        help="Output genomic fasta file."
    )

    parser.add_argument(
        "-m", "--mito",
        type=argparse.FileType('w'),
        required=True,
        help="Output mitochondrial fasta file."
    )

    parser.add_argument(
        "-t", "--table",
        type=argparse.FileType('w'),
        default=None,
        help=(
            "Write a tab separated file showing the alignment stats."
        )
    )

    parser.add_argument(
        "-c", "--coverage",
        type=float,
        default=0.95,
        help=(
            "The coverage of the contig sequence required."
        )
    )

    parser.add_argument(
        "-e", "--percentile",
        type=float,
        default=0.99,
        help=(
            "The median percentile of the read coverage required. "
            "Generally, you'll want this value to be quite high "
            "(e.g 0.99 - 0.999)."
        )
    )

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s " + __version__
    )

    parsed = parser.parse_args(args)
    return parsed


def main():
    args = cli(prog=sys.argv[0], args=sys.argv[1:])

    seqs = SeqIO.parse(
        args.fasta,
        format="fasta",
        alphabet=generic_dna
    )

    paf = parse_paf(args.input)
    pbcovs = coverages(parse_pbcov(args.per_base_cov))

    total_pbcov = chain.from_iterable(pbcovs.values())
    percentile_thres = percentile(total_pbcov, args.percentile)

    pbcov_medians = coverage_medians(pbcovs)

    matches = filter_paf(paf, args.coverage)

    to_exclude = {
        match.query
        for match
        in matches
        if match.is_filtered and pbcov_medians[match.query] > percentile_thres
    }

    genomic_seqs = []
    mito_seqs = []
    for seq in seqs:
        if seq.id in to_exclude:
            mito_seqs.append(seq)
        else:
            genomic_seqs.append(seq)

    SeqIO.write(genomic_seqs, args.genomic, format="fasta")
    SeqIO.write(mito_seqs, args.mito, format="fasta")

    columns = [
        "query",
        "target",
        "query_len",
        "target_len",
        "query_start",
        "query_end",
        "target_start",
        "target_end",
        "query_coverage",
        "query_missing",
        "is_filtered",
        "pbcov_median",
        "pbcov_threshold"
    ]

    print("\t".join(columns), file=args.table)
    for match in matches:
        args.table.write(
            "\t".join([str(getattr(match, col)) for col in columns[:-3]]),
        )
        args.table.write(
            "\t{}\t{}\t{}\n".format(
                match.query in to_exclude,
                pbcov_medians[match.query],
                percentile_thres,
            )
        )
    return


if __name__ == "__main__":
    main()
