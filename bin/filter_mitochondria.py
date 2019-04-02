#!/usr/bin/env python3

import sys
import argparse
from collections import namedtuple
from collections import defaultdict
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


def cli(prog, args):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Renames fasta sequences using a sane format.",
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
            "The coverage of the shorter sequence required."
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
    matches = filter_paf(paf, args.coverage)

    to_exclude = {match.query for match in matches if match.is_filtered}

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
    ]

    print("\t".join(columns), file=args.table)
    for match in matches:
        print(
            "\t".join([str(getattr(match, col)) for col in columns]),
            file=args.table
        )
    return


if __name__ == "__main__":
    main()
