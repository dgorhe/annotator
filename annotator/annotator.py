#!/usr/bin/env python
"""Command-line entry point for annotating BED files.

Parses CLI arguments and delegates to annotate_bed.annotate_bed, using
built-in or user-supplied transcript and gene annotation priority lists.
"""

from argparse import ArgumentParser
import sys

from . import annotate_bed


TRANSCRIPT_PRIORITY = [
    ['miRNA', 'exon'],
    ['protein_coding', 'CDS'],
    ['protein_coding', 'start_codon'],
    ['protein_coding', 'stop_codon'],
    ['protein_coding', '3utr'],
    ['protein_coding', 'three_prime_utr'],
    ['protein_coding', '5utr'],
    ['protein_coding', 'five_prime_utr'],
    ['protein_coding', 'unclassified_utr'],
    ['protein_coding', 'proxintron500'],
    ['protein_coding', 'distintron500'],
    ['non_coding', 'exon'],
    ['protein_coding', 'Selenocysteine'],
    ['non_coding', 'proxintron500'],
    ['non_coding', 'distintron500'],
    ['non_coding', 'transcript'],
    ['non_coding', 'gene'],
    ['non_coding', 'Selenocysteine'],
]

GENE_PRIORITY = [
    ['miRNA', 'exon'],
    ['protein_coding', 'CDS'],
    ['protein_coding', 'start_codon'],
    ['protein_coding', 'stop_codon'],
    ['protein_coding', '3utr'],
    ['protein_coding', 'three_prime_utr'],
    ['protein_coding', '5utr'],
    ['protein_coding', 'five_prime_utr'],
    ['protein_coding', 'unclassified_utr'],
    ['non_coding', 'exon'],
    ['protein_coding', 'proxintron500'],
    ['non_coding', 'proxintron500'],
    ['protein_coding', 'distintron500'],
    ['non_coding', 'distintron500'],
    ['protein_coding', 'Selenocysteine'],
    ['non_coding', 'transcript'],
    ['non_coding', 'gene'],
    ['non_coding', 'Selenocysteine'],
]


def parse_annotation_priority(priority_file, delim=','):
    """Load an annotation priority list from a delimited text file.

    Each non-empty line must contain at least two fields: a transcript or
    gene type and a feature type, in that order.

    Args:
        priority_file (str): Path to the priority file.
        delim (str): Field delimiter. Defaults to comma.

    Returns:
        list: List of [type, featuretype] pairs, one per line.
    """
    priority = []
    with open(priority_file, 'r') as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split(delim)
            if len(fields) < 2:
                raise ValueError(
                    "Priority file {} requires at least two fields per line: {!r}".format(
                        priority_file, line.rstrip()
                    )
                )
            priority.append([fields[0], fields[1]])
    return priority


def main():
    """Run the annotator CLI.

    Parses command-line arguments, resolves default or file-based priority
    lists, and invokes annotate_bed for each input BED file.
    """
    parser = ArgumentParser()

    parser.add_argument(
        "--output",
        dest="output",
        help="output file(s)",
        required=False,
        nargs='+'
    )
    parser.add_argument(
        "--input",
        dest="input",
        help="input bed6 file(s)",
        required=True,
        nargs='+'
    )
    parser.add_argument(
        "--gtfdb",
        dest="gtfdb",
        help="gtf database file create from gffutils (Currently only "
             "GTF files are allowed)",
        required=True
    )
    parser.add_argument(
        "--transcript-priority-file",
        dest="transcript_priority_file",
        required=False,
        help='For each transcript (grouped by genes), return the transcript'
             ' using this priority',
        default=None
    )
    parser.add_argument(
        "--gene-priority-file",
        dest="gene_priority_file",
        required=False,
        help='For each gene, return the gene region using this priority',
        default=None
    )
    parser.add_argument(
        "--unstranded",
        dest="unstranded",
        required=False,
        action='store_true',
        help='ALLOW unstranded - will still search for correctly'
             ' stranded features first (searches positive first if a BED3)',
        default=False
    )
    parser.add_argument(
        "--append-chr",
        dest="append_chr",
        required=False,
        action='store_true',
        help='If your bedfile is formatted as something like [chrI] but the '
             'gtf database expects [I], we can retroactively append [chr] '
             'for you.',
        default=False
    )
    parser.add_argument(
        "--format",
        dest="gtf_format",
        required=False,
        help='sets default GTF nomenclature to either [gencode] or ensembl/ncbi',
        default='gencode'
    )
    parser.add_argument(
        "--limit-chroms-to",
        dest="limit_chroms_to",
        required=False,
        help="limit to annotating to these chromosomes "
             "only, saves time/memory. When omitted, chromosomes are "
             "inferred from the input BED unless --all-chroms is set.",
        nargs='+',
        default=[]
    )
    parser.add_argument(
        "--all-chroms",
        dest="auto_chroms",
        action="store_false",
        help="Load annotation for every chromosome in the database/index "
             "instead of inferring contigs from the input BED.",
    )
    parser.add_argument(
        "--no-triage-buckets",
        dest="triage_buckets",
        action="store_false",
        help="When using .annotator.idx, load all buckets for selected "
             "chromosomes instead of only buckets overlapping the BED.",
    )
    parser.add_argument(
        "--lazy-buckets",
        dest="lazy_buckets",
        action="store_true",
        help="Load index buckets on demand with an LRU cache (single-core "
             "only). Useful for sparse BEDs when peak memory matters.",
    )
    parser.add_argument(
        "-c", "--cores",
        dest="cores",
        type=int,
        default=1,
        help="number of CPU cores to use for annotation (default: 1)",
    )
    parser.add_argument(
        "--no-index",
        dest="use_index",
        action='store_false',
        help="Do not load a pre-built .annotator.idx companion file even if "
             "one exists next to the gtfdb.",
    )
    args = parser.parse_args()

    if args.cores < 1:
        parser.error("--cores must be a positive integer")

    input_bed_files = args.input
    output_annotated_files = args.output
    if output_annotated_files is None:
        output_annotated_files = [
            bed_file + '.annotated' for bed_file in input_bed_files
        ]
    elif len(output_annotated_files) != len(input_bed_files):
        parser.error("number of --output files must match number of --input files")

    gtfdb_file = args.gtfdb
    chroms = args.limit_chroms_to
    stranded = not args.unstranded
    gtf_format = args.gtf_format
    append_chr = args.append_chr
    cores = args.cores
    fuzzy = 0

    if args.transcript_priority_file is not None:
        transcript_priority = parse_annotation_priority(
            args.transcript_priority_file)
    else:
        transcript_priority = [list(item) for item in TRANSCRIPT_PRIORITY]

    if args.gene_priority_file is not None:
        gene_priority = parse_annotation_priority(args.gene_priority_file)
    else:
        gene_priority = [list(item) for item in GENE_PRIORITY]

    annotate_bed.annotate_bed(
        gtfdb_file, input_bed_files, output_annotated_files, stranded, chroms,
        transcript_priority, gene_priority, gtf_format, append_chr,
        fuzzy, cores, use_index=args.use_index,
        auto_chroms=args.auto_chroms,
        triage_buckets=args.triage_buckets,
        lazy_buckets=args.lazy_buckets,
    )


if __name__ == "__main__":
    main()
