#!/usr/bin/env python
"""Build a gffutils FeatureDB from GTF/GFF annotation files.

Creates SQLite annotation databases for use with the annotator package and
optionally builds a companion .annotator.idx index for faster BED lookup.
"""

# these two are really a minimum

# from __future__ import print_function
# from __future__ import division

# uncomment from this compatibility import list, as py3/py2 support progresses

# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
# from future.utils import iteritems

from argparse import ArgumentParser

import gffutils

from .index_store import build_annotator_index, index_path_for_db


def build_db(
        annotation_file, db_file, force=True, disable_infer_genes=False,
        disable_infer_transcripts=False, build_index=False, index_chroms=None,
        gtf_format='gencode', append_chr=False, index_file=None,
):
    """Create a gffutils database from a GTF or GFF annotation file.

    Args:
        annotation_file (str): Path to the input GTF or GFF file.
        db_file (str): Output path for the gffutils SQLite database.
        force (bool): If True, overwrite an existing database at db_file.
        disable_infer_genes (bool): If True, do not infer gene features
            missing from the annotation.
        disable_infer_transcripts (bool): If True, do not infer transcript
            features missing from the annotation.
        build_index (bool): If True, also build a .annotator.idx companion
            index after the database is created.
        index_chroms (list, optional): Chromosomes to include in the index.
            Empty or None includes all chromosomes.
        gtf_format (str): GTF nomenclature used when building the index.
        append_chr (bool): If True, prepend 'chr' to contig names when
            building the index.
        index_file (str, optional): Explicit path for the companion index.
            Defaults to ``index_path_for_db(db_file)``.

    Returns:
        gffutils.FeatureDB: The created (or overwritten) database handle.
    """
    db = gffutils.create_db(
        annotation_file, dbfn=db_file, force=force,
        keep_order=True, merge_strategy='merge', sort_attribute_values=True,
        disable_infer_genes=disable_infer_genes,
        disable_infer_transcripts=disable_infer_transcripts
    )
    if build_index:
        build_annotator_index(
            db_file,
            index_file=index_file,
            chroms=index_chroms or [],
            gtf_format=gtf_format,
            append_chr=append_chr,
            force=True,
        )
    return db


def build_index_from_db(
        db_file, index_file=None, index_chroms=None, gtf_format='gencode',
        append_chr=False, force=True):
    """Build a companion ``.annotator.idx`` from an existing gffutils database.

    Args:
        db_file (str): Path to an existing gffutils SQLite database.
        index_file (str, optional): Output index path. Defaults to the
            sidecar path next to ``db_file``.
        index_chroms (list, optional): Chromosomes to include. Empty or None
            indexes all chromosomes (full annotate-only index). A single
            chromosome produces a chrom-sharded slim index.
        gtf_format (str): GTF nomenclature for attribute keys.
        append_chr (bool): If True, prepend ``chr`` when hashing contigs.
        force (bool): If True, overwrite an existing index file.

    Returns:
        str: Path to the written index file.
    """
    return build_annotator_index(
        db_file,
        index_file=index_file or index_path_for_db(db_file),
        chroms=index_chroms or [],
        gtf_format=gtf_format,
        append_chr=append_chr,
        force=force,
    )


def main():
    """Run the gffutils database builder CLI.

    Parses command-line arguments and calls build_db with optional annotator
    index generation, or builds an index-only artifact from an existing DB.
    """
    # Setup argument parser
    parser = ArgumentParser()

    parser.add_argument(
        "--annotation_file",
        dest="annotation_file",
        help="gff or gtf file (required unless --index-only)",
        required=False
    )
    parser.add_argument(
        "--db_file",
        dest="db_file",
        help="database file output (or existing DB when --index-only)",
        required=True
    )
    parser.add_argument(
        "--force",
        action='store_true',
        help="if true, override any existing database file if it exists "
             "in the same path as specified output.",
        required=False,
        default=False
    )
    parser.add_argument(
        "--disable_infer_genes",
        action='store_true',
        help="If True, gffutils won't infer genes (default: False).",
        required=False,
        default=False
    )
    parser.add_argument(
        "--disable_infer_transcripts",
        action='store_true',
        help="If True, gffutils won't infer transcripts (default: False).",
        required=False,
        default=False
    )
    parser.add_argument(
        "--build_annotator_index",
        action='store_true',
        help="After building the gffutils DB, also build a .annotator.idx "
             "companion index for fast BED annotation startup.",
        required=False,
        default=False
    )
    parser.add_argument(
        "--index-only",
        dest="index_only",
        action='store_true',
        help="Build or rebuild only the .annotator.idx from an existing "
             "--db_file (does not create a gffutils DB). Use with "
             "--limit-chroms-to for chrom-sharded slim indexes.",
        required=False,
        default=False
    )
    parser.add_argument(
        "--index-file",
        dest="index_file",
        required=False,
        default=None,
        help="Explicit output path for the companion index (default: "
             "<db_file>.annotator.idx). Useful for chrom-sharded artifacts.",
    )
    parser.add_argument(
        "--format",
        dest="gtf_format",
        required=False,
        help="GTF nomenclature used when building the annotator index.",
        default='gencode'
    )
    parser.add_argument(
        "--append-chr",
        dest="append_chr",
        action='store_true',
        help="Append 'chr' to contig names when building the annotator index.",
        required=False,
        default=False
    )
    parser.add_argument(
        "--limit-chroms-to",
        dest="limit_chroms_to",
        required=False,
        help="Limit annotator index to these chromosomes only (slim / "
             "chrom-sharded annotate-only index).",
        nargs='+',
        default=[]
    )
    args = parser.parse_args()

    if args.index_only:
        if not args.db_file:
            parser.error("--db_file is required with --index-only")
        build_index_from_db(
            args.db_file,
            index_file=args.index_file,
            index_chroms=args.limit_chroms_to,
            gtf_format=args.gtf_format,
            append_chr=args.append_chr,
            force=True,
        )
        return

    if not args.annotation_file:
        parser.error("--annotation_file is required unless --index-only is set")

    build_db(
        args.annotation_file, args.db_file, args.force,
        args.disable_infer_genes, args.disable_infer_transcripts,
        build_index=args.build_annotator_index,
        index_chroms=args.limit_chroms_to,
        gtf_format=args.gtf_format,
        append_chr=args.append_chr,
        index_file=args.index_file,
    )


if __name__ == "__main__":
    main()
