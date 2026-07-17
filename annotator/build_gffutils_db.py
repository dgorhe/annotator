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
import sys

import gffutils

from .index_store import build_annotator_index


def build_db(
        annotation_file, db_file, force=True, disable_infer_genes=False,
        disable_infer_transcripts=False, build_index=False, index_chroms=None,
        gtf_format='gencode', append_chr=False
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
            chroms=index_chroms or [],
            gtf_format=gtf_format,
            append_chr=append_chr,
            force=True,
        )
    return db


def main():
    """Run the gffutils database builder CLI.

    Parses command-line arguments and calls build_db with optional annotator
    index generation.
    """
    # Setup argument parser
    parser = ArgumentParser()

    parser.add_argument(
        "--annotation_file",
        dest="annotation_file",
        help="gff or gtf file",
        required=True
    )
    parser.add_argument(
        "--db_file",
        dest="db_file",
        help="database file output",
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
        help="Limit annotator index to these chromosomes only.",
        nargs='+',
        default=[]
    )
    args = parser.parse_args()

    annotation_file = args.annotation_file
    db_file = args.db_file
    force = args.force
    disable_infer_genes = args.disable_infer_genes
    disable_infer_transcripts = args.disable_infer_transcripts

    build_db(
        annotation_file, db_file, force,
        disable_infer_genes, disable_infer_transcripts,
        build_index=args.build_annotator_index,
        index_chroms=args.limit_chroms_to,
        gtf_format=args.gtf_format,
        append_chr=args.append_chr,
    )


if __name__ == "__main__":
    main()
