#!/usr/bin/env python
"""Convert tabular files labeled with gene names to gene IDs.

Uses a gffutils database and/or a custom id-to-name mapping file to look up
gene IDs and append them to an input table.
"""

# from __future__ import print_function
# from __future__ import division
#
# # uncomment from this compatibility import list, as py3/py2 support progresses
#
# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# # from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
# from future.utils import iteritems

from argparse import ArgumentParser
import sys

import gffutils
from collections import defaultdict


def _require_pandas():
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "pandas is required for gene_name2id. "
            "Install with: pip install 'annotator[utils]'"
        ) from exc
    return pd


def get_id_from_namerow(row, name2id, name_col):
    """Look up gene ID(s) for a single row's gene name.

    Args:
        row (pandas.Series): Row containing a gene name column.
        name2id (dict): Mapping from gene name to list of gene IDs.
        name_col (str): Column name holding the gene name.

    Returns:
        str: Pipe-delimited gene ID(s), or an empty string if the name is missing.
    """
    try:
        return '|'.join(name2id[str(row[name_col])])
    except KeyError:
        print('Cant find {} in the database!'.format(str(row[name_col])))
        return ''


def add_id_to_dataframe(df, name2id_dict, name_col):
    """Append a ``gene id`` column to a dataframe.

    Args:
        df (pandas.DataFrame): Input table with gene names.
        name2id_dict (dict): Mapping from gene name to list of gene IDs.
        name_col (str): Column name holding the gene name.

    Returns:
        pandas.DataFrame: Copy of ``df`` with a new ``gene id`` column.
    """
    assert 'gene id' not in df.columns
    df['gene id'] = df.apply(get_id_from_namerow, args=[
                             name2id_dict, name_col], axis=1)
    return df


def build_name_to_id_dict(gffdb_file, custom_file):
    """Build a gene-name-to-ID mapping from a gffutils DB and/or custom file.

    Names with more than one associated ID are stored as a list; callers join
    multiple IDs with ``|``.

    Args:
        gffdb_file (str or None): Path to a gffutils database. The database
            must contain ``gene`` features with ``gene_name`` and ``gene_id``
            attributes.
        custom_file (str or None): Tab-delimited file of ``gene_id<TAB>gene_name``
            rows for genes not present in the gffutils database.

    Returns:
        dict: Mapping ``{gene_name: [gene_id, ...]}``.
    """
    # delim = '|'
    name2id = defaultdict(list)
    if gffdb_file is not None:
        try:
            db = gffutils.FeatureDB(gffdb_file)
            for transcript in db.features_of_type('gene'):
                for gene_name in transcript.attributes['gene_name']:
                    for gene_id in transcript.attributes['gene_id']:
                        name2id[gene_name].append(gene_id)

        except Exception as e:
            print(e)
            print("Could not use the gffdb file to annotate precursor ids!")
    if custom_file is not None:
        try:
            with open(custom_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):  # comment chars
                        gene_id, gene_name = line.rstrip(
                            '\n').split('\t')  # ignore other columns
                        name2id[gene_name].append(gene_id)
        except Exception as e:
            print(e)
            print("Could not use the custom file to annotate precursor ids!")
    return name2id


def read_and_append_id(in_file, sep, name_col, gffdb_file, custom_file, delete_original_column):
    """Read a delimited file and append matching gene IDs.

    Args:
        in_file (str): Path to the input delimited file.
        sep (str): Field delimiter.
        name_col (str): Column name containing gene names.
        gffdb_file (str or None): Path to a gffutils database file.
        custom_file (str or None): Path to a custom id-to-name mapping file.
        delete_original_column (bool): If True, remove ``name_col`` and index
            by ``gene id``; otherwise keep ``name_col`` as the index.

    Returns:
        pandas.DataFrame: Table with gene IDs appended and index set as
        described above.
    """
    pd = _require_pandas()
    df = pd.read_table(in_file, sep=sep)

    name2id_dict = build_name_to_id_dict(gffdb_file, custom_file)
    df = add_id_to_dataframe(df, name2id_dict, name_col)

    if delete_original_column:
        del df[name_col]
        df.set_index('gene id', inplace=True)
    else:
        df.set_index(name_col, inplace=True)
    return df


def convert(in_file, sep, name_col, out_file, gtfdb_file, custom_file, delete_original_column):
    """Read an input file, append gene IDs, and write the result.

    Args:
        in_file (str): Path to the input delimited file.
        sep (str): Field delimiter.
        name_col (str): Column name containing gene names.
        out_file (str): Path to the output file.
        gtfdb_file (str or None): Path to a gffutils database file.
        custom_file (str or None): Path to a custom id-to-name mapping file.
        delete_original_column (bool): If True, drop the original name column
            before writing.
    """
    df = read_and_append_id(
        in_file, sep, name_col, gtfdb_file, custom_file,
        delete_original_column
    )
    df.to_csv(out_file, sep=str(sep), index=True, header=True)


def main():
    """Parse CLI arguments and run gene-name-to-ID conversion."""
    # Setup argument parser
    parser = ArgumentParser()

    parser.add_argument(
        "--output",
        dest="output",
        help="output file",
        required=False
    )
    parser.add_argument(
        "--input",
        dest="input",
        help="input file",
        required=True
    )
    parser.add_argument(
        "--sep",
        dest="sep",
        help="delimiter, default is tabbed (\'\\t\'), but can be specified here",
        required=False,
        default='\t'
    )
    parser.add_argument(
        "--name_col",
        dest="name_col",
        help="name of the column that contains the miRNA name",
        required=False,
        default=None
    )
    parser.add_argument(
        "--gffdb",
        dest="gffdb",
        help="gff database file create from gffutils (Currently only "
             "GFF files are allowed)",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--custom",
        dest="custom",
        help="custom tabbed file (id\tname) for genes not in the gff database",
        required=False,
        default=None
    )
    parser.add_argument(
        "--delete_original_column",
        dest="delete_original_column",
        help="if True, delete the original column (name_col)",
        required=False,
        action='store_true',
        default=False
    )

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    in_file = args.input
    sep = args.sep
    name_col = args.name_col
    out_file = args.output
    gffdb_file = args.gffdb
    custom_file = args.custom
    delete_original_column = args.delete_original_column

    convert(in_file, sep, name_col, out_file, gffdb_file,
            custom_file, delete_original_column)


if __name__ == "__main__":
    main()
