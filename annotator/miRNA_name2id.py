#!/usr/bin/env python
"""Convert tabular files labeled with miRNA names to miRNA accessions.

Uses a gffutils database and/or a custom mapping file to look up precursor
(and optionally mature) miRNA IDs and append them to an input table.
"""

# these two are really a minimum

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
            "pandas is required for miRNA_name2id. "
            "Install with: pip install 'annotator[utils]'"
        ) from exc
    return pd


def get_precursor_id_from_namerow(row, name2id, name_col):
    """Look up precursor miRNA ID(s) for a single row's miRNA name.

    Args:
        row (pandas.Series): Row containing a miRNA name column.
        name2id (dict): Mapping from miRNA name to list of precursor IDs.
        name_col (str): Column name holding the miRNA name.

    Returns:
        str: Pipe-delimited precursor ID(s), or an empty string if missing.
    """
    try:
        return '|'.join(name2id[row[name_col]])
    except KeyError:
        print('Cant find {} in the database!'.format(row[name_col]))
        return ''


def get_mature_id_from_idrow(row, precursor2mature, name_col):
    """Look up mature miRNA ID(s) for a single row's precursor ID.

    Args:
        row (pandas.Series): Row containing a precursor ID column.
        precursor2mature (dict): Mapping from precursor ID to list of mature IDs.
        name_col (str): Column name holding the precursor ID.

    Returns:
        str: Pipe-delimited mature ID(s), or an empty string if missing.
    """
    try:
        return '|'.join(precursor2mature[row[name_col]])
    except KeyError:
        print('Cant find {} in the database!'.format(row[name_col]))
        return ''


def add_precursor_id_to_dataframe(df, name2id_dict, name_col):
    """Append a ``miRNA id`` column with precursor IDs to a dataframe.

    Args:
        df (pandas.DataFrame): Input table with miRNA names.
        name2id_dict (dict): Mapping from miRNA name to list of precursor IDs.
        name_col (str): Column name holding the miRNA name.

    Returns:
        pandas.DataFrame: Copy of ``df`` with a new ``miRNA id`` column.
    """
    assert 'miRNA id' not in df.columns
    df['miRNA id'] = df.apply(get_precursor_id_from_namerow, args=[
                              name2id_dict, name_col], axis=1)
    return df


def add_mature_id_to_dataframe(df, precursor2mature_dict, name_col):
    """Append a ``mature`` column with mature miRNA IDs to a dataframe.

    Args:
        df (pandas.DataFrame): Input table with precursor IDs.
        precursor2mature_dict (dict): Mapping from precursor ID to list of
            mature IDs.
        name_col (str): Column name holding the precursor ID.

    Returns:
        pandas.DataFrame: Copy of ``df`` with a new ``mature`` column.
    """
    df['mature'] = df.apply(get_mature_id_from_idrow, args=[
                            precursor2mature_dict, name_col], axis=1)
    return df


def build_precursor_id_to_mature_id_dict(gffdb_file, custom_file):
    """Build a precursor-ID-to-mature-ID mapping from a gffutils DB and/or custom file.

    Args:
        gffdb_file (str or None): Path to a gffutils database containing
            ``miRNA`` features with ``ID`` and ``Derives_from`` attributes.
        custom_file (str or None): Tab-delimited file with precursor ID,
            unused column, and pipe-delimited mature IDs.

    Returns:
        dict: Mapping ``{precursor_id: [mature_id, ...]}``.
    """
    precursor2mature = defaultdict(list)
    if gffdb_file is not None:
        try:
            db = gffutils.FeatureDB(gffdb_file)
            for transcript in db.features_of_type('miRNA'):
                for mature_id in transcript.attributes['ID']:
                    for parent in transcript.attributes['Derives_from']:
                        precursor2mature[parent].append(mature_id)
        except Exception as e:
            print(e)
            print("Could not use the gffdb file to annotate mature ids!")
    if custom_file is not None:
        try:
            with open(custom_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):  # comment chars
                        id, _, mature_ids = line.rstrip('\n').split('\t')
                        mature_ids = mature_ids.split('|')
                        if id in precursor2mature:
                            precursor2mature[id].extend(mature_ids)
                        else:
                            precursor2mature[id] = list(mature_ids)
        except Exception as e:
            print(e)
            print("Could not use the custom file to annotate mature ids!")

    return precursor2mature


def build_name_to_precursor_id_dict(gffdb_file, custom_file):
    """Build a miRNA-name-to-precursor-ID mapping from a gffutils DB and/or custom file.

    Names with more than one associated ID are stored as a list; callers join
    multiple IDs with ``|``.

    Args:
        gffdb_file (str or None): Path to a gffutils database. The database
            must contain ``miRNA_primary_transcript`` features with ``Name``
            and ``ID`` attributes.
        custom_file (str or None): Tab-delimited file of ``id<TAB>name`` rows
            for miRNAs not present in the gffutils database.

    Returns:
        dict: Mapping ``{name: [id, ...]}``.
    """
    # delim = '|'
    precursor_name2id = defaultdict(list)
    if gffdb_file is not None:
        try:
            db = gffutils.FeatureDB(gffdb_file)
            for transcript in db.features_of_type('miRNA_primary_transcript'):
                for name in transcript.attributes['Name']:
                    precursor_name2id[name].append(transcript.id)

        except Exception as e:
            print(e)
            print("Could not use the gffdb file to annotate precursor ids!")
    if custom_file is not None:
        try:
            with open(custom_file, 'r') as f:
                for line in f:
                    if not line.startswith('#'):  # comment chars
                        id, name = line.rstrip('\n').split(
                            '\t')[:2]  # ignore other columns
                        precursor_name2id[name].append(id)
        except Exception as e:
            print(e)
            print("Could not use the custom file to annotate precursor ids!")
    return precursor_name2id


def read_and_append_id(in_file, sep, name_col, gffdb_file, custom_file, add_mature, delete_original_column):
    """Read a delimited file and append matching miRNA precursor IDs.

    Args:
        in_file (str): Path to the input delimited file.
        sep (str): Field delimiter.
        name_col (str): Column name containing miRNA names.
        gffdb_file (str or None): Path to a gffutils database file.
        custom_file (str or None): Path to a custom id-to-name mapping file.
        add_mature (bool): If True, also append a ``mature`` column with mature
            miRNA IDs.
        delete_original_column (bool): If True, remove ``name_col`` and index
            by ``miRNA id``; otherwise keep ``name_col`` as the index.

    Returns:
        pandas.DataFrame: Table with miRNA IDs appended and index set as
        described above.
    """
    pd = _require_pandas()
    df = pd.read_table(in_file, sep=sep)

    precursor_name2id_dict = build_name_to_precursor_id_dict(
        gffdb_file, custom_file)
    df = add_precursor_id_to_dataframe(df, precursor_name2id_dict, name_col)

    if add_mature:
        precursor2mature_dict = build_precursor_id_to_mature_id_dict(
            gffdb_file, custom_file)
        df = add_mature_id_to_dataframe(df, precursor2mature_dict, 'miRNA id')

    if delete_original_column:
        del df[name_col]
        df.set_index('miRNA id', inplace=True)
    else:
        df.set_index(name_col, inplace=True)
    return df


def convert(in_file, sep, name_col, out_file, gtfdb_file, custom_file, add_mature, delete_original_column=True):
    """Read an input file, append miRNA IDs, and write the result.

    Args:
        in_file (str): Path to the input delimited file.
        sep (str): Field delimiter.
        name_col (str): Column name containing miRNA names.
        out_file (str): Path to the output file.
        gtfdb_file (str or None): Path to a gffutils database file.
        custom_file (str or None): Path to a custom id-to-name mapping file.
        add_mature (bool): If True, also append mature miRNA IDs.
        delete_original_column (bool): If True, drop the original name column
            before writing.
    """
    df = read_and_append_id(
        in_file, sep, name_col, gtfdb_file, custom_file,
        add_mature,
        delete_original_column
    )

    df.to_csv(out_file, sep=sep, index=True, header=True)


def main():
    """Parse CLI arguments and run miRNA-name-to-ID conversion."""
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
        help="default is tabbed (\'\\t\'), but can be specified here",
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
        "--add_mature",
        action='store_true',
        help="add the mature sequences, | delimited, as another column",
        required=False,
        default=False
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
        help="custom tabbed file (id\tname) for mirs not in the gff database",
        required=False,
        default=None
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
    add_mature = args.add_mature

    convert(in_file, sep, name_col, out_file,
            gffdb_file, custom_file, add_mature)


if __name__ == "__main__":
    main()
