#!/usr/bin/env python
"""Calculate average and total genomic region lengths from a gffutils database.

Reports CDS, UTR, and intron lengths (total and average) for use in metagene
plots and region-distribution summaries.
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
from . import annotation_functions as af
from . import create_region_bedfiles
from collections import defaultdict


def calculate_total_cds_length(db, keys):
    """Calculate the total non-overlapping CDS length across all transcripts.

    Args:
        db (gffutils.FeatureDB): Annotation database.
        keys (dict): Species-specific GTF column naming scheme; maps logical
            region names (e.g. ``'cds'``) to feature type strings.

    Returns:
        int: Total CDS coverage in base pairs.
    """

    ### This method counts every CDS feature, may double count overlapping features ###

    # total_cds_feature_tx_lengths = 0
    # for cds_feature in db.features_of_type(chr19_keys['cds']):
    #     total_cds_feature_tx_lengths += (cds_feature.end - cds_feature.start)
    # print('cds lengths (total): {}'.format(total_cds_feature_tx_lengths))

    ### This method counts only nonoverlapping CDS regions ###
    cds = create_region_bedfiles.create_cds_region_bedfile(
        db, keys, None, True)
    cds = cds.sort()
    print('cds lengths (total): {}'.format(cds.total_coverage()))
    return cds.total_coverage()


def calculate_avg_cds_length(db, keys):
    """Calculate the average total CDS length per transcript.

    Args:
        db (gffutils.FeatureDB): Annotation database.
        keys (dict): Species-specific GTF column naming scheme.

    Returns:
        float: Mean CDS length across transcripts, in base pairs.
    """

    cds_feature_tx_lengths = defaultdict(list)

    # get all cds regions that belong to each transcript
    for cds_feature in db.features_of_type(keys['cds']):
        for transcript in cds_feature.attributes[keys['transcript_id']]:
            cds_feature_tx_lengths[transcript].append(
                cds_feature.end - cds_feature.start)

    # calculate for each transcript, the total cds length
    cds_features = []
    for transcript in cds_feature_tx_lengths.keys():
        cds_features.append(sum(cds_feature_tx_lengths[transcript]))
    print('cds lengths (avg): {}'.format(
        sum(cds_features) / float(len(cds_features))))
    return sum(cds_features) / float(len(cds_features))


def calculate_total_utr_lengths(db, cds_dict, keys):
    """Calculate total non-overlapping 5' and 3' UTR lengths.

    Args:
        db (gffutils.FeatureDB): Annotation database.
        cds_dict (collections.defaultdict): Mapping
            ``{transcript_id: {'start': int, 'end': int}}`` of CDS bounds.
        keys (dict): Species-specific GTF column naming scheme.

    Returns:
        tuple: ``(total_five_prime_utr_length, total_three_prime_utr_length)``
            in base pairs from merged non-overlapping UTR BED coverage.
    """

    utr5_out = None  # no need to write to file
    utr3_out = None  # no need to write to file
    by_transcript = True  # dont merge per gene, keep transcript separate

    utr5, utr3 = create_region_bedfiles.create_utr_region_bedfiles(
        db, keys, cds_dict, utr3_out, utr5_out, by_transcript
    )
    utr5 = utr5.sort()
    utr3 = utr3.sort()

    total_five_prime_utr_tx_lengths = utr5.total_coverage()
    total_three_prime_utr_tx_lengths = utr3.total_coverage()
    print('five_prime_utr_lengths (total): {}'.format(total_five_prime_utr_tx_lengths))
    print('three_prime_utr_lengths (total): {}'.format(total_three_prime_utr_tx_lengths))
    return total_five_prime_utr_tx_lengths, total_three_prime_utr_tx_lengths


def calculate_avg_utr_lengths(db, cds_dict, keys):
    """Calculate average 5' and 3' UTR lengths per transcript.

    Args:
        db (gffutils.FeatureDB): Annotation database.
        cds_dict (collections.defaultdict): Mapping
            ``{transcript_id: {'start': int, 'end': int}}`` of CDS bounds.
        keys (dict): Species-specific GTF column naming scheme.

    Returns:
        None: Results are printed to stdout.
    """

    # {transcript_id: [len(tx1), len(tx2), ...]
    five_prime_utr_tx_lengths = defaultdict(list)
    three_prime_utr_tx_lengths = defaultdict(list)

    # Use CDS dict to determine length for all UTRs, whether it's 3/5'
    for utr_feature in db.features_of_type(keys['utr']):
        classified_utr = af.classify_utr(
            utr_feature, cds_dict, keys['transcript_id']
        )
        if classified_utr == '5utr':
            for transcript in utr_feature.attributes[keys['transcript_id']]:
                five_prime_utr_tx_lengths[transcript].append(
                    utr_feature.end - utr_feature.start
                )
        elif classified_utr == '3utr':
            for transcript in utr_feature.attributes[keys['transcript_id']]:
                three_prime_utr_tx_lengths[transcript].append(
                    utr_feature.end - utr_feature.start
                )

    # list of sums of each transcript utr length
    five_prime_utr_lengths = []
    three_prime_utr_lengths = []

    for transcript in three_prime_utr_tx_lengths.keys():
        three_prime_utr_lengths.append(
            sum(three_prime_utr_tx_lengths[transcript])
        )
    for transcript in five_prime_utr_tx_lengths.keys():
        five_prime_utr_lengths.append(
            sum(five_prime_utr_tx_lengths[transcript])
        )

    # divide by total length by number of transcripts
    avg_three_prime_utr_length = sum(
        three_prime_utr_lengths) / float(len(three_prime_utr_lengths))
    avg_five_prime_utr_length = sum(
        five_prime_utr_lengths) / float(len(five_prime_utr_lengths))

    print(
        'three_prime_utr_lengths (avg): {}'.format(
            avg_three_prime_utr_length
        )
    )
    print(
        'five_prime_utr_lengths (avg): {}'.format(
            avg_five_prime_utr_length
        )
    )


def calculate_total_intron_lengths(db, exons_dict, transcripts_dict, keys):
    """Calculate total proximal, distal, and combined intron lengths.

    Args:
        db (gffutils.FeatureDB): Annotation database.
        exons_dict: Transcript-to-exon mapping from
            ``annotation_functions.get_all_exons_dict``.
        transcripts_dict: Transcript feature mapping from
            ``annotation_functions.get_all_transcripts_dict``.
        keys (dict): Species-specific GTF column naming scheme.

    Returns:
        None: Results are printed to stdout.
    """

    proxintron_out = None  # don't write to any file
    distintron_out = None  # don't write to any file
    allintron_out = None  # don't write to any file
    by_transcript = True  # don't merge by gene (keep transcripts separate)

    proxintrons, distintrons, allintrons = create_region_bedfiles.create_intron_region_bedfiles(
        db, exons_dict, transcripts_dict, keys,
        proxintron_out, distintron_out, allintron_out, by_transcript
    )
    proxintrons = proxintrons.sort()
    distintrons = distintrons.sort()
    allintrons = allintrons.sort()

    print('proximal intron lengths (total): {}'.format(
        proxintrons.total_coverage()))
    print('distal intron lengths (total): {}'.format(
        distintrons.total_coverage()))
    print('all intron lengths (total): {}'.format(allintrons.total_coverage()))


def calculate_avg_lengths(db_file, species):
    """Print average CDS and UTR lengths for a gffutils database.

    Can be used to generate ratios for metagene plots.

    Args:
        db_file (str): Path to a gffutils database file.
        species (str): Species identifier used to select GTF column nomenclature.

    Returns:
        None: Results are printed to stdout.
    """

    keys = af.get_keys(species)
    db = gffutils.FeatureDB(db_file)

    cds_dict = af.get_all_cds_dict(db, keys['cds'], keys['transcript_id'])

    calculate_avg_cds_length(db, keys)
    calculate_avg_utr_lengths(db, cds_dict, keys)


def calculate_total_lengths(db_file, species):
    """Print total CDS, UTR, and intron base counts for a gffutils database.

    Can be used to generate ratios for region-distribution pie charts.

    Args:
        db_file (str): Path to a gffutils database file.
        species (str): Species identifier used to select GTF column nomenclature.

    Returns:
        None: Results are printed to stdout.
    """
    keys = af.get_keys(species)
    db = gffutils.FeatureDB(db_file)

    cds_dict = af.get_all_cds_dict(db, keys['cds'], keys['transcript_id'])
    exons_dict = af.get_all_exons_dict(
        db,
        exon_key=keys['exon'],
        transcript_id_key=keys['transcript_id']
    )
    transcripts_dict = af.get_all_transcripts_dict(
        db,
        transcript_key=keys['transcript'],
        transcript_id_key=keys['transcript_id']
    )

    calculate_total_cds_length(db, keys)
    calculate_total_utr_lengths(db, cds_dict, keys)
    calculate_total_intron_lengths(db, exons_dict, transcripts_dict, keys)


def main():
    """Parse CLI arguments and print average and total region lengths."""
    # Setup argument parser
    parser = ArgumentParser()

    parser.add_argument(
        "--db_file",
        dest="db_file",
        help="gtfdb file",
        required=True
    )
    parser.add_argument(
        "--species",
        dest="species",
        help="sets the species gtf nomenclature",
        required=True
    )
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    db_file = args.db_file
    species = args.species

    calculate_avg_lengths(db_file, species)
    calculate_total_lengths(db_file, species)


if __name__ == "__main__":
    main()
