"""Annotation helper functions for GTF-based transcript and gene analysis.

Provides utilities for mapping GTF feature keys across annotation formats,
building transcript and gene dictionaries from gffutils databases, computing
gene and mRNA lengths, classifying UTR regions, and deriving intron coordinates.
"""

import gffutils
import sys
from collections import defaultdict

from .intervals import Interval

try:
    import pybedtools
except ImportError:
    pybedtools = None
# from future.utils import iteritems

MAXVAL = 1000000000
MINVAL = 0


def get_keys(gtf_format):
    """Return feature-type key names for a given GTF annotation format.

    Different annotation sources use different feature and attribute names in
    their GTF files. Defaults to the GENCODE standard used for hg19, mm9, mm10,
    and hg38 when ``gtf_format`` is unrecognized.

    Args:
        gtf_format: Annotation source identifier. One of ``'ensembl'``,
            ``'ce10'``, ``'ce11'``, ``'wormbase'``, ``'dm6_genbank'``,
            ``'genbank'``, ``'refseq'``, ``'ucsc'``, or any other value
            (which defaults to GENCODE).

    Returns:
        A dictionary mapping logical feature names to GTF feature/attribute
        keys. Keys include ``'cds'``, ``'utr3'``, ``'utr5'``, ``'utr'``,
        ``'gene'``, ``'gene_name'``, ``'transcript'``, ``'transcript_id'``,
        ``'transcript_type'``, ``'exon'``, ``'gene_id'``, and ``'gene_type'``.
    """
    if gtf_format == 'ensembl':
        cds_key = 'CDS'
        utr3_key = 'three_prime_utr'
        utr5_key = 'five_prime_utr'
        utr_key = None
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_biotype'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_biotype'
    elif gtf_format == 'ce10':
        cds_key = 'CDS'
        utr3_key = None  # 'three_prime_UTR'
        utr5_key = None  # 'five_prime_UTR'
        utr_key = 'UTR'
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_biotype'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_biotype'
    elif gtf_format == 'ce11' or gtf_format == 'wormbase':
        cds_key = 'CDS'
        utr3_key = 'three_prime_utr'
        utr5_key = 'five_prime_utr'
        utr_key = None
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_biotype'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_biotype'
    elif gtf_format == 'dm6_genbank' or gtf_format == 'genbank':
        cds_key = 'Genbank CDS'
        utr3_key = None
        utr5_key = None
        utr_key = None
        gene_key = 'Genbank gene'
        gene_name_key = 'gene'
        transcript_key = 'transcript_id'
        transcript_id_key = 'transcript_id'
        type_key = 'gene_biotype'
        exon_key = 'Genbank exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_biotype'
    elif gtf_format == 'refseq':
        cds_key = 'CDS'
        utr3_key = None
        utr5_key = None
        utr_key = None
        gene_key = 'gene'
        gene_name_key = 'gene'
        transcript_key = 'ID'
        transcript_id_key = 'ID'
        type_key = 'gene_biotype'
        exon_key = 'exon'
        gene_id_key = 'GeneID'
        gene_type_key = 'gene_biotype'
    elif gtf_format == 'ucsc':
        cds_key = 'CDS'
        utr3_key = '3UTR'
        utr5_key = '5UTR'
        utr_key = None
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'ID'
        transcript_id_key = 'transcript_id'
        type_key = 'gene_biotype'
        exon_key = 'exon'
        gene_id_key = 'GeneID'
        gene_type_key = 'gene_type'
    else:  # default to gencode
        cds_key = 'CDS'
        utr3_key = None  #
        utr5_key = None  # in human/mice, this key doesn't exist
        utr_key = 'UTR'
        gene_key = 'gene'
        gene_name_key = 'gene_name'
        transcript_key = 'transcript'
        transcript_id_key = 'transcript_id'
        type_key = 'transcript_type'
        exon_key = 'exon'
        gene_id_key = 'gene_id'
        gene_type_key = 'gene_type'

    keys = {
        'cds': cds_key,
        'utr3': utr3_key,
        'utr5': utr5_key,
        'utr': utr_key,
        'gene': gene_key,
        'gene_name': gene_name_key,
        'transcript': transcript_key,
        'transcript_id': transcript_id_key,
        'transcript_type': type_key,
        'exon': exon_key,
        'gene_id': gene_id_key,
        'gene_type': gene_type_key
    }
    return keys


def get_gene_to_transcript_dict(db, gene_id_key, transcript_id_key):
    """Build a mapping from gene IDs to their transcript IDs.

    Args:
        db: gffutils FeatureDB containing transcript features.
        gene_id_key: Attribute name for the gene identifier in the GTF.
        transcript_id_key: Attribute name for the transcript identifier in
            the GTF.

    Returns:
        A ``defaultdict(list)`` mapping each gene ID to a list of associated
        transcript IDs, e.g. ``{gene_id: [transcript1, transcript2, ...]}``.
    """

    genes_dict = defaultdict(list)
    for gene_feature in db.features_of_type('transcript'):
        for gene_id in gene_feature.attributes[gene_id_key]:
            for transcript_id in gene_feature.attributes[transcript_id_key]:
                genes_dict[gene_id].append(transcript_id)
    return genes_dict


def get_longest_transcripts(genes_dict, transcripts_dict):
    """Select the longest transcript for each gene by genomic span.

    Args:
        genes_dict: Mapping from gene ID to a list of transcript IDs.
        transcripts_dict: Mapping from transcript ID to a dict with ``'start'``
            and ``'end'`` coordinates.

    Returns:
        A ``defaultdict(dict)`` mapping each gene ID to the transcript ID with
        the largest ``end - start`` span.
    """
    longest_genes_dict = defaultdict(dict)
    for gene, transcripts in genes_dict.items():
        max_transcript_len = -1
        max_transcript = ""
        for transcript in transcripts:
            transcript_len = transcripts_dict[transcript]['end'] - \
                transcripts_dict[transcript]['start']
            if transcript_len > max_transcript_len:
                max_transcript_len = transcript_len
                max_transcript = transcript
        longest_genes_dict[gene] = max_transcript
    return longest_genes_dict


def most_upstream_downstream_positions(genes_dict, transcripts_dict):
    """Compute the genomic span covering all transcripts of each gene.

    For each gene, finds the minimum transcript start and maximum transcript
    end across all associated isoforms.

    Args:
        genes_dict: Mapping from gene ID to a list of transcript IDs.
        transcripts_dict: Mapping from transcript ID to a dict with ``'start'``
            and ``'end'`` coordinates.

    Returns:
        A ``defaultdict(dict)`` mapping each gene ID to
        ``{'start': min_start, 'end': max_end}``.
    """
    d = defaultdict(dict)
    for gene, transcripts in genes_dict.items():
        # as long as we don't have any chromosomes larger than 1 billion... no support for axolotls!
        min_transcript_pos = 1000000000
        max_transcript_pos = -1
        for transcript in transcripts:
            if transcripts_dict[transcript]['end'] > max_transcript_pos:
                max_transcript_pos = transcripts_dict[transcript]['end']
            if transcripts_dict[transcript]['start'] < min_transcript_pos:
                min_transcript_pos = transcripts_dict[transcript]['start']

        d[gene] = {'start': min_transcript_pos, 'end': max_transcript_pos}
    return d


def get_premrna_lengths(gene_id, upstream_downstream):
    """Compute the pre-mRNA length for a gene from its genomic span.

    Args:
        gene_id: Gene identifier to look up in ``upstream_downstream``.
        upstream_downstream: Mapping from gene ID to a dict with ``'start'``
            and ``'end'`` coordinates (1-based inclusive).

    Returns:
        The pre-mRNA length in base pairs: ``end - start + 1``.
    """
    start = upstream_downstream[gene_id]['start']
    end = upstream_downstream[gene_id]['end']
    return end - start + 1    # (coords are 1-based inclusive still)


def get_mrna_lengths(gene_id, exons_dict, genes_dict):
    """Compute the merged mRNA length for a gene across all isoforms.

    Collects exons from every transcript of the gene, merges overlapping
    intervals (not strand-specific), and sums the merged exonic length.

    Args:
        gene_id: Gene identifier to look up in ``genes_dict``.
        exons_dict: Mapping from transcript ID to a list of exon dicts, each
            with ``'chrom'``, ``'start'``, ``'end'``, and ``'strand'``.
        genes_dict: Mapping from gene ID to a list of transcript IDs.

    Returns:
        Total merged exonic length in base pairs for the gene.

    Raises:
        ImportError: If ``pybedtools`` is not installed.
    """
    if pybedtools is None:
        raise ImportError("pybedtools is required for get_mrna_lengths")
    exons_list = []
    total_mrna_length = 0
    for transcript in genes_dict[gene_id]:
        for exon in exons_dict[transcript]:
            exons_list.append(pybedtools.create_interval_from_list(
                [exon['chrom'],
                 # since we're converting to bedtool, use 0-based
                 str(exon['start'] - 1),
                 str(exon['end']),
                 transcript,
                 '0',
                 exon['strand']]
            ))
    exons_list = pybedtools.BedTool(exons_list)
    exons_list = exons_list.sort().merge()  # not strand specific
    for exon in exons_list:
        total_mrna_length += (exon.end - exon.start)
    return total_mrna_length


def get_all_cds_dict(db, cds_key, transcript_id_key='transcript_id'):
    """Build a per-transcript CDS coordinate dictionary.

    For every CDS-annotated transcript, records the lowest start and highest
    end across all CDS features belonging to that transcript.

    Args:
        db: gffutils FeatureDB containing CDS features.
        cds_key: GTF feature type name for CDS entries (format-specific).
        transcript_id_key: Attribute name for the transcript identifier in
            the GTF. Defaults to ``'transcript_id'`` (GENCODE).

    Returns:
        A ``defaultdict(dict)`` mapping each transcript ID to
        ``{'low': min_start, 'hi': max_end}``.
    """
    # chr19_cds_dict = defaultdict(lambda: {'low': MAXVAL, 'hi': MINVAL})
    cds_dict = defaultdict(dict)
    for cds_feature in db.features_of_type(cds_key):
        for transcript_id in cds_feature.attributes.get(transcript_id_key, []):
            # if cds_feature.start <= chr19_cds_dict[transcript_id]['low']:
            if cds_feature.start <= cds_dict[transcript_id].get("low", MAXVAL):
                cds_dict[transcript_id]['low'] = cds_feature.start
            # if cds_feature.end >= chr19_cds_dict[transcript_id]['hi']:
            if cds_feature.end >= cds_dict[transcript_id].get("hi", MINVAL):
                cds_dict[transcript_id]['hi'] = cds_feature.end
    return cds_dict


def get_all_exons_dict(db, exon_key, transcript_id_key):
    """Build a per-transcript exon coordinate dictionary.

    Args:
        db: gffutils FeatureDB containing exon features.
        exon_key: GTF feature type name for exon entries (format-specific).
        transcript_id_key: Attribute name for the transcript identifier in
            the GTF.

    Returns:
        A ``defaultdict(list)`` mapping each transcript ID to a list of exon
        dicts, each with ``'chrom'``, ``'start'``, ``'end'``, and ``'strand'``.
    """

    exons_dict = defaultdict(list)
    for exon_feature in db.features_of_type(exon_key):
        for transcript_id in exon_feature.attributes[transcript_id_key]:
            exons_dict[transcript_id].append(
                {
                    'chrom': exon_feature.seqid,
                    'start': exon_feature.start,
                    'end': exon_feature.end,
                    'strand': exon_feature.strand,
                }
            )
    return exons_dict


def get_all_transcripts_dict(db, transcript_key, transcript_id_key):
    """Build a per-transcript genomic coordinate dictionary.

    Args:
        db: gffutils FeatureDB containing transcript features.
        transcript_key: GTF feature type name for transcript entries
            (format-specific).
        transcript_id_key: Attribute name for the transcript identifier in
            the GTF.

    Returns:
        A ``defaultdict(dict)`` mapping each transcript ID to
        ``{'start': start, 'end': end}``.
    """
    transcripts_dict = defaultdict(dict)
    for transcript_feature in db.features_of_type(transcript_key):
        for transcript_id in transcript_feature.attributes[transcript_id_key]:
            transcripts_dict[transcript_id] = {
                'start': transcript_feature.start,
                'end': transcript_feature.end
            }
    return transcripts_dict


def chromosome_set(db):
    """Return the set of chromosome (seqid) names in a FeatureDB.

    Args:
        db: gffutils FeatureDB to query.

    Returns:
        A ``set`` of unique seqid strings present in the database.
    """
    ret = db.execute("SELECT seqid FROM features").fetchall()
    all_chromosomes = [r['seqid'] for r in ret]
    return set(all_chromosomes)


def gene_id_to_name(db, gene_name_key, gene_id_key='gene_id', strict=False):
    """Build a gene ID to gene name mapping from a FeatureDB.

    Adapted from ``gscripts.region_helpers``. Attribute key names may differ
    if ``gene_id_key`` or the supplied ``gene_name_key`` are absent from the
    source GTF.

    Args:
        db: gffutils FeatureDB created via ``gffutils.create_db()``.
        gene_name_key: Attribute name for the gene display name, typically
            ``'gene_name'``.
        gene_id_key: Attribute name for the gene identifier in the GTF.
        strict: If True, raise ``KeyError`` when required attributes are
            missing. If False, skip genes with missing keys after printing
            a warning.

    Returns:
        A dict mapping each gene ID to its gene name.

    Raises:
        KeyError: If ``strict`` is True and a gene lacks required attributes.
    """

    genes = db.features_of_type('gene')
    gene_name_dict = {}

    for gene in genes:
        try:
            gene_id = gene.attributes[gene_id_key][0] if type(
                gene.attributes[gene_id_key]
            ) == list else gene.attributes[gene_id_key]
            gene_name_dict[gene_id] = gene.attributes[gene_name_key][0]
        except KeyError:
            print(gene.attributes.keys())
            print("Warning. Key not found for {}".format(gene))
            if strict:
                raise
    return gene_name_dict


def classify_utr(utr_feature, cds_dict, transcript_id_key='transcript_id'):
    """Classify a UTR feature as 5' or 3' based on CDS positions.

    Determines whether a UTR lies upstream (5') or downstream (3') of the
    coding sequence, accounting for transcript strand.

    Args:
        utr_feature: gffutils Feature already classified as UTR (within a
            coding transcript but outside CDS regions).
        cds_dict: Mapping from transcript ID to CDS bounds with ``'low'``
            (minimum start) and ``'hi'`` (maximum end).
        transcript_id_key: Attribute name for the transcript identifier in
            the GTF. Defaults to ``'transcript_id'`` (GENCODE).

    Returns:
        One of ``'5utr'``, ``'3utr'``, or ``'unclassified_utr'``. Returns
        ``'unclassified_utr'`` when no CDS is found for the transcript.
    """
    three_prime_utr = False
    five_prime_utr = False

    for transcript_id in utr_feature.attributes.get(transcript_id_key, []):
        try:
            if utr_feature.strand == '+':
                if cds_dict[transcript_id]['low'] > utr_feature.end:
                    five_prime_utr = True
                if cds_dict[transcript_id]['hi'] < utr_feature.start + 1:
                    three_prime_utr = True
            elif utr_feature.strand == '-':
                if cds_dict[transcript_id]['low'] > utr_feature.end:
                    three_prime_utr = True
                if cds_dict[transcript_id]['hi'] < utr_feature.start + 1:
                    five_prime_utr = True
        except KeyError as e:
            # no CDS found for this transcript so cannot positionally assign utr
            return 'unclassified_utr'
    if five_prime_utr:
        return '5utr'
        # return 'five_prime_utr'
    elif three_prime_utr:
        return '3utr'
        # return 'three_prime_utr'
    else:
        return 'unclassified_utr'


def find_introns(transcript, exons):
    """Derive intron coordinates from transcript bounds and exon positions.

    Computes the gaps between exons within a transcript. Coordinates are
    1-based inclusive, consistent with the exon input format.

    Args:
        transcript: Dict with ``'start'`` and ``'end'`` defining the transcript
            span.
        exons: List of exon dicts, each with ``'start'`` and ``'end'``.

    Returns:
        A list of intron dicts, each with ``'start'`` and ``'end'``. Returns
        an empty list if no exons are found.
    """
    positions = []  # list of
    introns = []

    pos_append = positions.append
    intron_append = introns.append

    for exon in exons:
        pos_append(exon['start'] - 1)
        pos_append(exon['end'] + 1)
    positions = sorted(positions)
    try:
        # feature doesn't start with an intron
        if positions[0] < transcript['start']:
            positions.pop(0)
        else:
            positions.insert(transcript['start'])
        # feature doesn't end at an intron
        if positions[-1] > transcript['end']:
            positions.pop(-1)
        else:
            pos_append(transcript['end'])
        for i in range(0, len(positions) - 1, 2):
            intron_append({'start': positions[i], 'end': positions[i + 1]})
    except IndexError:
        print("No exons found: {}".format(transcript), positions)
        return []
    return introns


def get_proxdist_from_intron(interval, distance=500):
    """Split an intron interval into proximal and distal sub-intervals.

    Divides an intron into flanking proximal regions (``distance`` bp at each
    end) and a central distal region. Intervals use 0-based half-open
    coordinates matching BED format.

    Args:
        interval: Intron ``Interval`` with ``chrom``, ``start``, ``end``,
            ``score``, and ``strand`` attributes.
        distance: Size in base pairs of each proximal flank. Defaults to 500.

    Returns:
        A dict with keys ``'prox'`` and ``'dist'``, each mapping to a list of
        ``Interval`` objects. If the intron is shorter than ``2 * distance``,
        the entire interval is returned as a single proximal region and the
        distal list is empty.
    """
    d = distance * 2
    score = str(interval.score)
    prox_name = 'proxintron{}'.format(distance)
    dist_name = 'distintron{}'.format(distance)

    if interval.end - interval.start <= d:
        prox = Interval(
            interval.chrom,
            interval.start,
            interval.end,
            prox_name,
            score,
            interval.strand,
        )
        return {'prox': [prox], 'dist': []}

    prox_left = Interval(
        interval.chrom,
        interval.start,
        interval.start + distance,
        prox_name,
        score,
        interval.strand,
    )
    prox_right = Interval(
        interval.chrom,
        interval.end - distance,
        interval.end,
        prox_name,
        score,
        interval.strand,
    )
    dist = Interval(
        interval.chrom,
        interval.start + distance,
        interval.end - distance,
        dist_name,
        score,
        interval.strand,
    )
    return {'prox': [prox_left, prox_right], 'dist': [dist]}
