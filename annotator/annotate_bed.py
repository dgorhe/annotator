#!/usr/bin/env python
"""BED file annotation engine.

Provides functions to annotate BED6 intervals against a gffutils annotation
database by spatial hashing, overlap detection, and priority-based
gene/transcript selection. Supports single-core and parallel annotation
workflows.
"""

import copy
import csv
import multiprocessing
import os
import sys
from collections import OrderedDict, defaultdict

import gffutils
from tqdm import trange

from . import annotation_functions as af
from .intervals import Interval
from .index_store import index_exists, load_annotator_index, index_is_current, index_path_for_db

HASH_VAL = 1000000
MAXVAL = 1000000000
MINVAL = 0
FEATURES_DICT = {}


def db_chrom_from_bed(chrom, append_chr):
    """Map a BED contig name to the gffutils seqid used under ``append_chr``.

    When ``append_chr`` is True the BED uses ``chr``-prefixed names while the
    database stores unprefixed seqids (e.g. ``chrI`` → ``I``).

    Args:
        chrom: Contig name from a BED interval.
        append_chr: Whether annotation is running with ``--append-chr``.

    Returns:
        Contig name suitable for ``db.region(seqid=...)`` / chrom limits.
    """
    if append_chr and chrom.startswith('chr'):
        return chrom[3:]
    return chrom


def bed_chrom_from_db(chrom, append_chr):
    """Map a gffutils seqid to the features_dict / index lookup key.

    Args:
        chrom: Database contig / seqid.
        append_chr: Whether annotation is running with ``--append-chr``.

    Returns:
        Contig name matching BED queries and index ``features.chrom`` values.
    """
    if append_chr and not chrom.startswith('chr'):
        return 'chr' + chrom
    return chrom


def _clone_feature_for_region(element, start, end, featuretype):
    """Shallow-copy a feature and set coordinates/type for prox/dist introns.

    Args:
        element: gffutils Feature to copy.
        start: 1-based start coordinate for the cloned feature.
        end: 1-based end coordinate for the cloned feature.
        featuretype: Feature type string (e.g. ``proxintron``, ``distintron``).

    Returns:
        A shallow copy of ``element`` with updated start, end, and featuretype.
    """
    feature = copy.copy(element)
    feature.start = start
    feature.end = end
    feature.featuretype = featuretype
    return feature


def create_definitions(db_file, chroms=None, gtf_format='gencode',
                       append_chr=False, fuzzy=0, use_index=True,
                       buckets=None, lazy_buckets=False,
                       max_cached_buckets=256):
    """Create definitions and populate dictionaries used for annotation.

    Previously initialized in ``Annotator()``. When a pre-built annotator
    index exists and ``use_index`` is True, loads from disk instead of
    scanning the database.

    Args:
        db_file: Path to the gffutils SQLite database file.
        chroms: List of chromosome names to restrict indexing. Empty list
            indexes all chromosomes in the database.
        gtf_format: Annotation syntax specifier. One of ``'gencode'``,
            ``'ensembl'``, or ``'refseq'``.
        append_chr: If True, prepend ``'chr'`` to chromosome names when
            building the feature hash.
        fuzzy: Fuzzy-matching parameter passed to feature hashing.
        use_index: If True and a pre-built index exists, load from it
            instead of rebuilding from the database.
        buckets: Optional iterable of ``(chrom, bucket)`` pairs limiting
            which spatial hash buckets are loaded from the index. Ignored
            on the legacy rebuild path.
        lazy_buckets: If True, load index buckets on demand with an LRU
            cache instead of hydrating ``buckets`` (or all chrom buckets)
            eagerly.
        max_cached_buckets: LRU size when ``lazy_buckets`` is True.

    Returns:
        Tuple of ``(exons_dict, transcripts_dict, cds_dict, features_dict,
        keys)``. When loaded from index, ``exons_dict`` and
        ``transcripts_dict`` are empty dicts.
    """

    if chroms is None:
        chroms = []

    keys = af.get_keys(gtf_format)

    if use_index and index_exists(db_file):
        index_file = index_path_for_db(db_file)
        if index_is_current(
                index_file, db_file, chroms=chroms, gtf_format=gtf_format,
                append_chr=append_chr):
            # Index stores BED-style contig keys (chr-prefixed when append_chr).
            index_chroms = None
            if chroms:
                index_chroms = [
                    bed_chrom_from_db(c, append_chr) for c in chroms
                ]
            cds_dict, features_dict, keys = load_annotator_index(
                index_file,
                chroms=index_chroms,
                gtf_format=gtf_format,
                buckets=buckets,
                lazy=lazy_buckets,
                max_cached_buckets=max_cached_buckets,
            )
            return {}, {}, cds_dict, features_dict, keys
        print(
            "Warning: annotator index is stale or incompatible; "
            "rebuilding from database.",
            file=sys.stderr,
        )

    db = gffutils.FeatureDB(db_file)
    featuretypes = [x.lower() for x in list(db.featuretypes())]

    progress = trange(4, desc='Initializing/creating defs')
    progress.set_description("Building gene ids -> gene names dictionary")

    progress.update(1)
    if len(chroms) != 0:  # if use specific chromosomes, otherwise hash all
        chromosomes = set(chroms)
    else:
        chromosomes = af.chromosome_set(db)

    exons_dict = {}  # holds all of the transcript_id:exon start/stop
    transcripts_dict = {}  # holds all of the transcript_id:start/stop
    cds_dict = {}  # holds all of the transcript_id:cds start/stop

    # If 'exon' features exist in db, catalog it.
    if 'exon' in featuretypes:
        progress.set_description("Adding all exon boundaries")
        exons_dict = af.get_all_exons_dict(
            db, keys['exon'], keys['transcript_id']
        )
    progress.update(1)
    # If 'transcript' features exist in db, catalog it.
    if 'transcript' in featuretypes:
        progress.set_description("Adding all transcript boundaries")
        transcripts_dict = af.get_all_transcripts_dict(
            db, keys['transcript'], keys['transcript_id']
        )
    progress.update(1)
    # If 'cds' features exist in db, catalog it.
    if 'cds' in featuretypes:
        progress.set_description("Adding all CDS boundaries")
        cds_dict = af.get_all_cds_dict(
            db, keys['cds'], keys['transcript_id']
        )
    progress.update(1)
    progress.set_description("Adding all features in DB to dictionary")
    features_dict, num_features = hash_features(
        db, chromosomes, append_chr, fuzzy,
        exons_dict, transcripts_dict,
        keys['transcript_id'], featuretypes
    )

    return exons_dict, transcripts_dict, cds_dict, features_dict, keys


def hash_features(db, chromosomes, append_chr, fuzzy,
                  exons_dict, transcripts_dict,
                  transcript_id_key, featuretypes):
    """Hash features by genomic position for fast overlap lookup.

    Builds a spatial index keyed by ``(chromosome, hash_bin, strand)`` where
    each bin spans ``HASH_VAL`` bases. Inferred introns and prox/dist intron
    regions are also stored when introns are not explicitly annotated.

    ``gene``-level GTF rows are skipped because the annotate path resolves
    overlaps through transcript-associated features only; indexing them added
    startup cost without changing output.

    Args:
        db: gffutils FeatureDB instance.
        chromosomes: Iterable of chromosome names to index.
        append_chr: If True, prepend ``'chr'`` to chromosome keys.
        fuzzy: Unused; retained for API compatibility.
        exons_dict: Mapping of transcript_id to exon boundary lists.
        transcripts_dict: Mapping of transcript_id to transcript boundaries.
        transcript_id_key: GTF attribute key for transcript identifiers.
        featuretypes: Lowercased list of feature types present in the database.

    Returns:
        Tuple of ``(features_dict, num_features)`` where ``features_dict`` is
        a ``defaultdict(list)`` mapping ``(chrom, pos/HASH_VAL, strand)`` to
        lists of gffutils Feature objects, and ``num_features`` is the count
        of features processed.
    """
    num_features = 0

    features_dict = defaultdict(list)

    progress = trange(
        len(chromosomes),
        leave=False,
        desc='Building feature index'
    )

    for chrom in chromosomes:
        fixed_chrom = 'chr' + chrom if append_chr else chrom

        progress.set_description("Adding {}...".format(chrom))

        for element in db.region(seqid=chrom):
            if element.featuretype == 'gene':
                continue

            introns = []
            # if we see at least one intron, we assume introns are annotated!
            if element.featuretype == 'transcript' and \
                    'intron' not in featuretypes:
                for transcript_id in element.attributes[transcript_id_key]:
                    if (transcript_id not in exons_dict or
                            transcript_id not in transcripts_dict):
                        continue
                    exons = exons_dict[transcript_id]
                    transcript = transcripts_dict[transcript_id]
                    # Find introns
                    introns = af.find_introns(transcript, exons)

            # Create multi-key hash and store elements
            start = int(element.start / HASH_VAL)
            end = int(element.end / HASH_VAL)

            for i in range(start, end + 1):
                append = features_dict[fixed_chrom, i, element.strand].append

                append(element)

                # append inferred introns to features dictionary
                for intron in introns:
                    intron_interval = Interval.from_list([
                        fixed_chrom,
                        str(int(intron['start']) - 1),
                        intron['end'],
                        'intron',
                        str(element.score),
                        element.strand,
                    ])
                    proxdist_dict = af.get_proxdist_from_intron(intron_interval)
                    for prox in proxdist_dict['prox']:
                        proxintron_feature = _clone_feature_for_region(
                            element,
                            prox.start + 1,
                            prox.end,
                            prox.name,
                        )
                        append(proxintron_feature)
                    for dist in proxdist_dict['dist']:
                        distintron_feature = _clone_feature_for_region(
                            element,
                            dist.start + 1,
                            dist.end,
                            dist.name,
                        )
                        append(distintron_feature)

            num_features += 1
        progress.update(1)
    # for feature in chr19_features_dict['chr7', 129, '+']:
    #     if feature.attributes['gene_id'][0] == 'ENSG00000128607.9':
    #         if feature.featuretype == 'CDS':
    #             print(feature)
    return features_dict, num_features


def get_all_overlapping_features_from_query_stranded(
        chrom, qstart, qend, strand, features_dict):
    """Return all same-strand features overlapping a query interval.

    Functions similarly to ``gffutils.FeatureDB.region()``, but uses the
    pre-hashed ``features_dict`` for faster lookup.

    Args:
        chrom: Chromosome or sequence identifier.
        qstart: Query start (0-based BED coordinate).
        qend: Query end (0-based BED coordinate).
        strand: Strand of the query (``'+'`` or ``'-'``).
        features_dict: Position-hashed dictionary of features.

    Returns:
        List of ``(feature, overlap_type)`` tuples where ``overlap_type``
        describes the spatial relationship (e.g. ``'query_contains_feature'``,
        ``'feature_contains_query'``, or ``'partial_by_N_bases'``).
    """
    features = []
    append = features.append

    start_key = int(qstart / HASH_VAL)
    end_key = int(qend / HASH_VAL)

    qstart = qstart + 1  # change 0-based bed to 1-based gff

    for i in range(start_key, end_key + 1):
        for feature in features_dict[chrom, i, strand]:
            if qstart <= feature.start and qend >= feature.end:
                append((feature, 'query_contains_feature'))
            elif qstart >= feature.start and qend <= feature.end:
                append((feature, 'feature_contains_query'))
            elif qstart <= feature.start and qend >= feature.start:
                append((
                    feature,
                    'partial_by_{}_bases'.format(qend - (feature.start - 1)),
                ))
            elif qstart <= feature.end and qend >= feature.end:
                append((
                    feature,
                    'partial_by_{}_bases'.format(feature.end - qstart),
                ))

    return features


def get_all_overlapping_features_from_query_unstranded(
        chrom, qstart, qend, strand, features_dict):
    """Return overlapping features without strict strand matching.

    When ``strand`` is None, returns features from both strands. When a strand
    is given, returns same-strand features first; if none are found, falls
    back to the opposite strand.

    Args:
        chrom: Query chromosome or sequence identifier.
        qstart: Query start (0-based BED coordinate).
        qend: Query end (0-based BED coordinate).
        strand: Query strand (``'+'``, ``'-'``, or None for both strands).
        features_dict: Position-hashed dictionary of features.

    Returns:
        List of ``(feature, overlap_type)`` tuples from the matching strand(s).
    """
    if strand is None:  # If no strand given, return features from both strands
        positive_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, '+', features_dict
        )
        negative_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, '-', features_dict
        )
        return positive_features + negative_features
    # return features on the same strand, only return features on
    # opposite strand if it can't find anything otherwise.
    else:

        priority_features = get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, strand, features_dict
        )
        if len(
                priority_features) == 0:  # found nothing on the proper strand
            if strand == '+':
                return get_all_overlapping_features_from_query_stranded(
                    chrom, qstart, qend, '-', features_dict
                )
            elif strand == '-':
                return get_all_overlapping_features_from_query_stranded(
                    chrom, qstart, qend, '+', features_dict
                )
            else:
                print("Strand not correct: {}".format(strand))
                return []
        else:
            return priority_features


def get_all_overlapping_features_from_query(
        chrom, qstart, qend, strand, features_dict, stranded=True):
    """Return overlapping features from a query interval.

    Dispatches to stranded or unstranded overlap logic based on the
    ``stranded`` flag.

    Args:
        chrom: Query chromosome or sequence identifier.
        qstart: Query start (0-based BED coordinate).
        qend: Query end (0-based BED coordinate).
        strand: Query strand (``'+'``, ``'-'``, or None).
        features_dict: Position-hashed dictionary of features.
        stranded: If True, return only same-strand overlaps. If False, use
            unstranded logic (same-strand preferred, opposite-strand fallback).

    Returns:
        List of ``(feature, overlap_type)`` tuples.
    """
    if stranded:
        return get_all_overlapping_features_from_query_stranded(
            chrom, qstart, qend, strand, features_dict
        )
    else:
        return get_all_overlapping_features_from_query_unstranded(
            chrom, qstart, qend, strand, features_dict
            # TODO: implement better. For now, returns both strands
        )


def parse_annotation_string(features_string):
    """Split a pipe-delimited annotation string into individual features.

    Args:
        features_string: Pipe-delimited string of colon-separated feature
            annotations.

    Returns:
        List of individual feature annotation strings.
    """
    features = features_string.split('|')
    return features


def is_protein_coding(transcript_type):
    """Check whether a transcript type is protein-coding.

    Args:
        transcript_type: GTF transcript_type attribute value.

    Returns:
        True if ``transcript_type`` equals ``'protein_coding'``, else False.
    """
    if transcript_type == 'protein_coding':
        return True

    return False


def classify_transcript_type(
        transcript_type, gene_type, distinct_nc_tx=None
):
    """Classify a transcript into a broad functional category.

    Deprecates ``is_protein_coding()`` with finer-grained non-coding
    classification.

    Args:
        transcript_type: GTF transcript_type attribute value.
        gene_type: GTF gene_type attribute value (currently unused).
        distinct_nc_tx: Optional list of non-coding transcript types to
            preserve as distinct categories rather than binning into
            ``'non_coding'``. For example, ``['lincRNA', 'pseudogene']``.

    Returns:
        Classification string: ``'protein_coding'``, ``'miRNA'``, a value
        from ``distinct_nc_tx``, or ``'non_coding'``.
    """

    # if the transcript type is protein coding, easy. return 'protein coding'.
    if transcript_type == 'protein_coding':
        return 'protein_coding'
    # TODO: Is this protein coding?
    # elif gene_type == 'protein_coding' and transcript_type == 'retained_intron':
    #     return 'protein_coding'
    # if we specify any non_coding type to pull out, returns the type here.
    # if transcript_type == 'snoRNA':
    #     return 'snoRNA'
    # if transcript_type == 'scaRNA':
    #     return 'scaRNA'
    if transcript_type == 'miRNA':
        return 'miRNA'
    if distinct_nc_tx is not None:
        for distinct_type in distinct_nc_tx:
            if transcript_type == distinct_type:
                return distinct_type
    # if not either protein_coding or anything in distinct_nc_tx, return non_coding.
    return 'non_coding'


def return_highest_priority_feature(formatted_features, priority):
    """Return the highest-priority feature from a list of formatted features.

    Groups features by ``(transcript_type_class, feature_type)`` and selects
    the group with the highest priority. On ties, returns the first match.

    Args:
        formatted_features: List of colon-delimited feature annotation
            strings.
        priority: List of ``[transcript_type_class, feature_type]`` pairs
            in descending priority order.

    Returns:
        Tuple of ``(feature_type_key, feature_string_list)`` for the
        highest-priority group, where ``feature_type_key`` is a
        ``(transcript_type_class, feature_type)`` tuple.
    """
    if not formatted_features:
        raise ValueError("No formatted features to prioritize")

    # Build dict
    combined_dict = defaultdict(list)
    for feature_string in formatted_features:
        parts = feature_string.split(':')
        if len(parts) != 10:
            raise ValueError(
                "Malformed annotation string (expected 10 colon-separated "
                "fields, got {}): {!r}".format(len(parts), feature_string[:120])
            )
        transcript, start, end, strand, feature_type, gene_id, \
            gene_name, transcript_type_list, gene_type_list, overlap = parts
        transcript_type_list = transcript_type_list.split(',')
        for transcript_type, gene_type in zip(transcript_type_list, gene_type_list):
            combined_dict[
                classify_transcript_type(transcript_type, gene_type),
                feature_type
            ].append(
                feature_string
            )
    # print("COMBINED DICT STEP 1: {}".format(combined_dict))
    # return the highest one
    combined_dict = OrderedDict(
        combined_dict
    )
    if not combined_dict:
        raise ValueError("No features after grouping for prioritization")

    priority_copy = [list(item) for item in priority]
    for key, _ in combined_dict.items():
        if list(key) not in priority_copy:
            priority_copy.append(list(key))
    combined_dict = sorted(
        combined_dict.items(),
        key=lambda x: priority_copy.index([x[0][0], x[0][1]])
    )
    return combined_dict[0]


def prioritize_transcript(unique_transcript_features, transcript_priority):
    """Select one transcript per gene by transcript-level priority.

    Args:
        unique_transcript_features: Mapping of transcript IDs to lists of
            colon-delimited feature strings overlapping the query (e.g.
            ``{'ENSTXYZ': ['CDS:...', 'exon:...']}``).
        transcript_priority: List of ``(featuretype, transcripttype)`` tuples
            in descending priority order.

    Returns:
        ``defaultdict(list)`` mapping gene IDs to lists of highest-priority
        feature strings, one per gene. Returns ``'intergenic'`` if no
        transcripts overlap.
    """
    unique_genes = defaultdict(list)
    # print("unique transcript features: {}".format(unique_transcript_features))
    if len(unique_transcript_features.keys()) == 0:
        return 'intergenic'
    # print('all transcript features: {}'.format(unique_transcript_features.keys()))
    for transcript in unique_transcript_features.keys():
        top_transcript = return_highest_priority_feature(
            unique_transcript_features[transcript],
            transcript_priority
        )[1][0]  # [0] contains the dictionary key
        gene_list = top_transcript.split(':')[5].split(',')
        for gene in gene_list:
            unique_genes[gene].append(top_transcript)

    for g in unique_genes.keys():
        unique_genes[g] = return_highest_priority_feature(
            unique_genes[g],
            transcript_priority
        )[1]
    return unique_genes


def prioritize_genes(unique_genes, gene_priority):
    """Select the highest-priority gene from a set of candidate genes.

    Args:
        unique_genes: Mapping of gene IDs to lists of colon-delimited
            feature strings (output of ``prioritize_transcript``).
        gene_priority: List of ``(featuretype, transcripttype)`` tuples
            in descending priority order.

    Returns:
        Sorted list of colon-delimited feature strings for the
        highest-priority gene(s). Feature type labels are adjusted for
        ``miRNA`` and ``non_coding`` categories.
    """
    # print("UNIQUE GENES", unique_genes)
    final_transcripts = []
    final_transcripts_append = final_transcripts.append
    for gene, transcripts in unique_genes.items():
        for transcript in transcripts:
            final_transcripts_append(transcript)

    feature_type, final_genes = return_highest_priority_feature(
        final_transcripts, gene_priority
    )
    if feature_type[0] == 'miRNA':  # TODO: fix. hacky
        ff = []
        for f in sorted(final_genes):
            f = f.replace(
                'exon', 'miRNA'
            )
            ff.append(f)
        return ff
    elif feature_type[0] == 'non_coding':  # TODO: fix. hacky
        ff = []
        for f in sorted(final_genes):
            f = f.replace(
                'exon', 'noncoding_exon'
            ).replace(
                'intron', 'noncoding_intron'
            )
            ff.append(f)
        return ff

    else:
        return sorted(final_genes)


def prioritize_transcript_then_gene(
        formatted_features, transcript_priority, gene_priority, parent='gene'):
    """Prioritize features by transcript, then by gene.

    Groups features by transcript, selects the highest-priority feature per
    transcript, groups transcripts by gene, then selects the highest-priority
    gene.

    Args:
        formatted_features: List of colon-delimited feature annotation
            strings.
        transcript_priority: List of ``(featuretype, transcripttype)`` tuples
            for transcript-level prioritization.
        gene_priority: List of ``(featuretype, transcripttype)`` tuples for
            gene-level prioritization.
        parent: Feature type to exclude from transcript grouping (usually
            ``'gene'``).

    Returns:
        List of colon-delimited feature strings for the highest-priority
        gene(s) after transcript-then-gene sorting.
    """
    unique_transcripts = defaultdict(list)

    # Filter gene features (we deal with priority on a per-transcript basis)
    for feature_string in formatted_features:
        if feature_string.split(':')[4] != parent:
            transcript = feature_string.split(':')[0]
            unique_transcripts[transcript].append(
                feature_string
            )
    # Get one transcript per gene whose featuretype has highest priority.
    unique_genes = prioritize_transcript(
        unique_transcripts, transcript_priority
    )
    # Get one gene whose featuretype has highest priority
    top_gene = prioritize_genes(unique_genes, gene_priority)
    return top_gene


def split_string(priority_string, delimiter=':'):
    """Extract region, gene_id, and gene_name from a priority string.

    Args:
        priority_string: Colon-delimited annotation string with region,
            gene_id, and gene_name in fields 4-6 (0-indexed). The literal
            string ``'intergenic'`` is handled as a special case.
        delimiter: Field separator character.

    Returns:
        Tuple of ``(region, gene_id, gene_name)``. All three values are
        ``'intergenic'`` when the input is intergenic.
    """
    if priority_string == 'intergenic':
        return 'intergenic', 'intergenic', 'intergenic'
    else:
        fields = priority_string.split(delimiter)
        return fields[4], fields[5], fields[6]


def annotate(
        chrom, start, end, name, score, strand,
        stranded, transcript_priority, gene_priority,
        features_dict, cds_dict, keys):
    """Annotate a single genomic interval.

    Finds overlapping features, formats them as colon-delimited annotation
    strings, and applies transcript-then-gene priority to select the best
    annotation.

    Args:
        chrom: Chromosome or sequence identifier.
        start: Interval start (0-based BED coordinate).
        end: Interval end (0-based BED coordinate).
        name: BED name field.
        score: BED score field.
        strand: BED strand field (``'+'``, ``'-'``, or ``'.'``).
        stranded: If True, require same-strand overlaps.
        transcript_priority: List of ``(featuretype, transcripttype)`` tuples
            for transcript-level prioritization.
        gene_priority: List of ``(featuretype, transcripttype)`` tuples for
            gene-level prioritization.
        features_dict: Position-hashed dictionary of features.
        cds_dict: Mapping of transcript_id to CDS boundaries for UTR
            classification.
        keys: GTF attribute key mapping from ``af.get_keys()``.

    Returns:
        Tuple of ``(chrom, start, end, name, score, strand, gene_id,
        gene_name, region, annotation_string)``. Intergenic intervals
        return ``'intergenic'`` for gene_id, gene_name, region, and
        annotation_string.
    """
    append_count = 0  # number of features that were found to overlap region

    overlapping_features = get_all_overlapping_features_from_query(
        chrom, start, end, strand, features_dict, stranded
    )
    annotation_parts = []
    transcript = defaultdict(list)
    for feature, overlap in overlapping_features:
        if keys['transcript_id'] in feature.attributes.keys():
            for transcript_id in feature.attributes[keys['transcript_id']]:
                transcript[transcript_id].append((feature, overlap))

    if len(transcript.keys()) == 0:
        return chrom, start, end, name, score, strand, \
            'intergenic', 'intergenic', \
            'intergenic', 'intergenic'

    append_count = 0
    for transcript_id, features in transcript.items():
        for feature, overlap in features:
            append_count += 1
            feature_type = feature.featuretype
            if feature_type is None:
                print("Warning: featuretype for feature {} is None".format(feature))
            elif feature_type == keys['utr']:
                feature_type = af.classify_utr(
                    feature, cds_dict, keys['transcript_id']
                )
            part = "{}:{}:{}:{}:{}:".format(
                transcript_id,
                feature.start - 1,
                feature.end,
                feature.strand,
                feature_type,
            )
            if keys['gene_id'] in feature.attributes.keys():
                part += ','.join(feature.attributes[keys['gene_id']]) + ':'
            else:
                part += '-:'
            if keys['gene_name'] in feature.attributes.keys():
                part += ','.join(feature.attributes[keys['gene_name']]) + ':'
            elif keys['transcript_id'] in feature.attributes.keys():
                part += ','.join(feature.attributes[keys['transcript_id']]) + ':'
            else:
                part += '-:'
            if keys['transcript_type'] in feature.attributes.keys():
                part += ','.join(feature.attributes[keys['transcript_type']]) + ':'
            else:
                part += '-:'
            if keys['gene_type'] in feature.attributes.keys():
                part += ','.join(feature.attributes[keys['gene_type']]) + ':'
            else:
                part += '-:'
            part += overlap + '|'
            annotation_parts.append(part)

    to_append = ''.join(annotation_parts)[:-1]
    if append_count == 1:  # if just one region, just return it.
        region, gene, rname = split_string(to_append)
    else:
        priority = prioritize_transcript_then_gene(
            parse_annotation_string(to_append),
            transcript_priority,
            gene_priority
        )
        region, gene, rname = split_string(priority[0])
        for highest_priority_gene in priority[1:]:
            current_region, current_gene, current_rname = split_string(
                highest_priority_gene
            )
            # there should really just be one highest priority region.
            if region != current_region:
                raise ValueError(
                    "Conflicting priority regions: {!r} vs {!r}".format(
                        region, current_region
                    )
                )
            # if the priority region contains multiple genes, return all
            if gene != current_gene:
                gene = gene + ',' + current_gene
            if rname != current_rname:
                rname = rname + ',' + current_rname
    return chrom, start, end, name, score, strand, \
        gene, rname, region, to_append


def split_bed_interval(bed_string):
    """Parse a single BED6 line into typed fields.

    Args:
        bed_string: Tab-delimited BED line string.

    Returns:
        Tuple of ``(chrom, start, end, name, score, strand)`` where
        ``start`` and ``end`` are integers and the remaining fields are
        strings.

    Raises:
        ValueError: If the line does not contain six tab-delimited fields
            or start/end are not integers.
    """
    fields = bed_string.rstrip('\n').split('\t')
    if len(fields) < 6:
        raise ValueError(
            "BED line must contain at least 6 tab-delimited fields, got {}: {!r}".format(
                len(fields), bed_string.rstrip()[:120]
            )
        )
    try:
        start = int(fields[1])
        end = int(fields[2])
    except ValueError as exc:
        raise ValueError(
            "BED start/end must be integers in line: {!r}".format(
                bed_string.rstrip()[:120]
            )
        ) from exc
    return fields[0], start, end, fields[3], fields[4], fields[5]


def regions_from_bed_lines(lines, hash_val=None):
    """Collect chromosomes and spatial hash buckets touched by BED lines.

    Skips blank lines and comment lines starting with ``#``.

    Args:
        lines: Iterable of BED6 line strings.
        hash_val: Bucket size in bases. Defaults to ``HASH_VAL``.

    Returns:
        Tuple of ``(chroms, buckets)`` where ``chroms`` is a set of contig
        names and ``buckets`` is a set of ``(chrom, bucket)`` pairs.
    """
    if hash_val is None:
        hash_val = HASH_VAL
    chroms = set()
    buckets = set()
    for line in lines:
        if not line or line.startswith('#'):
            continue
        stripped = line.strip()
        if not stripped:
            continue
        chrom, start, end, _name, _score, _strand = split_bed_interval(line)
        chroms.add(chrom)
        start_bucket = int(start / hash_val)
        end_bucket = int(end / hash_val)
        for bucket in range(start_bucket, end_bucket + 1):
            buckets.add((chrom, bucket))
    return chroms, buckets


def read_bed_files(bed_files):
    """Read BED files and aggregate chromosomes / buckets across them.

    Args:
        bed_files: Sequence of BED6 file paths.

    Returns:
        Tuple of ``(lines_per_file, chroms, buckets)`` where
        ``lines_per_file`` is a list of line lists parallel to ``bed_files``,
        ``chroms`` is the union of contig names, and ``buckets`` is the union
        of ``(chrom, bucket)`` keys.
    """
    lines_per_file = []
    chroms = set()
    buckets = set()
    for bed_file in bed_files:
        with open(bed_file, 'r') as bed_handle:
            lines = bed_handle.readlines()
        lines_per_file.append(lines)
        file_chroms, file_buckets = regions_from_bed_lines(lines)
        chroms |= file_chroms
        buckets |= file_buckets
    return lines_per_file, chroms, buckets


def annotate_bed_single_core(line, stranded, transcript_priority,
                             gene_priority, cds_dict, features_dict,
                             keys, progress=None):
    """Annotate a single BED line on one core.

    Args:
        line: Tab-delimited BED6 line string.
        stranded: If True, require same-strand overlaps.
        transcript_priority: List of ``(featuretype, transcripttype)`` tuples
            for transcript-level prioritization.
        gene_priority: List of ``(featuretype, transcripttype)`` tuples for
            gene-level prioritization.
        cds_dict: Mapping of transcript_id to CDS boundaries.
        features_dict: Position-hashed dictionary of features.
        keys: GTF attribute key mapping.
        progress: Optional ``tqdm`` progress bar to increment after
            annotation.

    Returns:
        Annotation result tuple from ``annotate()``.
    """
    chrom, start, end, name, score, strand = split_bed_interval(line)
    result = annotate(
        chrom, start, end, name, score, strand,
        stranded, transcript_priority, gene_priority,
        features_dict, cds_dict, keys
    )
    if progress is not None:
        progress.update(1)
    return result


_MP_STATE = {}


def _init_worker(state):
    """Initialize multiprocessing worker with shared annotation state.

    Args:
        state: Dictionary containing ``stranded``, ``transcript_priority``,
            ``gene_priority``, ``features_dict``, ``cds_dict``, and ``keys``.
    """
    global _MP_STATE
    _MP_STATE = state


def _annotate_chunk(chunk):
    """Annotate a chunk of BED lines using worker-global state.

    Args:
        chunk: List of tab-delimited BED6 line strings.

    Returns:
        List of annotation result tuples from ``annotate()``.
    """
    rows = []
    for line in chunk:
        chrom, start, end, name, score, strand = split_bed_interval(line)
        rows.append(annotate(
            chrom, start, end, name, score, strand,
            _MP_STATE['stranded'],
            _MP_STATE['transcript_priority'],
            _MP_STATE['gene_priority'],
            _MP_STATE['features_dict'],
            _MP_STATE['cds_dict'],
            _MP_STATE['keys'],
        ))
    return rows


def _chunk_lines(lines, cores):
    """Split BED lines into roughly equal chunks for parallel processing.

    Args:
        lines: List of BED line strings.
        cores: Number of worker processes.

    Returns:
        List of line chunks. Returns a single-element list containing all
        lines when ``cores <= 1`` or only one line is present.
    """
    if cores <= 1 or len(lines) <= 1:
        return [lines]
    chunk_size = max(1, (len(lines) + cores - 1) // cores)
    return [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]


def annotate_bed_parallel(lines, stranded, transcript_priority, gene_priority,
                          cds_dict, features_dict, keys, cores):
    """Annotate BED lines using single-core or multiprocessing execution.

    Uses fork-based multiprocessing on POSIX systems. Falls back to
    single-core execution when ``cores <= 1`` or on non-POSIX platforms.

    Args:
        lines: List of tab-delimited BED6 line strings.
        stranded: If True, require same-strand overlaps.
        transcript_priority: List of ``(featuretype, transcripttype)`` tuples
            for transcript-level prioritization.
        gene_priority: List of ``(featuretype, transcripttype)`` tuples for
            gene-level prioritization.
        cds_dict: Mapping of transcript_id to CDS boundaries.
        features_dict: Position-hashed dictionary of features.
        keys: GTF attribute key mapping.
        cores: Number of worker processes.

    Returns:
        List of annotation result tuples in the same order as ``lines``.
    """
    chunks = _chunk_lines(lines, cores)
    if cores <= 1 or len(chunks) <= 1:
        progress = trange(len(lines), desc='annotating', leave=False)
        results = [
            annotate_bed_single_core(
                line, stranded, transcript_priority, gene_priority,
                cds_dict, features_dict, keys, progress
            )
            for line in lines
        ]
        progress.close()
        return results

    if os.name != 'posix':
        print(
            "Warning: multiprocessing requires POSIX; running serially.",
            file=sys.stderr,
        )
        return annotate_bed_parallel(
            lines, stranded, transcript_priority, gene_priority,
            cds_dict, features_dict, keys, cores=1
        )

    state = {
        'stranded': stranded,
        'transcript_priority': transcript_priority,
        'gene_priority': gene_priority,
        'features_dict': features_dict,
        'cds_dict': cds_dict,
        'keys': keys,
    }
    ctx = multiprocessing.get_context('fork')
    with ctx.Pool(cores, initializer=_init_worker, initargs=(state,)) as pool:
        chunk_results = pool.map(_annotate_chunk, chunks)

    results = []
    for chunk in chunk_results:
        results.extend(chunk)
    return results


def annotate_bed(db_file, bed_files, out_files, stranded, chroms,
                 transcript_priority, gene_priority, gtf_format, append_chr, fuzzy,
                 cores, use_index=True, auto_chroms=True, triage_buckets=True,
                 lazy_buckets=False, max_cached_buckets=256):
    """Annotate BED6 files and write tab-delimited output.

    For each input BED file, appends annotation columns (gene_id, gene_name,
    region, and full pipe-delimited annotation string) and writes results to
    the corresponding output file.

    When ``chroms`` is empty and ``auto_chroms`` is True, chromosomes are
    inferred from the input BED intervals. When a companion index is used and
    ``triage_buckets`` is True, only 1 Mb spatial buckets overlapping the BED
    are loaded (or loaded lazily if ``lazy_buckets`` is True).

    Args:
        db_file: Path to the gffutils SQLite database file.
        bed_files: List of unannotated BED6 input file paths.
        out_files: List of output file paths (parallel to ``bed_files``).
        stranded: If True, require same-strand overlaps.
        chroms: List of chromosome names to restrict annotation indexing.
            Empty means auto-detect from BED when ``auto_chroms`` is True,
            otherwise load all chromosomes.
        transcript_priority: List of ``(featuretype, transcripttype)`` tuples
            for transcript-level prioritization.
        gene_priority: List of ``(featuretype, transcripttype)`` tuples for
            gene-level prioritization.
        gtf_format: Annotation syntax specifier. One of ``'gencode'``,
            ``'ensembl'``, or ``'refseq'``.
        append_chr: If True, prepend ``'chr'`` to database chromosome names
            (needed for assemblies using ``I``/``II``/``1`` instead of
            ``chrI``/``chr1``).
        fuzzy: Fuzzy-matching parameter passed to feature indexing.
        cores: Number of worker processes for parallel annotation.
        use_index: If True, load a pre-built annotator index when available.
        auto_chroms: If True and ``chroms`` is empty, restrict loading to
            contigs present in the input BED files.
        triage_buckets: If True and an index is used, load only buckets that
            intersect input BED intervals.
        lazy_buckets: If True, use an LRU disk-backed bucket cache instead of
            eagerly loading triaged buckets. Forced off when ``cores > 1``.
        max_cached_buckets: LRU size when ``lazy_buckets`` is True.

    Returns:
        Integer ``0`` on success.
    """

    lines_per_file, bed_chroms, bed_buckets = read_bed_files(bed_files)

    # Chrom limits for create_definitions / hash_features must be DB seqids.
    # BED auto-detection and index bucket keys use BED-style names.
    effective_chroms = list(chroms) if chroms else []
    if not effective_chroms and auto_chroms and bed_chroms:
        effective_chroms = sorted(
            db_chrom_from_bed(c, append_chr) for c in bed_chroms
        )
    if effective_chroms and bed_buckets:
        allowed_bed_chroms = {
            bed_chrom_from_db(c, append_chr) for c in effective_chroms
        }
        allowed_bed_chroms.update(effective_chroms)
        bed_buckets = {
            (chrom, bucket) for chrom, bucket in bed_buckets
            if chrom in allowed_bed_chroms
        }

    use_lazy = bool(lazy_buckets) and cores <= 1
    if lazy_buckets and cores > 1:
        print(
            "Warning: lazy bucket loading requires a single core; "
            "loading triaged buckets eagerly.",
            file=sys.stderr,
        )

    buckets = None
    if triage_buckets and use_index and not use_lazy:
        buckets = bed_buckets

    exons_dict, transcripts_dict, \
        cds_dict, features_dict, keys = create_definitions(
            db_file, chroms=effective_chroms, gtf_format=gtf_format,
            append_chr=append_chr, fuzzy=fuzzy, use_index=use_index,
            buckets=buckets, lazy_buckets=use_lazy,
            max_cached_buckets=max_cached_buckets,
        )
    try:
        for n in range(len(bed_files)):
            out_file = out_files[n]
            results = annotate_bed_parallel(
                lines_per_file[n], stranded, transcript_priority, gene_priority,
                cds_dict, features_dict, keys, cores
            )
            with open(out_file, 'w') as out_handle:
                writer = csv.writer(out_handle, delimiter='\t')
                writer.writerows(results)
    finally:
        close = getattr(features_dict, 'close', None)
        if callable(close):
            close()
    return 0
