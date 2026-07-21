"""Pre-built annotation index for fast BED lookup."""

import sqlite3
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from . import annotation_functions as af

INDEX_SUFFIX = ".annotator.idx"
INDEX_VERSION = "1"


@dataclass
class FeatureRecord:
    """Minimal feature representation used by the annotate hot path.

    Stores genomic coordinates and normalized attribute keys without the
    overhead of a full gffutils Feature object.

    Attributes:
        chrom (str): Chromosome or contig name.
        start (int): 1-based inclusive start coordinate (gffutils convention).
        end (int): 1-based inclusive end coordinate (gffutils convention).
        strand (str): Strand ('+' or '-').
        featuretype (str): GTF/GFF feature type (e.g. 'exon', 'CDS').
        attributes (dict): Mapping of normalized attribute keys to list
            values (e.g. transcript_id, gene_id).
    """

    chrom: str
    start: int
    end: int
    strand: str
    featuretype: str
    attributes: dict = field(default_factory=dict)


def index_path_for_db(db_file):
    """Return the default companion index path for a gffutils database.

    Args:
        db_file (str): Path to the gffutils SQLite database file.

    Returns:
        str: Path formed by appending INDEX_SUFFIX to db_file.
    """
    return db_file + INDEX_SUFFIX


def index_exists(db_file):
    """Check whether a companion annotator index exists for a database.

    Args:
        db_file (str): Path to the gffutils SQLite database file.

    Returns:
        bool: True if the index file at index_path_for_db(db_file) exists.
    """
    return Path(index_path_for_db(db_file)).exists()


def _read_index_meta(index_file):
    """Read metadata key/value pairs from an annotator index file.

    Args:
        index_file (str): Path to the ``.annotator.idx`` SQLite file.

    Returns:
        dict: Mapping of metadata keys to string values.
    """
    conn = _connect(index_file)
    meta = {
        row['key']: row['value']
        for row in conn.execute("SELECT key, value FROM meta")
    }
    conn.close()
    return meta


def index_is_current(index_file, db_file, chroms=None, gtf_format='gencode',
                     append_chr=False):
    """Return whether an index matches its source database and runtime settings.

    Args:
        index_file (str): Path to the ``.annotator.idx`` file.
        db_file (str): Path to the source gffutils SQLite database.
        chroms (list, optional): Chromosomes requested for this annotation run.
            An empty list means all chromosomes in the database.
        gtf_format (str): GTF nomenclature for the current run.
        append_chr (bool): Whether ``chr`` is prepended to contig names.

    Returns:
        bool: True if the index version, hash bucket size, source
            modification time, format settings, and chromosome coverage all
            match the requested run. Absolute ``source_db`` path is not
            required to match, so indexes copied alongside a database to a
            new location remain usable when the database mtime is preserved.
    """
    from . import annotate_bed as ab

    index_path = Path(index_file)
    db_path = Path(db_file)
    if not index_path.exists() or not db_path.exists():
        return False
    try:
        meta = _read_index_meta(index_file)
    except sqlite3.Error:
        return False
    if meta.get('version') != INDEX_VERSION:
        return False
    if meta.get('hash_val') != str(ab.HASH_VAL):
        return False
    # Do not require meta source_db path to match: sidecar indexes are often
    # copied to a new location with their database. Still require mtime match
    # so in-place database updates invalidate the index.
    expected_mtime = meta.get('source_db_mtime')
    if expected_mtime is not None:
        try:
            if str(int(db_path.stat().st_mtime)) != expected_mtime:
                return False
        except OSError:
            return False
    if meta.get('gtf_format') is not None and meta.get('gtf_format') != gtf_format:
        return False
    if meta.get('append_chr') is not None:
        if meta.get('append_chr') != ('1' if append_chr else '0'):
            return False
    if not _index_covers_chroms(meta, chroms or []):
        return False
    return True


def _index_covers_chroms(meta, runtime_chroms):
    """Return whether index metadata covers the chromosomes requested at runtime.

    Legacy indexes without chromosome metadata are treated as full-genome indexes
    for backwards compatibility.

    Args:
        meta (dict): Index metadata key/value pairs.
        runtime_chroms (list): Chromosomes for the current run; empty means all.

    Returns:
        bool: True if the index is sufficient for the requested chromosomes.
    """
    indexed_chroms = _indexed_chrom_set(meta)
    if not runtime_chroms:
        return indexed_chroms is None
    if indexed_chroms is None:
        return True
    return set(runtime_chroms) <= indexed_chroms


def _indexed_chrom_set(meta):
    """Parse indexed chromosome coverage from index metadata.

    Args:
        meta (dict): Index metadata key/value pairs.

    Returns:
        None if the index covers all chromosomes, else a set of contig names.
    """
    if meta.get('all_chromosomes') == '1':
        return None
    chroms_str = meta.get('indexed_chroms')
    if chroms_str is None:
        return None
    if chroms_str in ('', '*'):
        return None
    return set(chroms_str.split(','))


def _connect(path):
    """Open a SQLite connection configured for row-based access.

    Args:
        path (str): Filesystem path to the SQLite database.

    Returns:
        sqlite3.Connection: Connection with row_factory set to sqlite3.Row.
    """
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    return conn


def _create_schema(conn):
    """Create annotator index tables and indexes if they do not exist.

    Defines meta, cds, features, and bucket_features tables plus coordinate
    lookup indexes used during BED annotation.

    Args:
        conn (sqlite3.Connection): Open connection to the index database.
    """
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS meta (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS cds (
            transcript_id TEXT PRIMARY KEY,
            low INT NOT NULL,
            hi INT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS features (
            id INTEGER PRIMARY KEY,
            chrom TEXT NOT NULL,
            start INT NOT NULL,
            end INT NOT NULL,
            strand TEXT NOT NULL,
            featuretype TEXT NOT NULL,
            transcript_id TEXT NOT NULL,
            gene_id TEXT,
            gene_name TEXT,
            transcript_type TEXT,
            gene_type TEXT
        );
        CREATE INDEX IF NOT EXISTS idx_features_coords
            ON features(chrom, strand, start, end);
        CREATE TABLE IF NOT EXISTS bucket_features (
            chrom TEXT NOT NULL,
            bucket INT NOT NULL,
            strand TEXT NOT NULL,
            feature_id INT NOT NULL,
            PRIMARY KEY (chrom, bucket, strand, feature_id)
        ) WITHOUT ROWID;
        CREATE INDEX IF NOT EXISTS idx_bucket_lookup
            ON bucket_features(chrom, bucket, strand);
        """
    )


def _feature_to_record(feature, chrom, keys):
    """Convert a gffutils Feature into a FeatureRecord for index storage.

    Args:
        feature: gffutils Feature with start, end, strand, featuretype,
            and attributes.
        chrom (str): Chromosome name to assign to the record.
        keys (dict): Normalized GTF attribute key mapping from
            annotation_functions.get_keys.

    Returns:
        FeatureRecord or None: Populated record, or None if the feature
            lacks a transcript_id attribute.
    """
    attrs = feature.attributes
    transcript_ids = attrs.get(keys['transcript_id'], [])
    if not transcript_ids:
        return None
    return FeatureRecord(
        chrom=chrom,
        start=feature.start,
        end=feature.end,
        strand=feature.strand,
        featuretype=feature.featuretype,
        attributes={
            keys['transcript_id']: list(transcript_ids),
            keys['gene_id']: list(attrs.get(keys['gene_id'], [])),
            keys['gene_name']: list(attrs.get(keys['gene_name'], [])),
            keys['transcript_type']: list(attrs.get(keys['transcript_type'], [])),
            keys['gene_type']: list(attrs.get(keys['gene_type'], [])),
        },
    )


def build_annotator_index(
        db_file, index_file=None, chroms=None, gtf_format='gencode',
        append_chr=False, force=False):
    """Build a companion SQLite index from an existing gffutils database.

    Reads feature definitions via annotate_bed.create_definitions, writes
    CDS bounds and bucketed features to a .annotator.idx file for fast
    startup during BED annotation.

    Args:
        db_file (str): Path to the source gffutils SQLite database.
        index_file (str, optional): Output index path. Defaults to
            index_path_for_db(db_file).
        chroms (list, optional): Chromosomes to include. Empty list means
            all chromosomes present in the database.
        gtf_format (str): GTF nomenclature ('gencode' or ensembl/ncbi style).
        append_chr (bool): If True, prepend 'chr' to contig names when
            building definitions.
        force (bool): If True, overwrite an existing index file.

    Returns:
        str: Path to the written index file.

    Raises:
        FileExistsError: If index_file already exists and force is False.
    """
    from . import annotate_bed as ab

    index_file = index_file or index_path_for_db(db_file)
    path = Path(index_file)
    if path.exists():
        if force:
            path.unlink()
        else:
            raise FileExistsError(
                "Index already exists: {} (use force=True to rebuild)".format(
                    index_file
                )
            )

    _, _, cds_dict, features_dict, keys = ab.create_definitions(
        db_file,
        chroms=chroms or [],
        gtf_format=gtf_format,
        append_chr=append_chr,
        use_index=False,
    )

    db_path = Path(db_file)
    conn = _connect(index_file)
    _create_schema(conn)
    meta = conn.cursor()
    chrom_list = chroms or []
    meta.executemany(
        "INSERT INTO meta(key, value) VALUES (?, ?)",
        [
            ("version", INDEX_VERSION),
            ("hash_val", str(ab.HASH_VAL)),
            ("gtf_format", gtf_format),
            ("source_db", str(db_path.resolve())),
            ("source_db_mtime", str(int(db_path.stat().st_mtime))),
            ("append_chr", "1" if append_chr else "0"),
            ("all_chromosomes", "1" if not chrom_list else "0"),
            ("indexed_chroms", "*" if not chrom_list else ",".join(sorted(chrom_list))),
        ],
    )
    meta.executemany(
        "INSERT INTO cds(transcript_id, low, hi) VALUES (?, ?, ?)",
        [
            (tid, bounds['low'], bounds['hi'])
            for tid, bounds in cds_dict.items()
            if 'low' in bounds and 'hi' in bounds
        ],
    )

    feature_key_to_id = {}
    feature_rows = []
    bucket_rows = []

    def feature_identity(record):
        """Return a hashable tuple identifying a feature record.

        Args:
            record (FeatureRecord): Feature to identify.

        Returns:
            tuple: Tuple of chrom, start, end, strand, featuretype,
                transcript_id, and gene_id values.
        """
        return (
            record.chrom,
            record.start,
            record.end,
            record.strand,
            record.featuretype,
            record.attributes.get(keys['transcript_id'], [''])[0],
            tuple(record.attributes.get(keys['gene_id'], [])),
        )

    for bucket_key, feature_list in features_dict.items():
        chrom, bucket, strand = bucket_key
        for feature in feature_list:
            record = _feature_to_record(feature, chrom, keys)
            if record is None:
                continue
            identity = feature_identity(record)
            feature_id = feature_key_to_id.get(identity)
            if feature_id is None:
                feature_id = len(feature_rows) + 1
                feature_key_to_id[identity] = feature_id
                attrs = record.attributes
                feature_rows.append((
                    feature_id,
                    record.chrom,
                    record.start,
                    record.end,
                    record.strand,
                    record.featuretype,
                    attrs.get(keys['transcript_id'], [''])[0],
                    ','.join(attrs.get(keys['gene_id'], [])) or None,
                    ','.join(attrs.get(keys['gene_name'], [])) or None,
                    ','.join(attrs.get(keys['transcript_type'], [])) or None,
                    ','.join(attrs.get(keys['gene_type'], [])) or None,
                ))
            bucket_rows.append((chrom, bucket, strand, feature_id))

    conn.executemany(
        """
        INSERT INTO features(
            id, chrom, start, end, strand, featuretype,
            transcript_id, gene_id, gene_name, transcript_type, gene_type
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        feature_rows,
    )
    conn.executemany(
        """
        INSERT OR IGNORE INTO bucket_features(chrom, bucket, strand, feature_id)
        VALUES (?, ?, ?, ?)
        """,
        bucket_rows,
    )
    conn.commit()
    conn.close()
    return index_file


def load_annotator_index(index_file, chroms=None, gtf_format=None):
    """Load cds_dict and features_dict from a pre-built annotator index.

    Reconstructs in-memory structures compatible with annotate_bed without
    re-querying the full gffutils database.

    Args:
        index_file (str): Path to the .annotator.idx SQLite file.
        chroms (list, optional): If provided, only features on these
            chromosomes are loaded.
        gtf_format (str, optional): GTF nomenclature for attribute keys.
            Defaults to the value stored in index metadata.

    Returns:
        tuple: (cds_dict, features_dict, keys) where cds_dict maps
            transcript_id to {'low', 'hi'} CDS bounds, features_dict maps
            (chrom, bucket, strand) to lists of FeatureRecord objects, and
            keys is the attribute key mapping from get_keys.

    Raises:
        ValueError: If the index version or hash bucket size is unsupported.
    """
    from . import annotate_bed as ab

    conn = _connect(index_file)
    meta = {
        row['key']: row['value']
        for row in conn.execute("SELECT key, value FROM meta")
    }
    if meta.get('version') != INDEX_VERSION:
        conn.close()
        raise ValueError("Unsupported index version: {}".format(meta.get('version')))
    if meta.get('hash_val') != str(ab.HASH_VAL):
        conn.close()
        raise ValueError(
            "Index hash_val {} does not match runtime HASH_VAL {}".format(
                meta.get('hash_val'), ab.HASH_VAL
            )
        )

    gtf_format = gtf_format or meta.get('gtf_format', 'gencode')
    keys = af.get_keys(gtf_format)

    cds_dict = defaultdict(dict)
    for row in conn.execute("SELECT transcript_id, low, hi FROM cds"):
        cds_dict[row['transcript_id']] = {'low': row['low'], 'hi': row['hi']}

    features_by_id = {}
    for row in conn.execute("SELECT * FROM features"):
        attrs = {keys['transcript_id']: [row['transcript_id']]}
        if row['gene_id']:
            attrs[keys['gene_id']] = row['gene_id'].split(',')
        if row['gene_name']:
            attrs[keys['gene_name']] = row['gene_name'].split(',')
        if row['transcript_type']:
            attrs[keys['transcript_type']] = row['transcript_type'].split(',')
        if row['gene_type']:
            attrs[keys['gene_type']] = row['gene_type'].split(',')
        features_by_id[row['id']] = FeatureRecord(
            chrom=row['chrom'],
            start=row['start'],
            end=row['end'],
            strand=row['strand'],
            featuretype=row['featuretype'],
            attributes=attrs,
        )

    chrom_filter = set(chroms) if chroms else None
    bucket_sql = (
        "SELECT chrom, bucket, strand, feature_id FROM bucket_features"
    )
    params = []
    if chrom_filter:
        placeholders = ",".join("?" for _ in chrom_filter)
        bucket_sql += " WHERE chrom IN ({})".format(placeholders)
        params = list(chrom_filter)

    features_dict = defaultdict(list)
    for row in conn.execute(bucket_sql, params):
        feature = features_by_id.get(row['feature_id'])
        if feature is not None:
            features_dict[(row['chrom'], row['bucket'], row['strand'])].append(
                feature
            )

    conn.close()
    return cds_dict, features_dict, keys
