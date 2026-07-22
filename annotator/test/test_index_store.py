import time

import pytest

from annotator import annotate_bed as ab
from annotator.index_store import (
    build_annotator_index,
    load_annotator_index,
    query_bucket_features,
    _connect,
)


@pytest.fixture(scope="module")
def chr21_db():
    from annotator.test.db_fixtures import ensure_db
    return ensure_db('gencode.v19.annotation.chr21.46M-47M.gtf')


def _feature_set(features_dict, transcript_id_key='transcript_id'):
    result = {}
    for key, features in features_dict.items():
        result[key] = {
            (f.start, f.end, f.featuretype,
             tuple(f.attributes.get(transcript_id_key, [])))
            for f in features
        }
    return result


def test_index_load_matches_legacy_definitions(chr21_db, tmp_path_factory):
    index_file = str(tmp_path_factory.mktemp("idx") / "test.annotator.idx")

    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )

    _, _, cds_legacy, features_legacy, keys = ab.create_definitions(
        chr21_db, ['chr21'], 'gencode', use_index=False
    )
    cds_index, features_index, keys = load_annotator_index(
        index_file, chroms=['chr21'], gtf_format='gencode'
    )

    assert cds_legacy == cds_index
    assert _feature_set(features_legacy, keys['transcript_id']) == _feature_set(
        features_index, keys['transcript_id']
    )


def test_index_load_is_faster_than_legacy_build(chr21_db, tmp_path_factory):
    index_file = str(tmp_path_factory.mktemp("idx") / "test.annotator.idx")
    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )

    t0 = time.time()
    ab.create_definitions(chr21_db, ['chr21'], 'gencode', use_index=False)
    legacy_time = time.time() - t0

    t0 = time.time()
    load_annotator_index(index_file, chroms=['chr21'], gtf_format='gencode')
    index_time = time.time() - t0

    assert index_time < legacy_time


def test_annotated_output_uses_index(chr21_db, tmp_path):
    index_file = chr21_db + ".annotator.idx"
    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )

    bed = tmp_path / "input.bed"
    bed.write_text("chr21\t46553734\t46553788\tADARB1\t0\t+\n")
    out_legacy = tmp_path / "legacy.bed"
    out_index = tmp_path / "index.bed"

    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_legacy)],
        stranded=True, chroms=['chr21'],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=False,
    )
    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_index)],
        stranded=True, chroms=['chr21'],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=True,
    )
    assert out_legacy.read_text() == out_index.read_text()


def test_regions_from_bed_lines_collects_chroms_and_buckets():
    lines = [
        "# comment\n",
        "chr21\t46553734\t46553788\tADARB1\t0\t+\n",
        "chr1\t999999\t1000001\tx\t0\t-\n",
        "\n",
    ]
    chroms, buckets = ab.regions_from_bed_lines(lines, hash_val=1000000)
    assert chroms == {'chr21', 'chr1'}
    assert ('chr21', 46) in buckets
    assert ('chr1', 0) in buckets
    assert ('chr1', 1) in buckets


def test_append_chr_chrom_name_helpers():
    assert ab.db_chrom_from_bed('chrI', True) == 'I'
    assert ab.db_chrom_from_bed('I', True) == 'I'
    assert ab.db_chrom_from_bed('chr21', False) == 'chr21'
    assert ab.bed_chrom_from_db('I', True) == 'chrI'
    assert ab.bed_chrom_from_db('chrI', True) == 'chrI'
    assert ab.bed_chrom_from_db('chr21', False) == 'chr21'


def test_auto_chroms_append_chr_legacy_hashes_db_seqids(tmp_path):
    import shutil

    from annotator.test.db_fixtures import ensure_db

    src_db = ensure_db(
        'c_elegans.PRJNA13758.WS257.canonical_geneset.chrIII.25000.gtf'
    )
    db = str(tmp_path / "ce.db")
    shutil.copy2(src_db, db)
    bed = tmp_path / "input.bed"
    bed.write_text("chrIII\t2816\t3092\tcTel54X.1\t0\t-\n")
    out = tmp_path / "out.bed"
    priority = [['protein_coding', 'exon'], ['non_coding', 'exon']]

    ab.annotate_bed(
        db, [str(bed)], [str(out)],
        stranded=True, chroms=[],
        transcript_priority=priority, gene_priority=priority,
        gtf_format='ce11', append_chr=True, fuzzy=0, cores=1,
        use_index=False, auto_chroms=True,
    )
    text = out.read_text().strip()
    assert text
    assert 'intergenic' not in text.split('\t')[6:10]


def test_bucket_triage_with_append_chr_limit_chroms(tmp_path):
    import shutil

    from annotator.test.db_fixtures import ensure_db

    src_db = ensure_db(
        'c_elegans.PRJNA13758.WS257.canonical_geneset.chrIII.25000.gtf'
    )
    db = str(tmp_path / "ce.db")
    shutil.copy2(src_db, db)
    build_annotator_index(
        db,
        chroms=['III'],
        gtf_format='ce11',
        append_chr=True,
        force=True,
    )

    bed = tmp_path / "input.bed"
    bed.write_text("chrIII\t2816\t3092\tcTel54X.1\t0\t-\n")
    out = tmp_path / "out.bed"
    priority = [['protein_coding', 'exon'], ['non_coding', 'exon']]

    ab.annotate_bed(
        db, [str(bed)], [str(out)],
        stranded=True, chroms=['III'],
        transcript_priority=priority, gene_priority=priority,
        gtf_format='ce11', append_chr=True, fuzzy=0, cores=1,
        use_index=True, auto_chroms=False, triage_buckets=True,
    )
    text = out.read_text().strip()
    assert text
    assert 'intergenic' not in text.split('\t')[6:10]


def test_bucket_triage_loads_fewer_keys_than_full_chrom(chr21_db, tmp_path):
    index_file = str(tmp_path / "triage.annotator.idx")
    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )
    _, features_full, _ = load_annotator_index(
        index_file, chroms=['chr21'], gtf_format='gencode'
    )
    _, features_triaged, _ = load_annotator_index(
        index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        buckets={('chr21', 46)},
    )
    _, features_miss, _ = load_annotator_index(
        index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        buckets={('chr21', 0)},
    )
    assert set(features_triaged) <= set(features_full)
    assert all(bucket == 46 for _chrom, bucket, _strand in features_triaged)
    assert features_triaged
    assert not features_miss


def test_auto_chroms_and_bucket_triage_match_explicit(chr21_db, tmp_path):
    index_file = chr21_db + ".annotator.idx"
    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )
    bed = tmp_path / "input.bed"
    bed.write_text("chr21\t46553734\t46553788\tADARB1\t0\t+\n")
    out_auto = tmp_path / "auto.bed"
    out_explicit = tmp_path / "explicit.bed"
    priority = [['protein_coding', 'distintron500']]

    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_auto)],
        stranded=True, chroms=[],
        transcript_priority=priority, gene_priority=priority,
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=True, auto_chroms=True, triage_buckets=True,
    )
    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_explicit)],
        stranded=True, chroms=['chr21'],
        transcript_priority=priority, gene_priority=priority,
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=True, auto_chroms=False, triage_buckets=True,
    )
    assert out_auto.read_text() == out_explicit.read_text()


def test_lazy_buckets_match_eager_annotation(chr21_db, tmp_path):
    index_file = chr21_db + ".annotator.idx"
    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )
    bed = tmp_path / "input.bed"
    bed.write_text("chr21\t46553734\t46553788\tADARB1\t0\t+\n")
    out_eager = tmp_path / "eager.bed"
    out_lazy = tmp_path / "lazy.bed"
    priority = [['protein_coding', 'distintron500']]

    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_eager)],
        stranded=True, chroms=['chr21'],
        transcript_priority=priority, gene_priority=priority,
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=True, lazy_buckets=False,
    )
    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_lazy)],
        stranded=True, chroms=['chr21'],
        transcript_priority=priority, gene_priority=priority,
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=True, lazy_buckets=True,
    )
    assert out_eager.read_text() == out_lazy.read_text()


def test_query_bucket_features_returns_records(chr21_db, tmp_path):
    from annotator import annotation_functions as af

    index_file = str(tmp_path / "q.annotator.idx")
    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )
    keys = af.get_keys('gencode')
    conn = _connect(index_file)
    try:
        features = query_bucket_features(conn, 'chr21', 46, '+', keys)
    finally:
        conn.close()
    assert features
    assert all(f.chrom == 'chr21' for f in features)


def test_index_only_build_from_existing_db(chr21_db, tmp_path):
    from annotator.build_gffutils_db import build_index_from_db

    index_file = str(tmp_path / "from_db.annotator.idx")
    written = build_index_from_db(
        chr21_db,
        index_file=index_file,
        index_chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )
    assert written == index_file
    cds_dict, features_dict, keys = load_annotator_index(
        index_file, chroms=['chr21'], gtf_format='gencode'
    )
    assert cds_dict
    assert features_dict
    assert keys['transcript_id'] == 'transcript_id'
