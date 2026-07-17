import time

import pytest

from annotator import annotate_bed as ab
from annotator.index_store import build_annotator_index, load_annotator_index


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
