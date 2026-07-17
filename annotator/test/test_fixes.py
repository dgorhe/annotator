import os
import time
from pathlib import Path

import pytest

from annotator import annotate_bed as ab
from annotator.index_store import (
    build_annotator_index,
    index_is_current,
    load_annotator_index,
    _index_covers_chroms,
)
from annotator.miRNA_name2id import build_precursor_id_to_mature_id_dict


@pytest.fixture(scope="module")
def chr21_db():
    from annotator.test.db_fixtures import ensure_db
    return ensure_db('gencode.v19.annotation.chr21.46M-47M.gtf')


def test_index_is_current_false_after_db_changes(chr21_db, tmp_path):
    import shutil

    local_db = tmp_path / "chr21.db"
    shutil.copy2(chr21_db, local_db)
    index_file = str(tmp_path / "test.annotator.idx")
    build_annotator_index(
        str(local_db),
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )
    assert index_is_current(
        index_file, str(local_db), chroms=['chr21'], gtf_format='gencode'
    )

    local_db.touch()
    time.sleep(0.01)
    assert not index_is_current(
        index_file, str(local_db), chroms=['chr21'], gtf_format='gencode'
    )


def test_stale_index_falls_back_to_legacy(chr21_db, tmp_path, capsys):
    import shutil

    local_db = tmp_path / "chr21.db"
    shutil.copy2(chr21_db, local_db)
    sidecar = str(local_db) + ".annotator.idx"
    build_annotator_index(
        str(local_db),
        index_file=sidecar,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )
    local_db.touch()
    time.sleep(0.01)

    bed = tmp_path / "input.bed"
    bed.write_text("chr21\t46553734\t46553788\tADARB1\t0\t+\n")
    out_stale = tmp_path / "stale.bed"
    out_legacy = tmp_path / "legacy.bed"

    ab.annotate_bed(
        str(local_db), [str(bed)], [str(out_stale)],
        stranded=True, chroms=['chr21'],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=True,
    )
    ab.annotate_bed(
        str(local_db), [str(bed)], [str(out_legacy)],
        stranded=True, chroms=['chr21'],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=False,
    )

    assert "stale or incompatible" in capsys.readouterr().err
    assert out_stale.read_text() == out_legacy.read_text()


def test_load_rejects_hash_val_mismatch(chr21_db, tmp_path_factory):
    index_file = str(tmp_path_factory.mktemp("idx") / "bad.annotator.idx")
    build_annotator_index(
        chr21_db,
        index_file=index_file,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )

    import sqlite3
    conn = sqlite3.connect(index_file)
    conn.execute(
        "UPDATE meta SET value = ? WHERE key = 'hash_val'",
        ("9999999",),
    )
    conn.commit()
    conn.close()

    with pytest.raises(ValueError, match="hash_val"):
        load_annotator_index(index_file, chroms=['chr21'], gtf_format='gencode')


def test_mirna_custom_file_extends_mature_ids(tmp_path):
    custom = tmp_path / "custom.tsv"
    custom.write_text(
        "prec1\t_\tm1|m2\n"
        "prec1\t_\tm3\n"
    )
    mapping = build_precursor_id_to_mature_id_dict(None, str(custom))
    assert mapping["prec1"] == ["m1", "m2", "m3"]


def test_parse_annotation_priority_rejects_short_lines(tmp_path):
    from annotator.annotator import parse_annotation_priority

    bad = tmp_path / "bad_priority.txt"
    bad.write_text("only_one_field\n")
    with pytest.raises(ValueError, match="at least two fields"):
        parse_annotation_priority(str(bad))


def test_partial_index_rejected_for_full_genome(chr21_db, tmp_path, capsys):
    import shutil

    local_db = tmp_path / "chr21.db"
    shutil.copy2(chr21_db, local_db)
    sidecar = str(local_db) + ".annotator.idx"
    build_annotator_index(
        str(local_db),
        index_file=sidecar,
        chroms=['chr21'],
        gtf_format='gencode',
        force=True,
    )

    meta = {'all_chromosomes': '0', 'indexed_chroms': 'chr21'}
    assert not _index_covers_chroms(meta, [])
    assert _index_covers_chroms(meta, ['chr21'])
    assert not index_is_current(
        sidecar, str(local_db), chroms=[], gtf_format='gencode'
    )

    bed = tmp_path / "input.bed"
    bed.write_text("chr21\t46553734\t46553788\tADARB1\t0\t+\n")
    out_partial = tmp_path / "partial.bed"
    out_legacy = tmp_path / "legacy.bed"

    ab.annotate_bed(
        str(local_db), [str(bed)], [str(out_partial)],
        stranded=True, chroms=[],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=True,
    )
    ab.annotate_bed(
        str(local_db), [str(bed)], [str(out_legacy)],
        stranded=True, chroms=[],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
        use_index=False,
    )

    assert "stale or incompatible" in capsys.readouterr().err
    assert out_partial.read_text() == out_legacy.read_text()


def test_split_bed_interval_rejects_short_lines():
    with pytest.raises(ValueError, match="at least 6"):
        ab.split_bed_interval("chr1\t100\t200\n")


def test_malformed_annotation_string_raises():
    with pytest.raises(ValueError, match="Malformed annotation string"):
        ab.return_highest_priority_feature(["too:few:fields"], [])


def test_gene_id_to_name_skips_missing_keys_by_default():
    import gffutils
    from annotator import annotation_functions as af

    class FakeGene:
        def __init__(self, attrs):
            self.attributes = attrs

    class FakeDB:
        def features_of_type(self, _featuretype):
            return [
                FakeGene({'gene_id': ['G1'], 'gene_name': ['NAME1']}),
                FakeGene({'gene_id': ['G2']}),
            ]

    mapping = af.gene_id_to_name(FakeDB(), 'gene_name')
    assert mapping == {'G1': 'NAME1'}


def test_core_import_without_optional_deps():
    import annotator.annotate_bed
    import annotator.annotator
    import annotator.index_store


def test_region_bedfiles_import_requires_extra():
    with pytest.raises(ImportError, match="pybedtools"):
        from annotator import create_region_bedfiles
        create_region_bedfiles._require_pybedtools()


@pytest.mark.skipif(
    os.environ.get("ANNOTATOR_BENCHMARK_BED") is None
    and not Path(__file__).resolve().parents[2].joinpath(
        "benchmark", "Ctr_LARP1.bed"
    ).exists(),
    reason="benchmark BED not available",
)
def test_benchmark_bed_serial_runs(tmp_path):
    from annotator import annotator as ann

    repo_root = Path(__file__).resolve().parents[2]
    input_bed = os.environ.get(
        "ANNOTATOR_BENCHMARK_BED",
        str(repo_root / "benchmark" / "Ctr_LARP1.bed"),
    )
    gtf_db = repo_root / "benchmark/gencode.v50.primary_assembly.annotation.gtf.db"
    if not gtf_db.exists():
        pytest.skip("benchmark GTF database not built yet")

    out = tmp_path / "Ctr_LARP1.annotated.bed"
    ab.annotate_bed(
        str(gtf_db), [input_bed], [str(out)],
        stranded=True, chroms=[],
        transcript_priority=[list(x) for x in ann.TRANSCRIPT_PRIORITY],
        gene_priority=[list(x) for x in ann.GENE_PRIORITY],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
    )
    lines = out.read_text().strip().splitlines()
    assert len(lines) == 648
