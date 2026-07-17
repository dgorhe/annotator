import pytest

from annotator import annotate_bed as ab
from annotator.intervals import Interval


def test_split_bed_interval_parses_bed6():
    chrom, start, end, name, score, strand = ab.split_bed_interval(
        "chr1\t100\t200\tpeak\t0\t+\n"
    )
    assert chrom == "chr1"
    assert start == 100
    assert end == 200
    assert name == "peak"
    assert score == "0"
    assert strand == "+"


def test_split_bed_interval_rejects_non_integer_coords():
    with pytest.raises(ValueError, match="must be integers"):
        ab.split_bed_interval("chr1\tx\t200\tpeak\t0\t+\n")


def test_cores_flag_rejects_invalid_value(capsys):
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "-m", "annotator.annotator", "-c", "0",
         "--input", "x.bed", "--gtfdb", "x.db"],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0


def test_parallel_matches_serial(tmp_path, chr21_db):
    bed = tmp_path / "input.bed"
    bed.write_text("chr21\t46553734\t46553788\tADARB1\t0\t+\n")
    out_serial = tmp_path / "serial.bed"
    out_parallel = tmp_path / "parallel.bed"

    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_serial)],
        stranded=True, chroms=['chr21'],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=1,
    )
    ab.annotate_bed(
        chr21_db, [str(bed)], [str(out_parallel)],
        stranded=True, chroms=['chr21'],
        transcript_priority=[['protein_coding', 'distintron500']],
        gene_priority=[['protein_coding', 'distintron500']],
        gtf_format='gencode', append_chr=False, fuzzy=0, cores=2,
    )
    assert out_serial.read_text() == out_parallel.read_text()


def test_proxdist_intervals_match_legacy_shape():
    from annotator import annotation_functions as af

    intron_interval = Interval.from_list(
        ['chr1', '0', '10', 'intron', '0', '+']
    )
    proxdist_dict = af.get_proxdist_from_intron(intron_interval, distance=5)
    assert len(proxdist_dict['prox']) == 1
    assert proxdist_dict['prox'][0] == Interval(
        'chr1', 0, 10, 'proxintron5', '0', '+'
    )
    assert proxdist_dict['dist'] == []

    intron_interval = Interval.from_list(
        ['chr1', '0', '11', 'intron', '0', '+']
    )
    proxdist_dict = af.get_proxdist_from_intron(intron_interval, distance=5)
    assert len(proxdist_dict['prox']) == 2
    assert proxdist_dict['prox'][0] == Interval(
        'chr1', 0, 5, 'proxintron5', '0', '+'
    )
    assert proxdist_dict['prox'][1] == Interval(
        'chr1', 6, 11, 'proxintron5', '0', '+'
    )
    assert proxdist_dict['dist'][0] == Interval(
        'chr1', 5, 6, 'distintron5', '0', '+'
    )
