"""Lightweight interval type for prox/dist intron splitting."""

from dataclasses import dataclass


@dataclass(frozen=True)
class Interval:
    """A BED6-style genomic interval.

    Immutable value type used when splitting introns into proximal and
    distal regions during annotation.

    Attributes:
        chrom (str): Chromosome or contig name.
        start (int): Zero-based start coordinate (BED convention).
        end (int): One-based end coordinate (BED convention).
        name (str): Feature or region name.
        score (str): Score field from the source BED record.
        strand (str): Strand ('+' or '-').
    """

    chrom: str
    start: int
    end: int
    name: str
    score: str
    strand: str

    @classmethod
    def from_list(cls, fields):
        """Construct an Interval from a BED6 field list.

        Args:
            fields (list): Sequence of at least six fields in BED6 order:
                chrom, start, end, name, score, strand. Start and end are
                coerced to int.

        Returns:
            Interval: Parsed interval with typed coordinates.
        """
        return cls(
            chrom=fields[0],
            start=int(fields[1]),
            end=int(fields[2]),
            name=fields[3],
            score=fields[4],
            strand=fields[5],
        )

    def __eq__(self, other):
        """Compare intervals for equality across all BED6 fields.

        Args:
            other: Object to compare against.

        Returns:
            bool: True if other is an Interval with identical field values;
                NotImplemented if other is not an Interval.
        """
        if not isinstance(other, Interval):
            return NotImplemented
        return (
            self.chrom == other.chrom
            and self.start == other.start
            and self.end == other.end
            and self.name == other.name
            and self.score == other.score
            and self.strand == other.strand
        )
