# API Reference

Documentation is generated automatically from Google-style docstrings in the
`annotator` package using [MkDocs](https://www.mkdocs.org/) and
[mkdocstrings](https://mkdocstrings.github.io/).

## Core annotation

- [`annotate_bed`](annotate_bed.md) — BED annotation pipeline, overlap lookup, and prioritization
- [`annotation_functions`](annotation_functions.md) — GTF key mapping, exon/CDS catalogs, UTR classification
- [`index_store`](index_store.md) — Pre-built `.annotator.idx` companion index
- [`intervals`](intervals.md) — Lightweight BED interval type

## CLI entry points

- [`annotator`](annotator.md) — Main `annotator` command-line interface
- [`build_gffutils_db`](build_gffutils_db.md) — Database and index builder

## Utilities

- [`create_region_bedfiles`](create_region_bedfiles.md) — Region BED export
- [`get_region_lengths`](get_region_lengths.md) — Region length statistics
- [`gene_name2id`](gene_name2id.md) — Gene name to ID conversion
- [`miRNA_name2id`](miRNA_name2id.md) — miRNA name to ID conversion
- [`create_AS_STRUCTURE`](create_AS_STRUCTURE.md) — AS_STRUCTURE GFF writer
