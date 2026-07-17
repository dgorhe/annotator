# Annotator

Annotator assigns genic annotations to BED intervals using a gffutils SQLite
database built from a GTF/GFF file. Each output row includes the original BED
fields plus assigned gene ID, gene name, genic region, region type, and a
pipe-delimited list of all overlapping annotations.

## Quick start

```bash
uv sync
uv run annotator \
  --input peaks.bed \
  --output peaks.annotated.bed \
  --gtfdb gencode.v19.annotation.gtf.db
```

## CLI tools

| Command | Description |
|---------|-------------|
| `annotator` | Annotate BED6 files with genic regions |
| `build_gffutils_db` | Build a gffutils database (optionally with `.annotator.idx`) |
| `create_region_bedfiles` | Export CDS, UTR, intron, gene, and exon BED files |
| `get_region_lengths` | Report average and total region lengths |
| `gene_name2id` | Map gene names to IDs in a tabular file |
| `miRNA_name2id` | Map miRNA names to precursor/mature IDs |
| `create_as_structure` | Write per-gene AS_STRUCTURE GFF entries |

## Documentation

See the [API Reference](api/index.md) for module-level documentation generated
from Google-style docstrings in the source code.

## Building docs locally

```bash
uv sync
uv run mkdocs serve    # live preview at http://127.0.0.1:8000
uv run mkdocs build    # static site in site/
```

## Agent-friendly documentation

Each `mkdocs build` also generates:

- `site/llms.txt` — curated index with links to Markdown versions of each page
- `site/llms-full.txt` — full documentation in one file for LLM context
- Per-page `.md` files alongside HTML (e.g. `site/api/annotate_bed/index.md`)

When serving locally, these are available at:

- http://127.0.0.1:8000/llms.txt
- http://127.0.0.1:8000/llms-full.txt
