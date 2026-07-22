# annotator

Annotate BED6 intervals with genic regions from a [gffutils](https://pythonhosted.org/gffutils/) SQLite database. Each output row keeps the original BED fields and adds assigned gene ID, gene name, genic region, region type, and a pipe-delimited list of all overlapping annotations.

Requires **Python ≥ 3.10**.

## Install

```bash
git clone https://github.com/byee4/annotator/
cd annotator
uv sync
```

Or with pip:

```bash
pip install -e .
```

Optional extras: `pip install -e ".[regions,utils]"` for `pybedtools`/`pandas` tools.

## Quick start

```bash
uv run annotator \
  --input peaks.bed \
  --output peaks.annotated.bed \
  --gtfdb gencode.v19.annotation.gtf.db
```

If a companion `.annotator.idx` sits next to the gtfdb, it is loaded automatically for fast startup. By default, chromosomes and 1 Mb spatial buckets are inferred from the input BED so only relevant slices are loaded. Use `--all-chroms` / `--no-triage-buckets` to load more, `--lazy-buckets` for an on-demand LRU (single-core), or `--no-index` to force the legacy in-memory rebuild.

## Build a database (and optional index)

```bash
uv run build_gffutils_db \
  --annotation_file annotation.gtf.gz \
  --db_file annotation.gtf.db \
  --disable_infer_genes \
  --disable_infer_transcripts \
  --build_annotator_index
```

That writes `annotation.gtf.db` and `annotation.gtf.db.annotator.idx`. See [`annotator/INDEX.md`](annotator/INDEX.md) for the index format, slim chrom-sharded indexes (`--index-only`), and [`annotator/ENGINE.md`](annotator/ENGINE.md) for why the annotate path stays on SQLite sidecar indexes rather than MySQL/Postgres.

### Pre-built databases

Use at your own risk — GTF attribute layouts differ by source.

- [Human hg19 Gencode v19](https://s3-us-west-1.amazonaws.com/genome-references/gencode.v19.annotation.gtf.db)
- [Mouse mm10 Gencode vM10](https://s3-us-west-1.amazonaws.com/genome-references/gencode.vM10.annotation.db)
- [C. elegans WS257 extended](https://s3-us-west-1.amazonaws.com/genome-references/c_elegans.PRJNA13758.WS257.canonical_geneset.extend_utr.gtf.db) (canonical geneset plus upstream/downstream annotations)

These older downloads do not include a `.annotator.idx`; rebuild one with `build_gffutils_db --build_annotator_index` if you want fast startup.

## Output

Original BED6 columns, plus:

| Column | Description |
|--------|-------------|
| Assigned Gene ID | Chosen gene after prioritization |
| Assigned Gene Name | Gene name for the assigned gene |
| Genic Region | Chosen feature region (e.g. CDS, 3utr) |
| Genic Region Type | Transcript class (e.g. protein_coding) |
| All overlapping annotations | `\|`-separated overlaps |

Each overlap looks like:

```text
transcript_id:region_start:region_stop:strand:region:gene_id:gene_name:transcript_type:overlap
```

## Default priority (hg19 GENCODE-style)

Example list (see `datasets/hg19.gencode19.transcript.priority.txt`):

```text
miRNA,exon
protein_coding,CDS
protein_coding,start_codon
protein_coding,stop_codon
protein_coding,3utr
protein_coding,5utr
protein_coding,unclassified_utr
protein_coding,proxintron500
protein_coding,distintron500
non_coding,exon
protein_coding,Selenocysteine
non_coding,proxintron500
non_coding,distintron500
non_coding,transcript
non_coding,gene
non_coding,Selenocysteine
```

## Methods

For each query interval, annotator first picks the highest-priority overlapping transcript region *within each gene*, then picks among genes with the same priority rules. Example: if a feature overlaps a 3′UTR of Transcript A and a CDS of Transcript B (same gene), CDS wins. If multiple genes overlap, each gene gets a chosen transcript first, then genes are ranked.

## Options

| Flag | Description |
|------|-------------|
| `--transcript-priority-file` | Comma-delimited `transcript_type,feature` list (line order = priority). Unlisted features are appended after explicit ones. See `datasets/*priority.txt`. |
| `--gene-priority-file` | Same format; chooses which gene to report. |
| `--format` | GTF attribute nomenclature: `gencode` (default), `ensembl`, `ce10`/`ce11`/`wormbase`, etc. Replaces the old `--species` flag. |
| `--unstranded` | Allow unstranded matches; stranded BED features still prefer matching strand first. |
| `--append-chr` | If the BED uses `chrI`-style names but the DB expects `I`, prepend `chr` when querying. |
| `--limit-chroms-to` | Restrict the feature hash / index load to listed chromosomes. When omitted, contigs are inferred from the input BED. |
| `--all-chroms` | Load every chromosome instead of inferring from the BED. |
| `--no-triage-buckets` | With an index, load all buckets for selected chroms instead of only BED-overlapping 1 Mb buckets. |
| `--lazy-buckets` | On-demand LRU bucket loads from the index (single-core only). |
| `-c` / `--cores` | Parallel annotation workers (default: 1). |
| `--no-index` | Ignore a neighboring `.annotator.idx` and rebuild from the gtfdb. |

`--format` only selects attribute key names (`gene_id`, `gene_name`, …). For a new GTF dialect, extend `annotation_functions.get_keys()`.

## CLI tools

| Command | Description |
|---------|-------------|
| `annotator` | Annotate BED6 files with genic regions |
| `build_gffutils_db` | Build a gffutils DB (optionally with `.annotator.idx`) |
| `create_region_bedfiles` | Export CDS, UTR, intron, gene, and exon BEDs |
| `get_region_lengths` | Report average and total region lengths |
| `gene_name2id` | Map gene names to IDs in a table |
| `miRNA_name2id` | Map miRNA names to precursor/mature IDs |
| `create_as_structure` | Write per-gene AS_STRUCTURE GFF entries |

### `create_region_bedfiles`

```bash
uv run create_region_bedfiles \
  --db_file gencode.vM10.annotation.db \
  --species mm10 \
  --cds_out mm10_cds.bed \
  --proxintron_out mm10_proxintrons.bed \
  --distintron_out mm10_distintrons.bed \
  --utr5_out mm10_five_prime_utrs.bed \
  --utr3_out mm10_three_prime_utrs.bed
```

Only specify the region outputs you need. For clip_analysis_legacy / analyze_motifs naming, use `${SPECIES}_cds.bed`, `${SPECIES}_distintron500.bed`, `${SPECIES}_proxintron500.bed`, `${SPECIES}_three_prime_utrs.bed`, `${SPECIES}_five_prime_utrs.bed`, `${SPECIES}_exons.bed`, `${SPECIES}_genes.bed`.

Region construction (per gene, merged unless `--by-transcript`):

- **CDS** — merge overlapping CDS features
- **Prox/dist introns** — infer introns from exons; proximal ≤ 500 nt from an exon, distal otherwise
- **5′/3′ UTRs** — classify UTR features relative to CDS

### `miRNA_name2id`

Append accession IDs to a table of miRNA names using a mirBase gffdb and/or a custom `PRECURSOR_ACCESSION\tNAME[\tMATURE]` map:

```bash
uv run miRNA_name2id \
  --input all_mirnas.csv \
  --sep , \
  --name_col miRNA \
  --custom ensembl2name.tsv \
  --gffdb mmu.mirBase_v21.GCm38.gff3.db \
  --output all_mirnas.with_ids.csv
```

Multiple accessions for one name are joined with `|`. Use `--add_mature` for an extra mature-ID column. Build mirBase gffdbs with `build_gffutils_db`.

## Documentation

```bash
uv sync
uv run mkdocs serve   # http://127.0.0.1:8000
uv run mkdocs build   # site/
```

API docs live under `docs/`. After `mkdocs build`, agent-oriented dumps are at `site/llms.txt` and `site/llms-full.txt`.
