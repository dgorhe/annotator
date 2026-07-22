# Annotator index format

This note describes the gffutils annotation database used by `annotator`, how
runtime search works today, and the companion `.annotator.idx` index that
speeds up repeated BED annotation runs.

## gffutils source database (`.gtf.db`)

Built by `build_gffutils_db` via `gffutils.create_db()`. Stored as SQLite.

| Table | Rows (GENCODE v50) | Purpose |
|-------|-------------------|---------|
| `features` | ~11.3M | All GTF intervals |
| `relations` | ~23.1M | Parent/child hierarchy (unused by annotate path) |
| `meta` | 1 | Dialect/version |

### `features` columns

| Column | Type | Notes |
|--------|------|-------|
| `id` | TEXT | Primary key (gffutils internal) |
| `seqid` | TEXT | Chromosome/contig |
| `featuretype` | TEXT | exon, CDS, UTR, transcript, gene, … |
| `start`, `end` | INT | 1-based closed GTF coordinates |
| `strand` | TEXT | `+`, `-`, `.` |
| `attributes` | TEXT | Semicolon-delimited GTF attributes (parsed at runtime) |
| `bin` | INT | UCSC bin field (unused by annotator) |

### Existing SQLite indexes (gffutils)

- `(seqid, start, end, strand)`
- `(seqid, start, end)`
- `(featuretype)`

These help gffutils range queries but **do not** match the annotator's 1 Mb
bucket lookup pattern.

## Runtime search (legacy path)

Each `annotator` invocation:

1. Scans all exons/transcripts/CDS via `features_of_type()` (full-table reads).
2. Iterates every feature on every chromosome via `db.region(seqid=chrom)`.
3. Builds an in-memory hash: `features_dict[(chrom, bucket, strand)] → [Feature]`
   where `bucket = position // 1_000_000`.
4. For each transcript, infers introns and materializes prox/dist intron
   intervals as cloned feature objects.
5. Annotates BED rows against the in-memory hash only (fast).

**Bottleneck:** step 2–4 dominate startup (~15 min on full GENCODE v50). Per-BED
annotation is comparatively cheap once the hash exists.

### Rows never used at annotation time

- **`gene` features** (~79k): bucketed but discarded because `annotate()`
  ignores features without a `transcript_id` attribute.
- **`relations` table**: never queried by the annotate CLI.
- Most **`transcript` rows**: used only to infer introns; the transcript
  interval itself may still overlap queries.

## Companion index (`.annotator.idx`)

A separate SQLite file built once alongside the gffutils DB. It stores exactly
what the annotate path needs, pre-bucketed:

| Table | Purpose |
|-------|---------|
| `meta` | `hash_val`, `gtf_format`, `append_chr`, `all_chromosomes`, `indexed_chroms`, source db path/mtime |
| `cds` | `transcript_id → (low, hi)` for UTR classification |
| `features` | Compact annotation records (integer PK) |
| `bucket_features` | `(chrom, bucket, strand, feature_id)` membership |

### Design choices

1. **Pre-bucket at build time** — eliminates per-run `hash_features()`.
2. **Drop `gene` rows** — they never contribute to output.
3. **Pre-compute prox/dist introns** — synthetic intervals stored as ordinary
   feature rows; no runtime cloning.
4. **Integer feature IDs + bucket table** — each logical interval stored once;
   buckets hold IDs only (avoids duplicating large objects in memory).
5. **Optional chromosome filter on load** — `--limit-chroms-to` loads a subset.
   Indexes record whether they cover all chromosomes or a fixed subset; a
   partial index cannot be used for full-genome annotation (runtime falls
   back to rebuilding from the gffutils DB).
6. **BED-driven triage** — by default, chromosomes are inferred from the input
   BED and only overlapping 1 Mb buckets are hydrated from the index. Use
   `--all-chroms` / `--no-triage-buckets` to disable. `--lazy-buckets` keeps an
   LRU of on-demand SQLite bucket queries (single-core).

Build (once):

```bash
build_gffutils_db \
  --annotation_file annotation.gtf.gz \
  --db_file annotation.gtf.db \
  --disable_infer_genes --disable_infer_transcripts \
  --build_annotator_index
```

This writes `annotation.gtf.db.annotator.idx`. Subsequent runs auto-detect and
load it when present (use `--no-index` to force the legacy rebuild path).

Expected tradeoff: index build adds time to the one-time DB prep step, but
repeat annotation runs should start in seconds instead of minutes.

## Slim annotate-only artifacts

Keep the full gffutils `.gtf.db` for utilities (`create_region_bedfiles`, etc.).
Ship or cache a smaller **annotate-only** `.annotator.idx` for BED runs:

### Rebuild index from an existing DB

```bash
build_gffutils_db --index-only --db_file annotation.gtf.db
```

### Chromosome-sharded indexes

```bash
build_gffutils_db --index-only --db_file annotation.gtf.db \
  --limit-chroms-to chr21 \
  --index-file annotation.gtf.db.chr21.annotator.idx
```

Point `--gtfdb` at the same `.gtf.db` and place/copy the shard as the default
sidecar (`annotation.gtf.db.annotator.idx`), or rebuild a shard matching the
BED’s contigs. Sharded indexes are smaller on disk and faster to load; they
cannot satisfy `--all-chroms` full-genome runs (annotator falls back to the
legacy rebuild).

The index already omits gene rows and raw GTF attribute blobs relative to the
source DB. Further feature-type pruning is not exposed: dropping types that
appear in priority lists would change annotation results.
