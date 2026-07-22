# Database engine evaluation

Notes from evaluating whether annotator should leave local SQLite / gffutils
for MySQL, Postgres, or another store.

## Current query model

BED annotation does **not** issue a SQL range query per interval. The hot path:

1. Hydrates a spatial hash from `.annotator.idx` (or builds it from gffutils).
2. Looks up `features_dict[(chrom, position // 1e6, strand)]`.
3. Linear-scans that bucket for overlaps.

So replacing SQLite with MySQL/Postgres **without** changing this model does
not speed up annotation: the cost is Python object materialization and full
(or large) hydrate, not the SQLite query planner.

## What we prototyped instead

Keep SQLite `.annotator.idx` and query **by bucket** from disk:

| Mechanism | Where | Role |
|-----------|--------|------|
| Eager BED triage | `annotate_bed` + `load_annotator_index(..., buckets=...)` | Load only `(chrom, bucket)` pairs touched by the input BED |
| Chrom filter on `features` | `load_annotator_index` | Avoid `SELECT * FROM features` when chroms are limited |
| `query_bucket_features` | `index_store` | One-shot `(chrom, bucket, strand)` SQL join |
| `LazyBucketFeaturesDict` | `index_store` | LRU of on-demand bucket queries (`--lazy-buckets`) |

These give the low-memory / cold-start wins people often expect from a “better
database,” while staying compatible with gffutils and the existing hash API.

## Why not migrate gffutils to Postgres/MySQL

- gffutils’s FeatureDB is SQLite-centric; a server backend fights the library.
- Local batch jobs pay connection and serialization overhead for no gain under
  a load-then-hash architecture.
- Concurrent multi-client OLTP is not the annotator workload.

## When another engine would make sense

Only as a **new annotate backend** that answers per-interval overlap queries
from disk without hydrating the genome, for example:

- Postgres `int4range` / GiST (or `btree_gist`) on `(chrom, range)`
- Tabix / CSI over a sorted annotation BED/GTF
- Packed NCList or similar interval index

For large genome-wide BEDs that touch most buckets, an in-memory hash after a
single load can still win on wall time. Prefer measuring BED coverage (unique
1 Mb buckets) before investing in a second backend.

**Decision:** stay on SQLite sidecar indexes; invest in BED-driven chrom/bucket
triage and optional lazy bucket LRU rather than migrating engines.
