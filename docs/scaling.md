# Scaling to Human Size (Prerequisite for Phase 3)

English | [日本語](scaling.ja.md)

2026-09-18. Measurements were taken on this development machine (macOS, Node 24.2).

---

## 1. Measurements

### 1.1 Data size (NCBI RefSeq GRCh38.p14, `GCF_000001405.40_GRCh38.p14_genomic.gff.gz`)

- 78MB compressed, **1.60GB and 4.93 million lines uncompressed**
- Line counts of the main features: exon 2.32 million, CDS 1.85 million (about 150,000 proteins), match 160,000, mRNA 145,000, biological_region 137,000, enhancer 115,000, cDNA_match 27,000
- **All lines with the same ID were adjacent** (across 1.73 million IDs spanning multiple lines, the maximum gap between lines was 1). On the other hand, there is only one `###` separator in the whole file.

### 1.2 Result of ingesting all of chromosome 3 with the current implementation (v0.1) (NCBI sviewer GFF3, 91MB, 300,000 lines)

| Item | Value |
|---|---|
| Processing time | 1.3–1.6 s |
| Memory (RSS) | **1.1–1.5 GB** (about 12–17 times the file size) |
| Output | 26,000 edges, 185,000 annotations |

- Before the fix, 662 warnings were produced. The cause was that `match` (RefSeqGene) with a reverse Target was not supported. After the fix, there were 0.

### 1.3 Projection to the whole genome

| Item | Projection | Assessment |
|---|---|---|
| CPU time | 20–30 s | No problem |
| Memory | About 20GB | ✗ |
| Reading | 1.6GB exceeds the V8 string limit (about 540 million characters) | ✗ Cannot be read |
| Genome sequence for self-validation (3.1GB) | Does not fit in memory | ✗ |

### 1.4 Performance of the candidate store (Node built-in `node:sqlite`, SQLite 3.50, with R*Tree)

Measured with synthetic data, storing 4 million blocks with two R*Trees, one for the src side and one for the tgt side.

| Operation | Value |
|---|---|
| Writing | **76 s** |
| Genomic interval lookup | **about 0.007 ms/query** |
| Protein-side lookup | about 0.045 ms/query |

---

## 2. Problems and countermeasures

| # | Problem | Countermeasure |
|---|---|---|
| P1 | The whole file is read as one string | **Read by streaming** (line by line; `.gz` is read while decompressing). GBFF is processed per record (`//`) |
| P2 | All lines, all features and all annotations are kept in memory | Change to **sequential output to a Sink**. For GFF3, use the approach "finalize the feature when the ID changes" (§3) |
| P3 | The duplicate check for proteins was quadratic | **Fixed** (O(1) check with a Set) |
| P4 | The whole genome sequence has to be loaded into memory for self-validation | **Random access from an indexed FASTA (`.fai`)**. Only the intervals needed for each CDS are read |
| P5 | Loading all edges for every query takes GB-scale memory and several seconds of startup with 4 million blocks | Retrieve only the blocks relevant to the query from a **store with an interval index** (SQLite and R*Tree) (§4) |
| P6 | JSON Lines output reaches several GB (mostly annotations) | Write directly to the store. By default, annotation attributes are limited to the ones needed |

---

## 3. Ingest design

```
.gff3(.gz) / .gbff(.gz)
  └→ read line by line → parser (GFF3: group by ID / GBFF: per record)
       └→ adapter (per-feature processing; current logic used as is)
            ├→ self-validation ← FaiSequenceSource (random access to genome FASTA via .fai)
            └→ Sink: sequence() / edge() / annotation() / warning()
                 ├→ MemorySink (current IngestResult; for tests and small data)
                 ├→ JsonlSink (current CLI output)
                 └→ SqliteSink (for large data; transactions per batch)
```

- **Grouping GFF3 features**: assuming lines with the same ID are adjacent, a feature is finalized when the ID changes. To handle general GFF3, finalization of the most recent N features (e.g. 1000) is deferred. If the same ID appears again after finalization, a warning is issued and counted, so that a broken assumption can be detected.
- **GBFF**: processed one record at a time per `//`. The sequence of a large chromosome record (250Mb each) is made readable from its position on disk, as with the `.fai` approach.
- **The adapters themselves are not changed**: only the places that append output to arrays are replaced with Sink calls. The current tests (matching against fixtures) confirm that the output is the same before and after.
- **Genome FASTA**: NCBI `.fna.gz` is plain gzip and does not allow random access. It is decompressed once at the start and a `.fai` is created (the same format as `samtools faidx`; the creation function is our own).
- **Building the R*Tree**: all rows are written first and the tree is built in bulk afterwards, aiming to shorten the 76 s.

Projection (to be measured after ingest): ingesting all of human in 1–2 minutes, with memory at 500MB or less.

---

## 4. Store and queries (foundation for phase 3)

```sql
sequence(id, ref UNIQUE, moltype, unit, length, topology, taxon, digest, provenance)
edge(id, kind, from_seq, to_seq, location, attributes JSON, provenance JSON, validation JSON)
block(id, edge, src_seq, src, tgt_seq, tgt, len, rev)
block_src  -- R*Tree (src_seq, src, src+len)
block_tgt  -- R*Tree (tgt_seq, tgt, tgt+len)
annotation(id, seq, type, location, attributes JSON)
annotation_idx  -- R*Tree (range of each segment)
```

**Query flow**
1. For each segment of the input Location ID, the overlapping blocks are retrieved with the R*Tree (`block_src` for the forward direction, `block_tgt` for the reverse direction).
2. A small `Mapping` is built from only the retrieved blocks and converted with the core's `mapLocation`. No change to the core is needed.
3. When going through multiple steps, blocks are retrieved and converted one step at a time. Path search uses only the list of edges (which sequences are connected to which).

DB access is projected at 1ms or less per step. The `Mapping` of frequently used edges is kept in an LRU cache.

**Projected size (human RefSeq)**

| Item | Projection |
|---|---|
| edge | About 400,000 (CDS 150,000, transcripts 200,000, alignments 190,000) |
| Blocks | About 4–5 million |
| DB | 1–2GB (with annotation attributes limited) |

---

## 5. Implementation order

1. The Sink abstraction and streaming parsers (with `.gz` support). Confirm with the current tests that the output does not change.
2. `FaiSequenceSource` and the `.fai` creation function.
3. `SqliteSink` and a read API (get blocks by interval, get edges by sequence). Add `--db` to the CLI.
4. **Measure with the human whole-genome GFF3 and genome FASTA** (ingest time, memory, DB size, self-validation results, query latency). Also check with mouse.

---

## 6. Decisions (2026-09-18)

- **Store**: `node:sqlite` (built into Node) was chosen. It requires no additional dependency packages, supports R*Tree, and can be distributed as a single file. However, it is experimental in Node 24 (stability 1.1).
- **Annotation scope**: exons are not stored (they can be derived from mRNA edges). They are stored if `--all-annotations` is given.

---

## 7. Measurements after implementation (2026-09-18)

Command: `togocoord-ingest --db OUT.sqlite --fasta genomic.fna --fasta protein.faa.gz genomic.gff.gz`

| | Human GRCh38.p14 | Mouse GRCm39 |
|---|---|---|
| Input | GFF3, 4.93 million lines (1.6GB) | GFF3, 3.11 million lines |
| Ingest time (including self-validation) | 70–73 s | 40 s |
| Peak memory (`--max-old-space-size=512`) | 0.77–1.1 GB (*) | 0.77 GB |
| DB size | 1.1 GB | 0.62 GB |
| edge / annotation | 505,000 / 714,000 | 279,000 / 372,000 |
| CDS self-validation | 142,549 matched. All 2,897 mismatches have an NCBI `/exception`. Unexplained mismatches: **0** | 97,327 matched. All 12 mismatches have an `/exception` |
| Creating `.fai` (3.3GB) | 2 s | — |

(*) Genome ingest (the streaming part) takes about 0.5GB, and the peak is in the final index-building stage (SQLite sorting and R*Tree). Decompressing the protein FASTA (gzip) into memory increases it further. The projected 500MB was not achieved, but as an amount used only once at build time, it was judged acceptable in practice.

**Check against published protein sequences** (`ingest/bench/verify-store.ts`. Protein residues are chosen at random, converted to the genome, and their codons are translated and compared with the published sequence. It also checks whether converting back from the genome returns to the original residue.)

| | Human (20,000 residues, 21,187 genomic mapped positions) | Mouse (10,000 residues) |
|---|---|---|
| Mismatches in CDSs without an exception | **0** | **0** |
| Mismatches in CDSs with an exception | 41 | 0 (1 with no genomic mapped position) |
| Rate of round trips returning to the original residue | **100%** | **100%** |
| Protein → genome (p50 / p99) | 0.05–0.06 ms / 0.5–1.0 ms | 0.04 / 0.12 ms |
| Genome → all mapped positions (p50 / p99) | 0.26–0.28 ms / 1.9–2.9 ms | 0.14 / 1.15 ms |

### 7.1 Issues found with real data and fixed

| Problem | Cause and fix |
|---|---|
| The heap grew to 2.5GB even with a Sink that keeps nothing | **V8 substrings kept holding the 1MB input chunk they were sliced from** (caused by long-lived reference keys). Long-lived keys are now copied (`ownString`). In addition, the readline async iterator was dropped in favor of line reading with a bound on buffer growth |
| `match` lines with Target `-` (RefSeqGene reverse to the chromosome) were not supported (662 on chromosome 3; 583 in the whole genome have a Gap) | Checked against real sequences and confirmed that they can be read with one rule: **Gap operations are written in the order the genome is traversed in the strand direction of the line** (spec-ingest §4) |
| `aa:Other` in `transl_except` (GFF3) could not be interpreted | Now interpreted case-insensitively (`OTHER` in GBFF) |
| A round trip of a codon split at an exon boundary gave `join(127..127c2,127c3)` | Pieces are now merged if contiguous on the target, even across input segments (spec-core v0.1.1 §5.4) |
| In intervals where multiple transcripts overlap, pieces to the same target were not merged | Merging is now done per target sequence. Property tests prevent regressions |
| Unrelated edges changed the truncation (`<`, `>`) decision | Truncation is now determined per target sequence (spec-core v0.1.1 §5.5) |
| `^` could not be converted for overlapping genes (ATP8 and ATP6 in different reading frames, etc.) | Now determined per target sequence (spec-core v0.1.1 §5.6) |
| Pseudogene CDSs (without protein_id) produced warnings | Now skipped without a warning, as with `/pseudo` in GBFF |

### 7.2 Remaining issues (phase 3 and later)

- **NCBI CDSs with `/exception`** (about 2,900 in human) have indels between the genome sequence and the protein or transcript, and cannot be mapped accurately from the genomic model. They need to be handled through the path "protein → RefSeq transcript (NM) → `cDNA_match` → genome". This requires ingesting RefSeq RNA GBFF (with CDSs on NM) and treating these edges as lower priority in path search (phase 3).
- **Self-validation of transcript and alignment edges** can be run if the transcript sequences (`rna.fna`) are given, but was not done this time (360,000 skipped).
- To bring **build-time memory** to 500MB or less, SQLite sort settings and a `.fai` for the protein FASTA are needed.
