# @togocoord/ingest

English | [日本語](README.ja.md)

The TogoCoord adapters (v0.2). They read GenBank/GenPept flat files and GFF3, and output sequences, mapping edges and annotations. Self-validation is done during ingest (checking translations and checking the number of identical bases in alignments).

- Rules: [../docs/spec-ingest.md](../docs/spec-ingest.md)
- Core: [../core](../core)

## CLI

```sh
# Small file: JSON Lines to standard output
node ingest/src/cli.ts NC_012920.1.gb > mt.jsonl

# Whole genome: to SQLite (reads .gz directly; FASTA is read only for the needed ranges via .fai)
node --max-old-space-size=512 ingest/src/cli.ts \
  --db human.sqlite --fasta GRCh38.fna --fasta GRCh38_protein.faa.gz GRCh38_genomic.gff.gz
# stderr: sequences 137512, edges 505277 {"skipped":359831,"ok":142549,"mismatch":2897}, ... (73 s)
#         mismatches: 0 unexplained, 2897 with an INSDC /exception
```

| Option | Description |
|---|---|
| `--db FILE` | Save to SQLite (otherwise output JSON Lines). An existing file is replaced with `--overwrite` |
| `--fasta FILE` | Sequences used for self-validation (genome, transcripts, proteins). Can be given any number of times. Files over 64MB or with a `.fai` are read with random access |
| `--all-annotations` | Also save exons as annotations |
| `--from-report FILE`, `--to-report FILE` | NCBI assembly reports of the source and target assemblies of a UCSC chain file (`.chain(.gz)`). Translates names such as `chr1` to RefSeq accessions (spec-ingest §14) |

In the JSON Lines output, the `record` of each line is one of `sequence` / `edge` / `annotation` / `warning`.

## Store validation

```sh
node ingest/bench/verify-store.ts human.sqlite GRCh38.fna GRCh38_protein.faa.gz 10000
# {"residues":10000,"genomicHits":10616,"agree":10594,"disagree":0,"disagreeOnExceptionEdges":22,...,"roundTrip":10616,"roundTripMissing":0,...}
# protein -> genome: p50 0.06 ms, p95 0.62 ms, p99 0.95 ms
```

Converts randomly chosen protein residues to the genome using the store, translates the codons and compares them with the published sequence. It also converts back from the genome and checks that the original residue is returned. For how to read the results, see [../docs/scaling.md](../docs/scaling.md) §7.

## Library

```ts
import { ingestGenBank, edgeMapping } from "@togocoord/ingest";
import { createContext, mapLocation, parseLocationId, formatLocationId } from "@togocoord/core";

const result = ingestGenBank(text);
const nd6 = result.edges.find((e) => e.from === "refseq:YP_003024037.1")!;
const units = new Map(result.sequences.map((s) => [s.ref, s.unit]));
const ctx = createContext({ units: (ref) => units.get(ref) });
mapLocation(parseLocationId("refseq:YP_003024037.1:174", ctx), edgeMapping(nd6), ctx)
  .targets.map((t) => formatLocationId(t.location, ctx));
// => ["refseq:NC_012920.1:complement(14152..14154)"]
```

## Layout

| File | Description |
|---|---|
| `src/gbff.ts`, `src/gff3.ts`, `src/fasta.ts` | Parsers |
| `src/adapter-gbff.ts`, `src/adapter-gff3.ts` | Adapters (spec-ingest §3, §4) |
| `src/adapter-fasta.ts`, `src/adapter-sifts.ts`, `src/adapter-mane.ts`, `src/adapter-chain.ts` | Adapters for FASTA, SIFTS, MANE and UCSC chain (spec-ingest §9–§11, §14) |
| `src/validate.ts` | Self-validation (§6) and protein length estimation (§5) |
| `src/sequence.ts` | Sequence retrieval, reverse complement, translation, refget digests |
| `src/stream.ts` | Streaming reads (with `.gz` support), `JsonlSink` |
| `src/fasta-index.ts` | Creating `.fai` and random access to FASTA |
| `src/store.ts` | Writing to SQLite (`SqliteSink`) and queries (`TogoCoordStore`) |
| `bench/verify-store.ts` | Checks against published protein sequences and query time measurement |
| `src/genetic-codes.ts` | Translation tables generated from NCBI `gc.prt` |
| `test/fixtures/` | Real NCBI data (downloaded 2026-09-18) |
