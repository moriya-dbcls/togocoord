# TogoCoord

English | [日本語](README.ja.md)

TogoCoord converts sequence coordinates between the layers of life-science data: genomes, transcripts, proteins and 3D structures. A position is written as a **Location ID** based on the INSDC location syntax (e.g. `refseq:NC_000003.12:complement(49358132..49358134)`, `uniprot:P07203:49`) and is also available as FALDO JSON-LD.

- Built on primary repositories (GenBank / GFF3 + FASTA) so that it works for any species; rich resources for human and mouse (Ensembl, UniProt, SIFTS, MANE, ...) are added as extensions.
- A correspondence is a mapping ("location = mapping"), handled as operations on block lists (map, invert, compose).
- Sequences with identical residues (refget digest) are treated as the same, and exact paths are preferred.
- Assemblies of one species (e.g. GRCh37 ↔ GRCh38) are crossed automatically when needed; other species are reached only when requested.
- Every correspondence is validated against the actual sequences at ingest time.

## Layout

| Directory | Contents |
|---|---|
| `core/` | Location ID parsing, canonical form and FALDO JSON-LD; mappings as block lists (no dependencies) |
| `ingest/` | Adapters for GenBank / GFF3 / FASTA / SIFTS / MANE / UCSC chain / PAF / BED / NCBI assembly reports, self-validation, SQLite (R*Tree) stores, CLI `togocoord-ingest` |
| `service/` | Path search across stores, REST API and web UI (`togocoord-serve`) |
| `docs/` | Design, specifications and data; Japanese versions as `*.ja.md` ([design.md](docs/design.md), [data.md](docs/data.md), [spec-core.md](docs/spec-core.md), [spec-ingest.md](docs/spec-ingest.md), [spec-service.md](docs/spec-service.md), [scaling.md](docs/scaling.md)) |
| `scripts/` | Downloading, building and serving the data ([data.md](docs/data.md)) |
| `poc/` | Proof-of-concept implementation (web UI and SPARQList; for reference) |

## Usage

Node.js 24 or later (TypeScript runs directly through type stripping).

```sh
npm install
npm test            # tests of core / ingest / service

# Build stores (data is not in the repository; get it from NCBI and other sources)
node ingest/src/cli.ts --db human.sqlite --assembly-report GRCh38.p14_assembly_report.txt \
  --fasta GRCh38.p14_genomic.fna --fasta GRCh38.p14_protein.faa.gz GRCh38.p14_genomic.gff.gz
node ingest/src/cli.ts --db human_uniprot.sqlite UP000005640_9606.fasta.gz

# REST API and web UI
node service/src/serve.ts --port 8080 human.sqlite human_uniprot.sqlite
# http://127.0.0.1:8080/  ·  /v1/convert?loc=uniprot:P07203:49&to=genome
# Across species (with UCSC chains loaded): /v1/convert?loc=uniprot:P07203:49&to=protein&db=uniprot&taxon=10090
# Assembly sequence names as input (with GRCh37 and chains loaded): /v1/convert?loc=hg19:chr7:140453136&to=genome&assembly=GRCh38
```

The demo data (human, mouse, Arabidopsis, Marchantia, SIFTS, fanta.bio, chains and alignments) is downloaded and built with `node scripts/data.ts build` and served with `node scripts/data.ts serve`; see [docs/data.md](docs/data.md).

To serve it under a subdirectory behind a reverse proxy, see [spec-service §6](docs/spec-service.md).

Inputs and CLI options: [ingest/README.md](ingest/README.md). API: [docs/spec-service.md](docs/spec-service.md).

## License

MIT ([LICENSE](LICENSE)).
