# Downloading and building data

English | [日本語](data.ja.md)

How to rebuild the demo stores (SQLite) from the primary data. The definitions and the runner are in `scripts/data.ts`, and the serving order is in `scripts/stores.txt`. The data itself is not kept on GitHub (`data/` is excluded by `.gitignore`).

## Requirements

- Node.js 23.6 or later (with `npm install` run in the repository)
- For the Marchantia v3.1 ↔ v7.1 alignment only: [minimap2](https://github.com/lh3/minimap2) (built with 2.31; `brew install minimap2`)
- About 30GB of disk (about 18GB for downloaded files and uncompressed genomes, about 9GB for the stores). About 8GB of memory (up to about 6GB for the human build)

## Usage

```sh
node scripts/data.ts list                  # list datasets, groups, and whether they are built (size)
node scripts/data.ts download [name...]    # download the required files (all if no name is given)
node scripts/data.ts build [name...]       # build stores (missing ones only; --force rebuilds)
node scripts/data.ts serve [--port 8080] [--base URL] [--host 0.0.0.0]
                                           # serve the built stores in the order of scripts/stores.txt
```

- A name is a store (`human_rna`) or a group (`human`). `build` downloads the required files first if they are missing.
- The location is `$TOGOCOORD_DATA` (default: `data/` in the repository). `raw/` holds downloaded files (with the source's file names), `work/` holds intermediate files (alignment PAFs), `stores/` holds the stores, and `logs/` holds build logs.
- For large gzipped FASTA files, ingest creates an uncompressed file and a `.fai` in `raw/` (spec-ingest §8).
- Building everything takes about 13 minutes from already downloaded files (8-core Mac, peak memory about 7GB). Download time depends on the connection; the total is about 6GB.
- Example: build and serve human only. `node scripts/data.ts build human && node scripts/data.ts serve` (stores listed in `stores.txt` but not yet built are skipped).

## Datasets

| Group | Store | Input (source) | Size |
|---|---|---|---|
| human | `human` | RefSeq GRCh38.p14 (GCF_000001405.40) GFF3, genome, RNA and protein sequences, assembly report (NCBI) | 1.6GB |
| | `human_rna` | RefSeq GRCh38.p14 RNA GenBank | 650MB |
| | `human_ensembl` | Ensembl release 116 GFF3 and protein sequences | 3.1GB |
| | `human_uniprot` | UniProt reference proteome UP000005640 (current_release) | 57MB |
| | `human_mane` | MANE v1.5 summary and RNA sequences | 23MB |
| | `grch38_names` | GRCh38.p14 assembly report (sequence names) | 0.3MB |
| | `human_uniprot_alignments` | Alignments of UniProt entries without an identical sequence to RefSeq and Ensembl proteins (T2; UniProt `HUMAN_9606_idmapping_selected.tab.gz`) | 5MB |
| grch37 | `grch37` | GRCh37.p13 (GCF_000001405.25) assembly report and genome sequence | 0.3MB |
| | `chain_hg19ToHg38`, `chain_hg38ToHg19` | UCSC liftOver chain | 1MB, 10MB |
| chm13 | `chm13`, `chain_hs1ToHg38`, `chain_hg38ToHs1` | The RefSeq annotation of T2T-CHM13v2.0 (GCF_009914755.1) and the UCSC chains both ways | 1.6GB, 119MB, 128MB |
| mouse | `mouse`, `mouse_rna`, `mouse_uniprot`, `grcm39_names`, `mouse_uniprot_alignments` | RefSeq GRCm39 (GCF_000001635.27), UniProt UP000000589 | 930MB, 390MB, 21MB, 0.1MB |
| grcm38 | `grcm38`, `chain_mm10ToMm39`, `chain_mm39ToMm10` | GRCm38.p6 (GCF_000001635.26), UCSC liftOver chain | 0.2MB, 0.3MB, 0.5MB |
| human_mouse | `chain_hg38ToMm39`, `chain_mm39ToHg38` | UCSC liftOver chain | 158MB, 156MB |
| arabidopsis | `arabidopsis`, `arabidopsis_rna`, `arabidopsis_uniprot`, `arabidopsis_uniprot_alignments` | RefSeq TAIR10.1 (GCF_000001735.4), UniProt UP000006548 | 300MB, 174MB, 14MB |
| tair10 | `tair10`, `chain_tair10.1ToTair10`, `tair10_to_tair10.1` | Assembly report and genome sequence of the previous RefSeq version TAIR10 (GCF_000001735.3), UCSC GenArk chain (TAIR10.1 → TAIR10 only), minimap2 alignment for the reverse direction | 0.1MB, 0.1MB, 0.1MB |
| marchantia | `marchantia_v71` | INSDC MpTak_v7.1 (GCA_039105155.1) GenBank | 72MB |
| | `marchantia` | INSDC MpTak v3.1 (GCA_003032435.1) GenBank | 80MB |
| | `marchantia_uniprot` | UniProt UP000244005 | 7MB |
| | `mp_v31_to_v71`, `mp_v71_to_v31` | PAF from a minimap2 alignment of the two genome sequences (created in `work/`) | 3MB, 3MB |
| | `marchantia_tak2`, `marchantia_v51` | GenBank of INSDC MpTak2_v7.1 (GCA_037833965.1, Tak-2) and ASM993635v2 (GCA_009936355.2, v5.1) | |
| | `marchantia_cmv12` | Assembly report and genome of cmMarPoly1.2 (GCA_965642975.2; no annotation) | |
| | `marchantia_mpv4` | GenBank of INSDC Mp_v4 (GCA_001641455.1), the assembly UniProt's UP000077202 is built on | 66MB |
| | `marchantia_uniprot_v4` | UniProt UP000077202 (the Mp_v4 proteome) | 7MB |
| | `mp_tak2_to_hub`, `mp_hub_to_tak2`, `mp_v51_to_hub`, `mp_hub_to_v51`, `mp_cmv12_to_hub`, `mp_hub_to_cmv12`, `mp_v4_to_hub`, `mp_hub_to_v4` | minimap2 alignments with the hub MpTak_v7.1 (spec-ingest §20); `asm20` for cmMarPoly1.2 only | |
| sifts | `sifts` | SIFTS `uniprot_segments_observed.tsv.gz` (PDBe), `pdb_seqres.txt.gz` (wwPDB), all of the UniProt data above | 180MB |
| fanta | `fanta_human_hg38`, `fanta_mouse_mm10` | fanta.bio CRE v1.2.1 BED | 197MB, 116MB |

The results for each group (counts, self-validation, conversion checks) are in spec-ingest and spec-service.

**Rebuild check on 2026-09-18**: All 29 stores were rebuilt with this script, and the numbers of sequences, edges and annotations and the self-validation results were confirmed to match the earlier hand-built stores. The only differences are in the 4 chain stores and come from later improvements (the sequences of the assemblies on both sides are recorded, and 5 to 277 more chains are ingested by translating the UCSC names of patches). Among the added chains, those that land on unplaced scaffolds not included in RefSeq (such as `KI270752.1`) cannot be checked against the RefSeq genome sequence, so their self-validation is `skipped`.

**TAIR10 and TAIR10.1**: The nuclear chromosomes and the chloroplast have the same accessions and the same sequences; only the mitochondrial genome differs (TAIR10 `NC_001284.2`, TAIR10.1 `NC_037304.1`). A nuclear position is a valid answer as is in either assembly (spec-service §2.2).

**Not included**: Col-CEN v1.2 (a Col-0 assembly contiguous through the centromeres; UCSC GenArk has a chain with TAIR10.1) has no INSDC accession and is distributed only on GitHub (schatzlab/Col-CEN), so its sequence keys cannot be determined, and it is not included.

**On versions**: NCBI assemblies, Ensembl (release-116), MANE (release_1.5) and fanta.bio (v1.2.1) are downloaded from version-pinned URLs. UniProt (current_release), SIFTS and PDB sequences are updated by their sources at the same URL, so their content depends on when they were downloaded.

## Records for reproduction

The store's `meta` records the following (shown under Loaded data in the Web UI).

- MD5 of the input files (`inputs_md5`) and MD5 of the sequence files given with `--fasta` (`sequences_md5`)
- The TogoCoord commit used for the build (`togocoord`; `+local changes` if there are uncommitted changes)
- How the input was made (`method`; for alignments, the minimap2 version and arguments and how the PAF was filtered)

Running `build --force` with the same input (same MD5) and the same commit produces the same store.

## Serving and deployment

`serve` runs `togocoord-serve --stores scripts/stores.txt --store-dir data/stores --skip-missing`. The order in `stores.txt` affects the results (if several stores record the same sequence, the value from the earlier store is used; put annotated assemblies first and alignments last).

To deploy on another server, copy the repository and `data/stores/` (the stores are enough; `raw/` and `work/` are not needed), and start it as follows.

```sh
node --no-warnings service/src/serve.ts --host 127.0.0.1 --port 8080 --base https://example.org/togocoord/ \
  --stores scripts/stores.txt --store-dir /path/to/stores
```

To serve under a subdirectory, see spec-service §6. Passing `*.sqlite` all at once also loads unneeded files and uses alphabetical order, so use `--stores`.

## Adding a new dataset

Add the name, group and `togocoord-ingest` arguments to `STORE_LIST` in `scripts/data.ts`. Files written as `raw(URL)` (for NCBI assemblies, `ncbi(assembly(accession, name), suffix)`) become targets of `download`. To serve the store, add its name at the appropriate position in `scripts/stores.txt`.

## ID syntax from TogoID

`service/src/togoid-patterns.ts` holds the ID syntax of each database, so that an input written without a database can be read (spec-service §6). It is generated from [TogoID](https://togoid.dbcls.jp/)'s `config/dataset.yaml`; regenerate it when TogoID adds or changes a pattern:

```sh
node scripts/togoid-patterns.ts
npm test -w service   # the cases in service/test/infer.test.ts guard the table
```
