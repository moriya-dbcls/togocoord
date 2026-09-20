# TogoCoord service layer specification (v0.2: path search, identical sequences, REST API)

English | [日本語](spec-service.ja.md)

2026-09-18. Rules for phases 3a to 3c and 3e (SIFTS). Implemented in `service/` (`@togocoord/service`).

---

## 1. Multiple stores (`StoreSet`)

- Multiple SQLite stores (one per species or per data source; for example, the human genome GFF3 and the RefSeq RNA GenBank) are handled as one set.
- Blocks are fetched from all stores using R*Tree, taking **only those that overlap the location being converted**. An edge is identified by `<store number>:<edge ID>`.
- Each store's caches (reference key to ID, units, edges) are bounded LRU caches (50,000 entries by default). The SQLite page cache is 8MB by default; reads beyond that are left to the OS file cache.
- A store is opened only if its schema version (`meta.schema`) matches.
- Sequence information is merged across stores. For each field, the first value found is used; a length of 0 means "unknown" and is replaced by a later value; tags are combined as a union. The merged result and the list of identical sequences are kept in a bounded cache.
- The context sequence length is the length from the store (0 means unknown). Therefore, input beyond the sequence length is an error.

## 2. Conversion (`convert`)

Input: a Location, `to` (the target), and `maxHops`.

| `to` | Meaning |
|---|---|
| `{ ref }` | A specific sequence |
| `{ category }` | `genome` / `gene_region` / `transcript` / `protein` / `structure` (determined from accession rules and molecule type; PDB chains are `structure`) |
| `{ namespace }` | A namespace (for example, `uniprot`). In the API and UI, the recommended form is a layer (such as `protein`) combined with the result filter `db` |
| Omitted | All sequences within `maxHops` (default 1) |

`prefer` (a list of tags; the API default is `MANE Select` and `MANE Plus Clinical`): among paths of the same cost, choose the path that passes through tagged sequences (the number of untagged intermediate sequences passed through is the next comparison key after cost). Results of the same cost list tagged targets first. The cost value itself is not changed. Example: from a genomic codon of GPX1 to PDB, there are several Ensembl proteins identical to UniProt, but the path through the MANE Select ENSP00000407375.1 is chosen.

**Search**: a cost-weighted shortest path search (Dijkstra's algorithm) with sequences as nodes and edges as graph edges.
- A state is the pair "a sequence and a location on it". Expansion follows only edges that have blocks overlapping the location, and applies the core `mapLocation` with that edge's blocks.
- Each sequence is reached only once, by the cheapest path.
- Once a sequence matching the target condition is reached, from there the search follows **only identical sequences and near-identical UniProt records linked by a T2 alignment**. This is so that records of other databases for the same residue (RefSeq, Ensembl, UniProt) are also returned as results; choosing "protein" includes UniProt. However, a sequence identical to the input (reached only through identity steps) is treated as the input itself, and even if it matches the target type, the search continues from it as usual (example: `uniprot:P07203:49` → identical NP_000572.2 → genome → another isoform).
- Within one path, the same edge is not used again.
- The default `maxHops` is 4 when a target is given and 1 when not.
- There is no early cutoff (the search continues until the neighborhood up to the maximum number of hops is fully explored). Therefore, the search scope is narrowed by the layer rules in §2.1.

### 2.1 Layer rules (limiting the search scope)

Sequence types have a loose hierarchy: **genome(0) − gene region(1) − transcript(2) − protein(3) − structure(4)**. A conversion between conceptually distant layers only needs one U-turn at some layer, so paths that go deep up and down are not searched. The rules are set, arbitrarily, as follows.

1. **Depth limit**: do not enter a layer below the lower (closer to structure) of the start and the target. When the target is given only as a namespace (the layer is unknown), there is no limit.
2. **At most one U-turn**: switching between going up (toward the genome) and going down happens at most once. Moves within the same layer (identity, alignments within the same layer) and moves to or from sequences of unknown layer (`other`) are not counted.
3. The same sequence in different states (last direction and number of U-turns) is treated as a different node. As targets, only the cheapest result per sequence is returned.

Example: paths kept are ATP8 → genome → ATP6 (up then down), UniProt → Ensembl → genome, and genome → protein → PDB. Paths cut are UniProt → PDB → UniProt → genome (goes down to structure even though the target is the genome) and protein A → genome → protein B → another genome (zigzag).

Effect (human, UniProt residue 100 → genome): the number of expanded sequences went from 492 to 27 for p53, and to 10 for hemoglobin α.

**Orthologs**: edges between orthologous proteins are treated as species-crossing edges under the scope rules in §2.2. No exception to layer rule 1 is made.

### 2.2 Scope rules (species and assembly)

A target is specified on two axes: the "sequence type" (`to`) and the "scope" (species `taxon`, assembly `assembly`). Species and assemblies are handled differently. Crossing species is a homology-based correspondence that the user should choose. Between assemblies of the same species (GRCh37 ↔ GRCh38), the difference is only in the coordinate system of the same DNA, and annotations often exist on only one of them (GRCh38).

1. **By default, the species stays the same as the input.** A step into another species (a liftOver chain, an identical sequence of another species) is taken only when a `taxon` different from the input is given, and only once, into that species.
2. **Assemblies of the same species are crossed when the path needs them, even without being specified, up to twice per species.** Alignments are made towards the annotated assembly of a species (spec-ingest §20), so two other assemblies are joined through it (Marchantia cmMarPoly1.2 → MpTak_v7.1 → v3.1). After crossing species, the search can also cross that species' assemblies (mm10 → mm39 → hg38 → hg19). Example: a GRCh37 position → chain → GRCh38 → CDS → protein. The path keeps the `liftover` step, so you can see where the assembly changed.
3. **`assembly` only determines in which assembly genome results are returned.** The default is the input's assembly if the input is on a genome of that species, otherwise the species' default assembly (the one with the most annotations; GRCh38 for human). When no target is given (directly connected sequences), results are not filtered unless specified. To move a GRCh37 position to GRCh38, use `to=genome&assembly=GRCh38`. `assembly` can be `GRCh37.p13`, `GRCh37`, or the UCSC name `hg19`.
4. **The current species and assembly are carried along the path.** Results are judged by the carried species, so sequences without their own taxon record (such as Ensembl) are handled correctly. A sequence's species is determined, in order, from the sequence's own record, the record of the store containing it (`--taxon`), and, if its identical sequences belong to a single species, that species. A genome sequence's assembly is determined from the `assembly` of the store containing it, or from an assembly report that lists the sequence (spec-ingest §15). A chain store records the species and lengths of the sequences of the assemblies on both sides. Therefore, species-crossing steps can be recognized even if the other species' store is not loaded.
5. **From an out-of-scope target (for example, a human protein when mouse is specified, or a GRCh38 position when GRCh37 is specified), the search follows only scope-crossing steps, steps toward the genome (before the U-turn), and identical sequences.** It does not spread sideways within the same species.
**Species boundaries**: within the same species (NCBI species rank) the difference is an "assembly difference"; between different species it is a "species difference". Genomes of subspecies, varieties, strains and individuals are treated as other assemblies of the same species (this is close to a difference in genome coordinate systems; the same holds for each haplotype of a pangenome). A taxon is merged into a species only when this is explicit: that is, when `--species-taxon` is recorded in the store, or when the NCBI scientific name explicitly has a rank below species after the binomial (`subsp.`, `var.`, `f.`, `str.`, `strain`, `substr.`, `serovar`, `biovar`, `pv.`, `cv.`) and that binomial species is loaded. Example: *Marchantia polymorpha* subsp. *ruderalis* (1480154, MpTak_v7.1) is merged into *Marchantia polymorpha* (3197, v3.1 and UniProt). Sharing only the first two words of the scientific name is not enough to merge (Human immunodeficiency virus 1 and 2 are different species).

**Default assembly**: among assemblies with annotations, the one with the newest release date (Date in the assembly report), then the one with the most annotations (MpTak_v7.1 is newer than Marchanta_polymorpha_v1; GRCh37 has no annotations, so GRCh38).

**Between assemblies or species without a chain**: when converting from a genome to a genome of an assembly or species different from the input, the protein layer is allowed as an exception to layer rule 1. Through identical proteins (genome → CDS → identical sequence → CDS → genome), only the coding regions of genes whose sequences did not change can be converted.

6. **Hop limit**: when crossing species, `CROSSING_HOPS` (2) is added. For a species with multiple assemblies, 1 is added when a target is given. When species are crossed without a target, the search stops on arrival.

**Assembly names in input**: `<assembly>:<sequence name>[:<position>]` (for example, `hg19:chr7:140453136`, `GRCh37:7:140453136`) is also accepted. The sequence name can be the Sequence-Name in the assembly report, the UCSC name, or the GenBank accession. The API rewrites it as `refseq:NC_000007.13:140453136` and returns the canonical form in the response. `/v1/location` returns the written name (`written`), the species and the assembly.

The only kind of species-crossing edge today is `liftover`. When ortholog correspondences are added, they are handled as edges of the same kind.

**Orientation**: for nucleotide results, it is the strand of the result's interval. For protein results (always written from N-terminus to C-terminus), it is `reverse` when the input corresponds to the reverse strand of the coding sequence. Each path step returns the output strand as is, so orientations must not be multiplied across steps. This used to be done, so "protein → minus-strand genome → another protein of the same gene" was shown as `reverse` (fixed 2026-09-18).

## 3. Edge costs and priority

If an exact path exists, it is always chosen. There is no code that decides which database to go through case by case (such as via Ensembl for UniProt); it is decided only by the following cost order.

| Step | Cost |
|---|---|
| **Identity** (refget digests match; §3.1) | 0 |
| Primary-data edge (CDS, `cDNA_match`, SIFTS, etc.), validated or not subject to validation | 1 |
| Whole-genome alignment (`liftover`, UCSC chain; §9) | 2 |
| Self-computed protein alignment (T2; spec-ingest §18) | 3 |
| Mismatch in self-validation (`mismatch`) | 1 + 10 |
| Has an NCBI `/exception` and is not ok by full sequence comparison (`basis: "full"`) | 1 + 10 |

### 3.1 Identical sequences (identity)

- At ingest, the sequence's refget digest (derived from SHA-512, used to judge identity) and MD5 (a key for matching external data such as UniParc; UniParc stores it in uppercase, but comparison is case-insensitive) are recorded. This applies to GBFF sequences and `/translation`, published protein sequences given with `--fasta`, and UniProt and Ensembl sequences ingested with the FASTA adapter.
- A location whose positions are all on one specific sequence **can move at cost 0, keeping the same coordinates, to another sequence with the same digest**. In the path it is recorded as a `kind: "identity"` step.
- UniParc's CRC64 is not used to judge identity, because collisions are known.
- Example: UniProt P07203 (GPX1) and RefSeq NP_000572.2 have matching digests (all 203 residues are identical, including the selenocysteine at 49). Therefore, `uniprot:P07203:49` is converted to `refseq:NM_000581.4:220..222` via an identity step and a CDS step (cost 1).

- **If verified by sequence, the exception flag is ignored.** For example, a transcript that differs from the genome only by base substitutions keeps its coordinates. On the other hand, an ok with `basis: "partial"` (which only checked that there is no stop codon) is not considered verified.
- **Alignments take priority**: when an alignment edge (such as `cDNA_match`) linking the same two sequences overlaps the location, the annotation edge built from the genomic model is not used. The alignment represents the transcript's own sequence, but the genomic model does not.
- **`approximate`**: set on conversion results whose path includes an edge subject to the cost penalties above. The position may be shifted due to insertions or deletions.

## 4. Results

Each result has the following.

- The target Location and its ID
- The target type (category)
- The cost
- `approximate`
- The orientation
- The path (`path`)

Each path step has the following.

- The edge used, and whether it was used forward or reverse
- The location that entered the step
- The **parts that could not be mapped** at that step
- The edge's attributes, validation result and provenance

---

## 5. Validation (2026-09-18, human GRCh38.p14)

**Stores**
- Genome GFF3 ingested with self-validation against the genome, RefSeq RNA and RefSeq protein sequences: 120 seconds, 1.1GB
- RefSeq RNA GenBank: 33 seconds, 0.6GB. All 136,794 CDS edges reproduced the translation from each mRNA's own sequence

**Method**: `service/bench/verify-search.ts`. Residues of published proteins are chosen at random and converted to the genome. Judged in two ways:

- **Single-residue check**: whether the target codon encodes the residue of the published sequence.
- **15-residue window check**: convert the window and check whether its translation matches the published window at 80% or more. Base substitutions barely lower the window match rate, but a coordinate shift makes the whole window mismatch, so this checks "whether the coordinates are correct".

**Set 1: proteins with CDSs that NCBI marked with `/exception` (2,000 residues from 1,543 proteins)**

| Path type / target | Single-residue check (match / mismatch) | Window check (correct / shifted) |
|---|---|---|
| Direct edge (verified by sequence) / primary chromosomes | 1,822 / 0 | 1,767 / 0 |
| Via transcript (NM → `cDNA_match`) / primary chromosomes | 129 / 1 | 125 / **0** |
| Direct edge / alternate sequences (NT_, NW_) | 1,789 / 0 | 1,758 / 0 |
| Via transcript / alternate sequences | 3,022 / 28 | 2,934 / 15 (*) |
| Edge with exception (last resort, `approximate`) | 205 / 159 | 185 / 151 |

* All are alternate haplotypes of the MHC region. Both steps of the path are verified by sequence, so these are thought to be amino acid differences between haplotypes, not coordinate shifts.

**Set 2: all proteins (5,000 residues)**

| Path type | Single-residue check | Window check |
|---|---|---|
| Direct edge | match 5,178 / mismatch 0 | shifted 0 |
| Via transcript | match 52 / mismatch 0 | shifted 0 |

**Query speed** (conversion to the genome): p50 0.3ms, p99 1 to 2.4ms

### 5.1 Problems found and fixed during validation

| Problem | Fix |
|---|---|
| **An mRNA in NCBI's GFF3 is a single line covering the whole gene range**, with exons written separately as child features. This was used as is for the transcript edge, so it was mapped to the unspliced range. Found only after validating against transcript sequences (190,000 mismatches) | Transcript-type features build their location from the exons linked by `Parent` (spec-ingest §4) |
| The **poly(A) tail** of RefSeq transcripts (not in the genome) caused them to be judged as length mismatches (4,697 cases) | If the transcript is longer than the model and the overhanging part is 90% or more A, it is treated as a poly(A) tail and comparison stops before it |
| For transcripts whose genomic model does not match the NM, the model's edge was chosen as the cheaper path | Transcript edges are also validated by sequence (ok if only base substitutions; mismatch if lengths differ). Together with this, the cost rules (§3) and the alignment priority rule were introduced |
| A validation ok had no distinction of strength | `basis` (`full` or `partial`) is recorded (the store schema was raised to version 2) |

### 5.2 Remaining issues

- Many alternate sequences (NT_, NW_) have no `cDNA_match`. For transcripts with exceptions, an `approximate` path is used as a last resort.
- The choice between a path that maps only part of the input and a path that maps all of it at a higher cost is decided only by per-sequence cost. The fraction of the input that was mapped is not used in the choice.

---

## 6. REST API (`togocoord-serve`)

`node service/src/serve.ts [--port 8080] [--host 127.0.0.1] [--base URL] STORE.sqlite...`. Implemented with Node's built-in http, with no dependency packages. CORS is fully open.

| Method | Path | Description |
|---|---|---|
| GET | `/v1/convert?loc=&to=&db=&taxon=&assembly=&maxHops=&codon=never&tag=` | Conversion. `db` (namespace, for example `uniprot`; multiple allowed) filters results to that database (unknown namespace is 400). `taxon` (NCBI taxon number, `taxon:10090`, or the scientific or common name of a loaded species, for example `Mus musculus`, `mouse`) and `assembly` (assembly for genome results, for example `GRCh37`, `hg19`) specify the target scope (§2.2). `loc` can also be written with an assembly sequence name (`hg19:chr7:140453136`). It can also be the ID of an annotation that can be looked up by ID (`fanta:FCHS_301358` is the region of that CRE, `refseq:NC_000003.12:181712289..181712497`; `/v1/location` returns the ID, type, name and link in `written.annotation`). Unknown species or assemblies are 400. `to` is a type (such as `genome`), a namespace (such as `uniprot`) or a sequence (`refseq:NC_000001.11`), and can be given multiple times. If omitted, all directly connected sequences are returned. With `tag` (for example `MANE Select`), only targets with that tag are returned. Results carry the target's tags (`tags`), species (`taxon`, `organism`) and, for genomes, the assembly (`assembly`); the response carries the input's species and assembly (`inputTaxon`, `inputAssembly`). The UI shows the species name on results whose species differs from the input |
| POST | `/v1/convert` | Batch conversion. `{"locations": [...], "to": ..., "db": ..., "taxon": ..., "assembly": ..., "maxHops": ..., "codon": ...}`. Up to 1000 items. Errors in individual inputs are returned as `error` on that element |
| GET | `/v1/location?loc=` | Canonical ID, IRI, segments (1-based; for proteins, the residue number and the position within the codon), species and assembly, and, if written with an assembly sequence name, that name (`written`) |
| GET | `/v1/location/faldo?loc=` | FALDO JSON-LD (`application/ld+json`) |
| GET | `/v1/sequences/{ref}` | Sequence information (merged from all stores) and the list of identical sequences |
| GET | `/v1/sequences/{ref}/edges` | Edges leaving and entering the sequence |
| GET | `/v1/annotations?loc=` | Annotations overlapping the interval. The index is looked up by the feature's outer bounds, but only features whose own intervals overlap are returned (at an intron position, that transcript is not returned). For genomic positions, **annotations on other assemblies of the same species are also returned**. The position is moved to that assembly by a chain or alignment (one step only) and looked up there, and each annotation gets `assembly`, the looked-up position (`via`) and the position moved back to the input assembly (`lifted`). Example: at an mm39 position, fanta.bio CREs that exist only on mm10 are visible. Annotations with IDs get `id` (`fanta:FCMM_194523`) and `link` |
| GET | `/v1/meta` | Per store: file name, metadata (display name, species, assembly, source files, build time), content summary (spec-ingest §12). `species` lists the loaded species (taxon, name, assemblies, default assembly), `assemblies` lists assemblies (name, UCSC name, accession), and `annotationNamespaces` lists the namespaces of annotations that can be looked up by ID |
| GET | `/<namespace>:<accession>:<location>` | IRI resolution (identifiers.org-like). Depending on `Accept`: a 303 to the Web UI (`/?loc=`) for browsers, FALDO JSON-LD for JSON-LD requests, and a 303 to `/v1/location` for JSON requests |
| GET | `/`, `/ui/*` | Web UI (§6.1) |

- `loc` may be just `namespace:accession`, which is treated as the whole sequence (`1..length`). The response's `input` returns the explicit range.
- Limits: the total interval length of the input is up to 5 million bases (or residues) (413 if exceeded). Results are up to 1000 per conversion (`truncated: true` if exceeded).
- An IRI without a location (`/refseq:NC_000001.11`) represents the sequence itself, not a range. It returns a 303 to `/v1/sequences/{ref}` for JSON requests and a 303 to the UI for browsers.

Errors are returned as `{"error": ..., "position"?: ...}` (400: syntax or semantic errors, 404, 405, 413).

**Deploying in a subdirectory**: the UI loads its styles and scripts with relative paths (`ui/…`) and derives the API URL from its own script URL. Therefore, it works as is when placed in a subdirectory behind a reverse proxy (for example, `https://example.org/togocoord/` → `http://127.0.0.1:8080/`). Configure the URL without the trailing `/` (`/togocoord`) to redirect to `/togocoord/` (nginx's `location /togocoord/ { proxy_pass http://127.0.0.1:8080/; }` does this automatically). Match the IRI prefix with `--base https://example.org/togocoord/`.

### 6.1 Web UI (3d)

`service/public/` (HTML, JavaScript and CSS with no build step; it keeps the colors and fonts of the concept version). A thin client that uses only the REST API. State is kept in the URL (`?loc=&to=&db=&taxon=&assembly=&codon=`), so the page URL can be shared as is.

- Enter a Location ID and choose the target (directly connected sequences, genome, gene region, transcript, protein, structure), the database (default "any database"; filters results) and the species (default "same species"; chosen from the loaded species). The assembly selector is shown only when multiple assemblies are loaded for one species (default "default assembly"; options also show the UCSC name, as in `GRCh37.p13 / hg19`, and the species, accession and whether annotations exist are shown on mouse hover). The input card shows the species, the assembly and, if written with an assembly sequence name, that name (`written as chr7`). Genome results get the assembly name when the species has multiple assemblies. Examples for data that is not loaded (such as GRCh37) are not shown. A few representative examples are shown and the rest appear under "more examples" (they are all kept, to show the range of conversions).
- **Options that have no effect are disabled** (the reason is shown on mouse hover). Databases not in the target layer (such as UniProt for genome; structure has only PDB, so the whole selector is disabled). The assembly can be chosen only when the target is genome or directly connected sequences and the relevant species has multiple assemblies, and the options are only that species' assemblies. MANE Select only is available only when the target is transcript or protein and the relevant species has MANE markers. Only species that have a genome alignment (liftOver chain) from the input species can be chosen; others cannot (the reason is shown in a tooltip; they may sometimes be reached only through identical sequences, but this is rare and would look as if conversion were possible; the API's `taxon` accepts any species). If no other species can be chosen, the whole selector is disabled. The input species itself is the same as "same species", so it is not listed. A disabled option is reset to its default, and the conversion uses that value. These checks use `crossings` (species and assemblies linked by chains) and `tags` (species of sequences with tags) from `/v1/meta`. The display of the position within the codon can be toggled.
- Input: canonical ID, each segment (residue number and position within the codon), FALDO JSON-LD, overlapping annotations.
- Results: target ID, type, cost, `approximate` (when the path went through an edge that could not be verified), orientation, path (edge kind and direction; mouse hover shows the edge's location, validation result, exception and provenance), and the parts that could not be mapped at each step.
- Clicking a result ID re-queries "directly connected sequences" from that position (equivalent to the concept version's operation of continuing the conversion).
- Annotations are sorted from the narrowest range, excluding features that cover the whole sequence (chromosome, region, etc.).
- Input errors are shown with `^` under the character at the error position.
- Target tags (such as MANE Select) are shown as markers, and results can be filtered with "MANE Select only".
- **Loaded data** (`?view=data`): shows the loaded stores grouped by species (taxon). Each species is collapsed (expand to see each store), and the heading shows the number of datasets, the taxon and the assemblies. Shown are the display name, assembly, source files, counts and build time. Clicking a per-store example tries the conversion on the spot.
- Display and operation were checked in Chrome with the full human stores (RefSeq, RefSeq RNA, Ensembl, UniProt, SIFTS).

## 7. FALDO JSON-LD

Only vocabulary confirmed in the FALDO definition (`faldo.ttl`) is used. `faldo:NegativeStrand`, `faldo:member` and `faldo:order`, used in the concept version, are not in the FALDO definition.

| Location | FALDO |
|---|---|
| Single base / single residue | `faldo:ExactPosition` |
| Range | `faldo:Region` (`faldo:begin` and `faldo:end`. On the reverse strand, the biological start is begin, so begin has the larger number. Confirmed to match the cheY example in the FALDO README) |
| Strand | `faldo:ForwardStrandPosition` / `faldo:ReverseStrandPosition` as the position type (not added for proteins) |
| `a^b` | `faldo:InBetweenPosition` (`faldo:after` / `faldo:before`, following the strand) |
| `a.b` | `faldo:InRangePosition` |
| `<`, `>` | The position at that end gets type `faldo:FuzzyPosition` (`faldo:position` is kept) |
| `join` / `order` | `faldo:ListOfRegions` (`rdf:Seq`) / `faldo:BagOfRegions` (`rdf:Bag`). Elements are ordered `rdf:_1`, `rdf:_2`… in biological order |
| Codon extension | Positions get `tgc:codonPosition` (1..3). If adopted by FALDO, it will move to `faldo:` |

- The top-level node is the location itself (`@id` is the location IRI). In IRIs, only `<`, `>` and `^` are percent-encoded.
- Sequences are referenced as `https://identifiers.org/<namespace>:<accession>`.
- Expanded to RDF with jsonld.js and confirmed to produce correct triples.
- The default base IRI is `https://togocoord.dbcls.jp/` (provisional; changed from `togocoord.example.org` on 2026-09-19; can be overridden with `--base`). The vocabulary (`tgc:`) is `https://togocoord.dbcls.jp/ontology#`.

---

## 8. Mapping to structures (3e, 2026-09-18, human)

**Stores** (all built with `togocoord-ingest`)

| Store | Input | Time | Self-validation |
|---|---|---|---|
| RefSeq genome | GRCh38.p14 GFF3, genome, RNA and protein sequences | 134 seconds | 0 unexplained mismatches in CDS and alignments (excluding 5,778 with exceptions) |
| RefSeq RNA | RNA GenBank | 37 seconds | All 136,794 CDSs ok |
| Ensembl | release 116 GFF3, Ensembl protein sequences, `--seqid-map` (NCBI assembly report) | 190 seconds | **All ~370,000 CDSs ok** |
| UniProt | Human reference proteome (20,652 canonical + 148,999 additional) | 2 seconds | — |
| SIFTS | `uniprot_segments_observed` (limited to human UniProt) | 9 seconds | 240,124 ok out of 240,000 edges |

**Fraction with identical sequences** (UniProt sequences that have an identical Ensembl or RefSeq protein)

| UniProt | Count | Ensembl | RefSeq | Either |
|---|---|---|---|---|
| canonical, with structure | 8,983 | 98.1% | 97.5% | **98.6%** |
| canonical, all | 11,669 | 89.9% | 86.3% | 90.5% |
| additional (isoforms, etc.) | 148,896 | 93.7% | 17.5% | 94.5% |

98.6% of UniProt entries with structures can reach the genome exactly without alignment (identity → CDS). The rest need T2 (self-computed alignment).

**Round trip with structures** (`service/bench/verify-structure.ts`, 3,000 residues from SIFTS intervals)
- Reached the genome: 2,957 (all by exact paths; 0 approximate). The 43 without identical sequences cannot be reached.
- The target codon encodes the UniProt residue: **2,957 / 2,957**
- Back from the genome to the same PDB residue: **2,957 / 2,957** (for structures where the same residue appears twice in the same chain, both positions are returned)
- The PDB SEQRES residue matches the UniProt residue: 2,944 (the rest are engineered mutations in the structure)

**Example**: the GPX1 selenocysteine codon `refseq:NC_000003.12:complement(49358132..49358134)` → Ensembl CDS (reverse) → UniProt P07203 (identity) → SIFTS → `pdb:2F8A.A:59` and `pdb:2F8A.B:59` (2F8A is the U49G mutant). Cost 2, not approximate.

**Speed**: UniProt residue → genome is p50 6.4ms, p95 25ms, p99 50ms (proteins with many structures are more likely to be sampled, so these are conservative values). Before the layer rules (§2.1), it was p50 12.5ms, p99 175ms, expanding all of the hundreds to thousands of PDB chains at the same cost level as the target. RefSeq protein → genome is p50 1.3ms, p99 12ms.

## 9. Mapping across species (UCSC liftOver chain, 2026-09-18, human ↔ mouse)

**Stores**: UCSC `hg38ToMm39.over.chain.gz` and `mm39ToHg38.over.chain.gz` were each ingested into a separate store (spec-ingest §14).

| Store | Chains | Blocks | Size | Time | Sample match rate |
|---|---|---|---|---|---|
| hg38 → mm39 | 80,818 | 31,772,095 | 158MB | about 31 seconds | 70.3% |
| mm39 → hg38 | 87,089 | 30,805,813 | 156MB | about 31 seconds | 70.4% |

- `liftover` edges **are directional** (`directional`). liftOver chains are deduplicated only on the source side, so they are not used in reverse; the reverse direction comes from the reverse file.
- A chain is followed only when a target species is specified (§2.2; for example `taxon=10090`).
- The cost is 2 (§3). Under the layer rules (§2.1) it is a genome-to-genome edge, so genome → genome → transcript → protein reaches proteins of the other species within one U-turn.
- Results are not marked `approximate` for chains whose self-validation (sample match rate ≥ 0.5) is ok. A cross-species correspondence is a homology-based position and does not guarantee sequence identity.

**Validation** (`service/bench/verify-liftover.ts`, 999 random residues from human MANE Select proteins, converted with `taxon=10090`)

| Item | Result |
|---|---|
| Human codon → mouse genome | 974 / 999 (97.5%) |
| 　That mouse codon encodes the same amino acid | 813 / 974 (83.5%; about 5% if positions were shifted) |
| Human residue → mouse protein | 941 / 999 (94.2%) |
| 　Protein with the same gene name | 855 / 941 (90.9%) |
| Speed (sum of conversion to genome and to protein, `taxon=10090`) | p50 15ms, p95 52ms, p99 100ms → after revising rule 5 and avoiding edge lookups in identity-only expansion, p50 9ms, p95 32ms, p99 59ms |

Most of the 86 with different gene names are differences in ortholog naming (ZNF14 → Zfp709, CYP2F1 → Cyp2f2, REG1A → Reg1) or multigene families (olfactory receptors, bitter taste receptors, histones).

**Example**: converting `uniprot:P07203:49` (human GPX1) with `to=protein&db=uniprot&taxon=10090` gives `uniprot:P11352:47` (mouse Gpx1) via the path identity → CDS → chain → CDS → identity (cost 4). The same holds for `refseq:NP_000572.2:49` → `refseq:NP_032186.2:47`. `refseq:NC_000075.7:106312500..106312550` (mouse chromosome 9) → `refseq:NC_000003.12:complement(join(51986735..51986765,51986769..51986774))` (a deletion in the middle becomes a join).

**Remaining issues**: ortholog correspondences (Ensembl Compara, etc.) can be added as species-crossing edges under the rules in §2.2. UniProt entries whose sequence differs from the reference genome have no identical sequence and cannot be reached. Example: mouse Nras Swiss-Prot `P08556` differs from the GRCm39 translation at 2 residues (168 L/M, 184 S/L). Therefore, from human NRAS (`uniprot:P01111`) only TrEMBL `A0A0G2JDN6`, identical to GRCm39, is reached. Differences that are only substitutions can be handled by a simple mapping between sequences of the same length (T2). Species pairs without chains are handled by T3 (importing GFA, self-computation).

## 10. Mapping between assemblies (GRCh37 ↔ GRCh38, 2026-09-18)

**Stores**: `grch37.sqlite` (GRCh37.p13 assembly report and genomic.fna, 297 sequences), `grch38_names.sqlite` (GRCh38.p14 sequence names), `chain_hg19ToHg38.sqlite` (1,278 chains, 53,950 blocks, 1.1MB), `chain_hg38ToHg19.sqlite` (25,369 chains, 185,410 blocks, 9.7MB). Both chains have a sample match rate of 99.6%.

**Validation** (`service/bench/verify-assembly.ts`, 995 random residues from human MANE Select proteins)

| Item | Result |
|---|---|
| Protein residue → GRCh37 codon (via GRCh38 and chain) | 992 / 995 (99.7%; unreached ones, such as GTPBP6 in the PAR, are not covered by the chain) |
| 　Same amino acid as the GRCh38 codon | 992 / 992 (100%) |
| 　GRCh37 codon → back to the same residue (chain → GRCh38 → CDS) | 989 / 992 (99.7%; the `hg19ToHg38` chain covers less than `hg38ToHg19`) |
| Speed (protein → GRCh37) | p50 6.0ms, p95 18ms, p99 32ms |

**Example**: the BRAF V600E position, GRCh37 `hg19:chr7:140453136` → GRCh38 `refseq:NC_000007.14:140753336` (matches the known mapping). Converting `uniprot:P15056:600` with `assembly=hg19` gives `refseq:NC_000007.13:complement(140453135..140453137)`. Converting `hg19:chr7:complement(140453135..140453137)` to protein gives `uniprot:P15056:600` (GRCh37 → chain → GRCh38 → CDS → identity).

**Mouse GRCm38 (mm10) ↔ GRCm39 (mm39)**: `grcm38.sqlite` (GRCm38.p6 assembly report and genomic.fna, 239 sequences), `grcm39_names.sqlite` (GRCm39 sequence names), `chain_mm10ToMm39.sqlite` (236 chains, 279KB), `chain_mm39ToMm10.sqlite` (910 chains, 512KB). Sample match rates are 98.0% and 97.7%. 1,000 random residues from mouse RefSeq proteins (`verify-assembly.ts mouse.sqlite …`): reach a GRCm38 codon 1,000/1,000, same amino acid 1,000/1,000, back to the same residue 1,000/1,000, p50 1.5ms, p99 11ms. Example: `mm10:chr9:108339451..108339453` (Gpx1) → `uniprot:P11352:47` (chain → CDS → identity); specifying human gives `uniprot:P07203:49`, and `taxon=9606&assembly=hg19` gives `refseq:NC_000003.11:complement(49395565..49395567)` (mm10 → mm39 → hg38 → hg19; the same as converting `uniprot:P07203:49` to hg19).

**Impact**: loading GRCh37 makes human conversions search one step deeper (§2.2 rule 6). The round-trip results with structures (§8) did not change; UniProt → genome went from p50 5.1 → 6.6ms and p99 39 → 53ms.

## 11. Between assemblies without a chain (Marchantia v3.1 ↔ v7.1, 2026-09-18)

**Stores**: `marchantia.sqlite` (GCA_003032435.1 Marchanta_polymorpha_v1 = MpTak v3.1, taxon 3197; rebuilt with the assembly report), `marchantia_v71.sqlite` (GenBank file of GCA_039105155.1 MpTak_v7.1, taxon 1480154; 10 chromosomes, mitochondrion, chloroplast; all 20,412 CDSs ok in self-validation).

| Identical sequences (refget digest) | Count |
|---|---|
| UniProt (UP000244005) entries identical to a v3.1 protein | 19,120 / 19,277 (99.2%) |
| UniProt entries identical to a v7.1 protein | 16,221 / 19,277 (84.1%) |
| v7.1 proteins identical to a v3.1 protein | 18,255 / 20,412 (89.4%) |

- No chain between v3.1 and v7.1 is distributed (marchantia.info has only a gene ID mapping table). The coding regions of the 89.4% of genes above can be converted through identical proteins. Example: `insdc:KZ772678.1:complement(1969369..1969371)` (v3.1) → `insdc:AP031344.1:complement(7591962..7591964)` (chr3 of v7.1; can also be written `MpTak_v7.1:chr3:...`).
- The default assembly is MpTak_v7.1 (the newer one). UniProt → genome returns the v7.1 position, and with `assembly=Marchanta_polymorpha_v1` also returns the v3.1 position.
- To also convert non-coding regions and genes whose sequences changed, whole-genome alignments in both directions were computed with minimap2 and ingested as PAF (spec-ingest §16).

**Alignment validation** (`service/bench/verify-genome-pair.ts`, 2,000 random positions each)

| Item | v3.1 → v7.1 | v7.1 → v3.1 |
|---|---|---|
| A 21-base genomic interval is moved | 1,843 / 2,000 (92.2%) | 1,724 / 2,000 (86.2%; some sequence exists only in v7.1: chrU from Tak-2, repeats resolved by HiFi) |
| 　The sequence at the destination is identical | 99.5% | 97.2% |
| The codon of a residue of a protein identical in both is also moved by the alignment | 99.4% | 99.1% |
| 　Same codon as the answer via proteins | 1,919 / 1,930 (99.4%) | 1,948 / 1,973 (98.7%) |
| Speed | p50 0.1ms, p99 0.9ms | p50 0.1ms, p99 0.3ms |

**Revised genes** (`service/bench/verify-revised-genes.ts`): 1,925 v7.1 proteins are not identical to any protein (their gene model or residues were revised). Converting their residues via v7.1 genome → alignment → v3.1 genome → v3.1 CDS → identical UniProt, 424 of 1,000 random residues (42.4%) reach UniProt, and 397 of those (93.6%) encode the same amino acid on the v3.1 genome as well. The differing ones are revisions that changed the reading frame (the UniProt-side position spans codons, `161c3..162c2`) and revisions of the residue itself (corrections of sequence errors). Unreached residues are in parts not in the v3.1 gene model (such as new exons). Example: `insdc:BFI18695.1:200` → `uniprot:A0A2R6X3H3:192` (the residue number is shifted).

The roughly 1% with differing answers were genes with multiple copies of proteins with the same sequence (tandem duplicates, paralogs). The path through proteins can land on a different copy, and the alignment, which maps by position, is more correct. So among paths of the same cost, **the path with fewer steps** is chosen (in order: cost, number of intermediate sequences without preferred tags, number of steps; for genome → genome, one alignment step beats three steps via proteins). This change did not alter the results of the round trip with structures (§8) or the MANE preference (§2).

## 12. Cis-regulatory elements (fanta.bio CRE, 2026-09-18)

Human (hg38) and mouse (mm10) CREs (promoters, enhancers) were loaded as genome annotations (spec-ingest §17).

- **From a CRE**: converting `fanta:FCHS_301358` (cp1@SOX2) to transcript gives `refseq:NM_003106.4:365..573`. Mouse `fanta:FCMM_194523` (cp2@Gpx1, mm10) gives `refseq:NC_000075.7:108216086..108216674` with `assembly=GRCm39`, and converted to protein gives `uniprot:P11352:<1..55c1` (the end of the CRE overlaps the start of the coding region).
- **From a position**: the annotations (Annotations in the UI) show CREs overlapping that position. Mouse CREs exist only on mm10, but from an mm39 position they are looked up by moving to mm10, moved back to the mm39 position, and shown as "from GRCm38.p6".
- UI examples: "CRE (fanta.bio) → transcript", "mouse CRE on mm10 → mm39 protein".

## 13. Proteins without identical sequences (T2, 2026-09-18)

UniProt entries without identical sequences were aligned with and linked to candidates chosen through ID relations (spec-ingest §18).

- **Path**: `uniprot:P08556:50` → `alignment` (T2) → CDS of `refseq:NP_035067.2` → `refseq:NC_000069.7:102967553..102967555`. The cost is 3 + 1 = 4, higher than a path through identity (0 + 1).
- **Residue differences**: when the path passes a position where the two aligned proteins have different residues, the result gets `differences` (for example `["168 L>M"]`: the residue number as seen from the input side, and input > target residues), and the UI shows "residue differs". `uniprot:P08556:168` → genome lands on a codon that encodes M in GRCm39. Conversely, `refseq:NP_035067.2:168` → `uniprot:P08556:168` gives `168 M>L`.
- **Spreading from the target**: once a sequence matching the target type is reached, T2 alignments are followed in addition to identical sequences. Converting human NRAS (`uniprot:P01111:168`) to mouse UniProt returns, in addition to the identical TrEMBL `A0A0G2JDN6`, Swiss-Prot `P08556` (`168 M>L`).
- **Effect** (`verify-structure.ts`, 3,000 residues from SIFTS intervals): residues that could not reach the genome with only the human stores went from 441 to 420. All 2,580 residues reached landed on codons encoding the UniProt residue, and the round trip with structures was preserved. Speed did not change.

## 14. What differs between input and result (2026-09-19)

Stores do not hold sequences, so differences are judged from the recorded positions and coordinates (spec-ingest §19).

| Result field | Content | Example |
|---|---|---|
| `differences` | Residues and bases that differ along the path. Substitutions in T2 protein alignments (residue number as seen from the input side, and input > target residues), and bases that differ in alignments between assemblies of the same species (`base sequence:position source>target`) | `168 L>M`, `base insdc:AP031344.1:3696566 G>T` |
| `cautions` | `frame differs`: a protein input aligned on residue boundaries landed in the middle of a codon of the target protein (the gene models have different reading frames). `orthologous position in another species`: the path went through a genome alignment (chain) between different species. The position is homologous, but the residue often differs | |

The UI shows `differences` as "residue differs" or "base differs", and `cautions` as a caution marker.

**Examples**: Marchantia v7.1 `insdc:BFI09307.1:713` (E) → UniProt (from v3.1) `A0A2R6XDE2:713` (D) has `base insdc:AP031344.1:3696566 G>T` (likely a correction of a sequence error in v3.1). `insdc:BFI18080.1:193` → `uniprot:A0A2R6VZB3:161c3..162c2` has `frame differs`. `uniprot:P07203:49` → mouse `uniprot:P11352:47` has `orthologous position in another species` (both are selenocysteine, but this is not checked). From a GRCh37 position to a protein, `base …` is added only when the path passes a base that differs between GRCh37 and GRCh38 (not for BRAF V600).
