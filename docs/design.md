# TogoCoord Design Document (v0.1 draft)

English | [日本語](design.ja.md)

2026-09-18

## 0. Background and Approach

TogoCoord is a service that converts sequence coordinates across the different layers of the life sciences: genome, transcript, protein, 3D structure, and so on. The existing proof-of-concept implementations (`poc/html/`, `poc/sparqlist/`) serve only as references; the service is rebuilt with the following approach.

| Issue (PoC) | Approach (this design) |
|---|---|
| Depends on external APIs such as Ensembl REST, TogoWS, and UniParc, which are rich for human and mouse | Build on primary repositories (GBFF, GFF3+FASTA) and generalize to all species |
| The delta-style correspondence tables (with special handling of negative coordinates and complement) are complex and a breeding ground for bugs | Replace them with the idea "location = mapping" and an algebra of block lists |
| SPARQList instances are chained over HTTP, which is slow and error-prone | Separate a pure core library from the services that use it |
| Cross-species conversion is not implemented | Handle it by chaining with compose and by adding alignment edges |
| Rich human and mouse resources | Allow them to be added later as an Enrichment layer on top of the generalized core |

Location IDs remain in use, with stricter rules (§3).

---

## 1. Terminology

| Term | Definition |
|---|---|
| Sequence | A string of residues that actually exists. Registered in the registry, with a namespace, accession.version, molecule type, length, and digest |
| Location ID | A string representing a position or a set of intervals on a sequence. Based on the INSDC location notation |
| Mapping | A correspondence from the coordinates of one sequence to those of another. Represented as a block list |
| Block | A gap-free one-to-one corresponding interval |
| edge | A mapping between sequences in the registry, with a kind and provenance attached |
| Annotation | Information placed on a sequence by a Location ID (domain, PTM, CRE, variant, etc.) |
| Adapter | A module that reads a data source and outputs sequences, edges, and annotations |

---

## 2. Overall Architecture

```
┌─ Enrichment (optional, species-specific, added later) ─────────────────┐
│ MANE/GENCODE, UniProt isoform, SIFTS, AlphaFold DB,                    │
│ UCSC chain, Ensembl Compara, HPRC pangenome, FANTOM CAGE,              │
│ ChIP-Atlas, TogoVar/VEP, domains, secondary structure …                │
└───────────────┬────────────────────────────────────────────────────────┘
                │ output in the same format (sequences / edges / annotations)
┌───────────────┴─ Core (all species, required) ─────────────────────────┐
│ Adapters: GBFF, GFF3+FASTA (INSDC / RefSeq)                            │
│ Normalized store: Sequence registry / Mapping edges / Annotations      │
└───────────────┬────────────────────────────────────────────────────────┘
┌───────────────┴─ Core library (pure, no I/O) ──────────────────────────┐
│ Location ID parsing and normalization / block operations / semantics   │
└───────────────┬────────────────────────────────────────────────────────┘
┌───────────────┴─ Service ──────────────────────────────────────────────┐
│ Path search, per-taxon profiles / REST API / workspaces                │
└────────────────────────────────────────────────────────────────────────┘
  Clients: Web UI, TogoStanza, SPARQList (only calls the API), CLI
```

- The **core library** is written in TypeScript (recommended), so that the same code runs in the browser, Node, SPARQList, and TogoStanza. It is also distributed as a CLI so that users with large data can use it in their own environment.
- **Dependencies go in one direction.** Core works fully without Enrichment. Enrichment sequences are always tied to a Core accession.version or a refget digest.
- No external APIs are called when running a conversion. External data is stored as a snapshot at ingest time (lazy fetching with caching is also acceptable).

---

## 3. Location ID Specification

### 3.1 Overall Format

```
<namespace>:<accession>[.<version>]:<location>
e.g. refseq:NM_014739.3:join(233..235,5958..6116)
    insdc:NC_000001.11:complement(join(201..300,401..500))
    uniprot:Q9BYF1-1:60
    refget:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2:1..100
```

- `namespace` is required. It is chosen from prefixes registered in bioregistry and similar registries (e.g. `insdc`, `refseq`, `ensembl`, `uniprot`, `pdb`, `refget`).
- **Version handling is decided per namespace** (§3.7). For versioned systems (INSDC, RefSeq, Ensembl), output always includes the version. If the input has no version, it is resolved to the latest version. UniProt and PDB are handled without versions.
- The accession does not contain `:`. Parsing treats everything up to the first `:` as the namespace, and up to the next `:` as the accession.
- **The assembly is not included in the ID**, because an INSDC or RefSeq sequence accession.version is unique by itself. The assembly is kept as metadata in the registry.
- **Coordinates are attached only to sequences that actually exist.** Layer specifiers such as `mrna#`, `-pre_mRNA`, and `-cDNA` are not used. The layer is determined by the molecule type in the registry. A pre-mRNA is represented as a location on the genome.

### 3.2 Location Grammar (a Subset of the INSDC Notation)

```ebnf
location   = complement | join | order | span ;
complement = "complement(" location ")" ;
join       = "join(" item { "," item } ")" ;
order      = "order(" item { "," item } ")" ;
item       = [ remote ":" ] location ;          (* remote = accession.version; namespace is the same as the outer one *)
span       = range | position | between | oneof ;
range      = [ "<" ] pos ".." [ ">" ] pos ;
position   = pos ;
between    = INT "^" INT ;                      (* between two adjacent residues *)
oneof      = INT "." INT ;                      (* any one residue within the range *)
pos        = INT [ codon ] ;
codon      = "c" ( "1" | "2" | "3" ) ;          (* extension: protein references only. §3.4 *)
INT        = integer >= 1 ;
```

- Negative coordinates, 0, and whitespace are not allowed.
- A span crossing the origin of a circular sequence is written as `join(4000..4641652,1..100)`.

### 3.3 Canonical Form

The same set of intervals is normalized to a single string. The canonical form is used when comparing or storing IDs.

1. Remove whitespace. Write `n..n` as `n`.
2. Omit remote references that point to the same reference as the outer one.
3. If all elements of a join are complements with the same reference, combine them into the form `complement(join(…))` (the element order is reversed according to INSDC semantics). Otherwise, complement is applied to each element.
4. Adjacent intervals are **not merged** (to keep exon boundary information). Equivalence is checked by a separate function (comparison as interval sets).
5. When a codon extension covers all of 1 to 3 (e.g. `12c1..12c3`), omit the extension.

### 3.4 Codon Position Extension (Optional)

- When a position derived from nt is exposed as an aa ID, the codon position is lost. This extension is used when it should be kept.
- The syntax is `<residue number>c<1|2|3>`. It can be used only in IDs that reference a protein, and cannot be combined with complement.

```
uniprot:Q9BYF1-1:60c2          2nd base of the codon of residue 60
uniprot:Q9BYF1-1:60c2..63c1    range
uniprot:Q9BYF1-1:60c2..63      if omitted, start is taken as c1 and end as c3
uniprot:Q9BYF1-1:60c1^60c2     between two bases within a codon
```

- In the canonical form, the omittable `c1` (start) and `c3` (end) are omitted.
- In nt units, this corresponds to position `3(N−1)+k` on the CDS. The codon position is uniquely determined regardless of which CDS encodes the protein.
- aa↔aa mappings (UniProt→PDB, between isoforms, between orthologs) carry the codon position over unchanged. Therefore, a round trip genome→UniProt→PDB→genome returns to the original single nucleotide.
- Rejected candidates: `60.2` (conflicts with INSDC's "one of the range"), `60+2` (confusable with HGVS intronic positions), `60:2` (conflicts with the reference separator), `60#2` (conflicts with IRI fragments).
- **Removing the extension part is guaranteed to yield a valid aa ID as is (a range covering the residue).** Processors that do not know the extension still read it correctly at aa resolution.
- In practice, for nt-derived data the nt-side ID (a location on the genome or mRNA) is authoritative, and the aa ID is a derived view (the same relationship as HGVS `c.`/`g.` and `p.`).

### 3.5 IRI (identifiers.org Style)

```
https://<togocoord-domain>/refseq:NM_014739.3:join(233..235,5958..6116)
                           └────────┬───────┘ └───────────┬───────────┘
                      resolvable by identifiers.org   location
```

- Uses the same `<prefix>:<accession>` form and prefix names as identifiers.org, with resolution done on TogoCoord's own domain. Content negotiation returns HTML, JSON, or JSON-LD.
- **Removing the location yields a compact identifier that identifiers.org can resolve as is.**
- In the future, registering a `togocoord` prefix with identifiers.org / bioregistry and redirecting from there will also be considered (whether a registration whose local IDs contain parentheses and colons is accepted needs to be checked).
- Encoding: RFC 3986 allows `( ) , : .` to appear as is in a path. **Only three characters are percent-encoded: `<` (%3C), `>` (%3E), and `^` (%5E).**
- Because IRIs are built from canonical-form IDs, "the IRIs match" and "the interval sets match" mean the same thing.
- In Turtle, IRIs containing parentheses cannot be written as prefixed names, so they are written in the full `<…>` form.

### 3.6 FALDO JSON-LD

Follows the official FALDO vocabulary (details in [spec-service.md](spec-service.md) §7). Strand is expressed by the position type (`faldo:ForwardStrandPosition` / `faldo:ReverseStrandPosition`), and on the reverse strand begin is the numerically larger side. `join` is `faldo:ListOfRegions` and `order` is `faldo:BagOfRegions`, with elements ordered by `rdf:_n`. `<` and `>` are `faldo:FuzzyPosition`, `^` is `faldo:InBetweenPosition`, and `.` is `faldo:InRangePosition`. (`NegativeStrand`, originally written in this section, is not FALDO vocabulary and has been corrected.)

- Location IDs and JSON-LD correspond 1:1.
- The codon extension will be proposed to FALDO as an addition of `faldo:codonPosition` (value 1..3) to `faldo:Position`. If it is not accepted into FALDO, it is defined as `tgc:codonPosition` in our own namespace. Processors that do not know this property read it correctly as an aa position. The current implementation outputs `tgc:codonPosition`.

```turtle
[] a faldo:ExactPosition ; faldo:position 60 ;
   faldo:reference <…/Q9BYF1-1> ; tgc:codonPosition 2 .
```

### 3.7 ID Rules per Namespace

For each namespace, the following four points are decided.

1. **Prefix**: follows identifiers.org / bioregistry prefixes (decided).
2. **Accession syntax**: a regular expression per namespace. The meaning of separator characters is also decided here.
3. **How the sequence is pinned**: which notation pins down a single sequence, and how a missing version is handled.
4. **The sequence that coordinates attach to, and its numbering scheme.**

| Namespace | Accepted input forms | Canonical output form | Coordinates |
|---|---|---|---|
| insdc / refseq / ensembl | `acc.ver`, `acc` (without version, resolved to the latest version) | `acc.ver` | Residue numbers of that sequence |
| uniprot | `P12883`, `P12883-2`, the canonical isoform written with a number such as `P12883-1`, etc. | No version. **The canonical isoform has no number** (`P12883`); others are like `P12883-2` | Residue numbers of the latest sequence |
| pdb | `4HHB.A` (chain). Extended PDB IDs (`pdb_00004hhb.A`) are also accepted | No version. `4HHB.A` | **label_seq_id** of that chain |
| refget (GFA segments, user data) | `SQ.<digest>` | Same | Residue numbers of that sequence |

**UniProt**
- Few users manage versions or retrieval dates, so entries are handled without versions and always resolved to the latest sequence.
- The canonical isoform is the one UniProt designates as canonical (not necessarily `-1`). If the input comes with a number (e.g. `P12883-1`) and it is the canonical isoform, it is treated the same as the unnumbered form. Which form is used internally is up to whatever is easier to implement.
- When a sequence is updated, the position the same ID points to may shift. For users who need reproducibility, the digest of the sequence used for resolution, the UniParc ID (UPI), and the UniProt release are returned as response metadata. They are not included in the ID itself.

**PDB**
- Handled without versions.
- Only in the PDB namespace, **as an exception, `.` is treated as the chain separator** (`4HHB.A`). The unit of reference is the chain.
- The chain ID is the **auth_asym_id** (the chain name shown in papers and viewers, e.g. the H/L chains of an antibody) (decided). RCSB's instance notation (`4HHB.A`) uses label_asym_id, and the two can differ. The correspondence with label_asym_id is kept in the registry.
- Coordinates are label_seq_id (sequential numbers starting at 1). Author numbering (such as `52A` with an insertion code) is not used for coordinates; it is returned as a correspondence table in the response.
- Observed and unresolved residues are kept as annotations on the chain.

**Common**
- The registry allows conversion in both directions between the "human-readable form (`ns:acc[.ver]:loc`)" and the "digest form (`refget:SQ…:loc`)".
- Be liberal in input, return the canonical form in output.

---

## 4. Coordinate Model and Mapping

### 4.1 Central Idea

> **An INSDC location string is itself a mapping from "feature sequence (1..L) → reference sequence".**

For each CDS, mRNA, and exon feature in GBFF, the location as written becomes an edge directly. No dedicated code is needed to build conversion tables. GFA paths can be represented in the same form (§7).

### 4.2 Internal Representation

```ts
// Everything internal is 0-based half-open. Conversion to/from 1-based closed happens only at the I/O boundary
type Block   = { src: number; tgt: number; len: number; rev: boolean };
type Mapping = { from: SeqRef; to: SeqRef; blocks: Block[]; unit: Unit; provenance: Provenance };
// Invariant: blocks are sorted by src ascending. Overlaps are allowed on both sides (to handle ribosomal slippage and duplications in alignments,
//            and because invert swaps src and tgt). When there are overlaps, map returns all matching positions
```

There are only two mapping rules: `y = tgt + (x - src)` when `rev = false`, and `y = tgt + len - 1 - (x - src)` when `rev = true`. No negative coordinates or complement-specific deltas are used.

Example: `complement(join(201..300,401..500))` (length 200)

| src | tgt | rev |
|---|---|---|
| [0,100) | [400,500) | ✓ |
| [100,200) | [200,300) | ✓ |

### 4.3 Operations

| Operation | Description |
|---|---|
| `fromLocation(loc)` / `toLocation(mapping)` | Convert between a Location ID and a block list |
| `map(m, interval)` | Map an interval. The result is split into the mapped part and the part that could not be mapped (unmapped) |
| `invert(m)` | Swap src and tgt. No separate implementation of the reverse conversion is needed |
| `compose(a→b, b→c)` | Build a→c by intersecting intervals. Every chained conversion is expressed with this |
| `scale(protein ↔ CDS)` | Map the i-th aa (0-based) to the CDS nt interval [3i+φ, 3i+φ+3). φ is the phase determined by `/codon_start` |

- The difference between aa and nt is handled only inside `scale`. All other mappings are handled as nt↔nt (or aa↔aa) with a unit ratio of 1.
- In chained conversions, intermediate results are composed without being extracted as IDs, so nt resolution is preserved.

### 4.4 Conversion Semantics

| Situation | Rule |
|---|---|
| Only part of an interval can be mapped | Split the result into the mapped pieces. Truncated ends are marked with `<` (start side) / `>` (end side). Unmapped parts are returned separately |
| nt → aa | Expand to the aa interval corresponding to the codons covering the interval. The codon positions at the ends are returned as metadata, and included in the ID with the codon extension (§3.4) on request |
| aa → nt | Map to all three bases of the codon. A codon spanning an exon boundary becomes a join |
| `^` (between residues) | Treated as a zero-length interval [k,k). Mapped only if the residues on both sides are adjacent in the target; otherwise unmapped |
| `a.b` (one of the range) | Mapped as a range, keeping an "uncertain" flag |
| `order()` | Map each element, preserving order |
| Stop codon | Included in the CDS, but not mapped to the protein |
| Partial CDS (`<`/`>`) | Mapped taking the phase into account. The end residues are treated as uncertain |
| The output consists of multiple intervals | Joined in the order of the source and output in the canonical form of §3.3 |

---

## 5. Data Model

### 5.1 Sequence registry

| Field | Description |
|---|---|
| `namespace`, `accession`, `version` | Components of the ID |
| `digest` | GA4GH refget sha512t24u. Used to determine whether sequences are identical across DBs |
| `moltype` | DNA / RNA / protein |
| `length`, `topology` | Length, and whether linear or circular |
| `taxon`, `assembly` | NCBI Taxonomy ID and the assembly it belongs to |
| `aliases` | List of alternative names, such as UCSC chromosome names (`chr1`) |

### 5.2 Mapping edge

| Field | Description |
|---|---|
| `from`, `to`, `blocks`, `unit` | As in §4.2 |
| `kind` | `annotation` / `identity` / `alignment` / `liftover` / `orthology` / `graph`, etc. Kinds can be added later |
| `provenance` | Data source, file, release, feature or qualifier, which tier (T0–T3), and the kind of origin (distributed / self-computed / user) |
| `validation` | Result of self-validation at ingest time (§6.3) |

### 5.3 Annotation

- Information placed on a sequence by a Location ID (domain, PTM site, CRE, variant, secondary structure, PDB observed residues, etc.) is kept as a layer separate from mappings.
- **Propagating an annotation only requires applying `map` to the annotation's Location ID.**

### 5.4 Storage

Stored per assembly in SQLite, DuckDB, Parquet, or similar. Block lists are held in array-typed columns. Publication as RDF (FALDO) is done separately.

---

## 6. Data Sources and Adapters

### 6.1 Adapter Interface

```ts
interface Adapter {
  meta: { id: string; layer: "core" | "enrichment"; taxa: number[] | "any";
          license: string; updateCycle: string; tier: "T0" | "T1" | "T2" | "T3" };
  ingest(input): AsyncIterable<SequenceRecord | Mapping | Annotation>;
}
```

Adding a human or mouse resource only requires adding one adapter, without touching the core code.

### 6.2 Tiers and Data Sources

Tiers are defined by "what granularity of correspondence is possible". Where the data comes from is treated as a separate axis.

| Tier (capability) | Distributed data | Self-computed | Brought in by users |
|---|---|---|---|
| **T0** Annotation and identical sequences | GBFF / GFF3 | ― | GFF3 + FASTA |
| **T1/T2** Correspondence between sequences | NCBI `cDNA_match`, SIFTS | Pairwise protein alignment, spliced alignment | PAF, etc. |
| **T3** Whole-genome correspondence | UCSC chain, Ensembl Compara, GRC alignments, HAL | minimap2, wfmash | GFA, chain, PAF, MAF |

What each tier enables:

- **T0** (all species)
  - Conversion between genome ↔ mRNA ↔ CDS ↔ protein ↔ exon within the same assembly (including annotated isoforms and alt. ORFs)
  - Correspondence between DBs with identical sequences (GenBank protein, RefSeq, UniProt, Ensembl, AlphaFold DB)
  - Correspondence between GCA/GCF with identical sequences
- **T1** (mainly model organisms)
  - Absorbing mismatches between RefSeq transcripts and the genome
  - Residue-level correspondence between UniProt and PDB
  - Liftover between assemblies, and cross-species conversion including non-coding regions between human and mouse
- **T2** (all species, low cost)
  - Correspondence between DBs whose sequences differ slightly (e.g. mismatches between MANE and UniProt)
  - Correspondence between isoforms
  - **Cross-species conversion limited to coding regions**, via protein alignment of ortholog pairs
  - Placing mRNAs not annotated on the genome onto the genome

  All of these are computed on demand and cached.
- **T3**
  - **Cross-species conversion including non-coding regions**, and conversion between assemblies for non-model organisms
  - Handled first with distributed data and user-provided data. Whether to compute it ourselves will be decided after measuring demand

### 6.3 Self-validation at Ingest Time

- Extract and translate the CDS sequence, and check it against `/translation` or the protein sequence. Mismatches are flagged and switched to alignment edges.
- Special cases that must be handled:
  - `/codon_start`
  - `/transl_except`
  - ribosomal slippage
  - RNA editing
  - partial CDS
  - `/exception`
- Identity edges are created only when the digests match.

### 6.4 Correspondence of Core and Enrichment (Elements of the Poster Figure)

| Figure element | Layer | Representation |
|---|---|---|
| Genome / mRNA / CDS / Protein / Exon | Core | The feature location is used directly as a mapping |
| Protein isoform, alt. Splicing | Core (UniProt isoforms are Enrichment) | A separate sequence node for each transcript and protein |
| alt. TSS | Enrichment (FANTOM CAGE, etc.) | Kept as annotations. When a new transcript model is built, it is added as a sequence node |
| alt. ORF | Core / Enrichment | One mRNA has multiple CDS mappings |
| Structure, Unresolved residues | Enrichment (SIFTS, mmCIF, AlphaFold DB) | The chain (`4HHB.A`) is a node, with label_seq_id as coordinates. Observed residues are kept as annotations on the chain |
| Domain, α helix / β sheet | Enrichment | Kept as annotations on the protein |
| Genome of other organism | Enrichment (chain, Compara) | `liftover` / `orthology` edges |
| Pangenome graph | Enrichment / user-provided | Kept as `graph` edges (§7) |
| CREs, Variant | Enrichment (ChIP-Atlas, TogoVar) | Kept as annotations on the genome and propagated through mappings |

---

## 7. GFA and User-provided Data

### 7.1 Handling GFA

A GFA path can be written directly with INSDC location remote references and complement.

```
P  hapA  s1+,s2-,s3+  *
→ hapA = join(s1:1..L1, complement(s2:1..L2), s3:1..L3)
```

- Segments are registered in the registry as actual sequences (with digests).
- A path (and likewise a W line) becomes a mapping edge "path sequence → segments".
- The correspondence between two paths is obtained with `compose(pathA→segs, invert(pathB→segs))`. No new operation is needed.
- Segments that appear in only one of the paths are unmapped. In the target, their position is indicated with `^`.
- Graphs with overlaps in L lines are trimmed and converted to blunt graphs at ingest time.
- Path names are parsed using PanSN naming (`sample#hap#contig`). They are also tied to known sequences by digest matches, so that annotations from public data can be propagated.
- For large graphs, combinations of all paths are not precomputed. Only the requested pairs are composed on the fly. If rGFA tags (SN/SO/SR) are present, they are used as a shortcut to coordinates on the reference sequence.

### 7.2 Workspaces

- Provided data is isolated per workspace. It is private by default, and its origin kind is recorded as `user`.
- User-specific sequences are referenced by `refget:` digests. The same sequence loaded elsewhere gets the same ID.
- At ingest time, data is validated by digest matches.
- Limits are set on the size of files that can be ingested and on the scale of compose per query. Data beyond these limits should be handled with the distributed CLI.

---

## 8. Path Search

- Search for weighted shortest paths on a graph with sequences as nodes and edges as edges. Weights are set per kind, roughly in the order identity < annotation < alignment (distributed) < alignment (self-computed) < liftover/orthology.
- **Per-taxon profiles** allow the priorities to be overridden by configuration. For example, MANE Select is preferred for human. For species without a configuration, the default policy (prefer annotation from primary repositories) applies.
- Multiple candidate paths can be returned, not just one (e.g. when all corresponding PDB structures are wanted).
- **Queries for which no conversion path was found** (target species, source and target layers) are logged and used as input for decisions on investing in T3.

---

## 9. API (Proposal)

| Method | Path | Description |
|---|---|---|
| GET/POST | `/v1/convert` | `loc` (Location ID), `to` (namespace, sequence, layer), `via`, `profile`, `codon`. POST for batch processing |
| GET | `/v1/location/parse` | Returns the canonical form and JSON in block-list form |
| GET | `/v1/location/faldo` | Returns FALDO JSON-LD |
| GET | `/v1/sequences/{ref}` | Registry information |
| GET | `/v1/sequences/{ref}/neighbors` | List of directly connected edges |
| GET | `/v1/annotations` | Returns annotations overlapping the given interval, propagated to the target coordinates |
| POST | `/v1/workspaces/{id}/datasets` | Users provide their own data |

Example response:

```json
{
  "input": "uniprot:Q9BYF1-1:60",
  "results": [{
    "location": "insdc:NC_000020.11:complement(join(…))",
    "pieces": [{ "source": "60", "target": "…", "codon": {"begin": 1, "end": 3} }],
    "unmapped": [],
    "path": [{ "from": "…", "to": "…", "kind": "identity", "provenance": {…} }],
    "faldo": {…}
  }]
}
```

---

## 10. Testing Strategy

1. **Use a naive implementation as the ground truth (oracle)**: prepare a separate implementation that expands every mapping into per-residue arrays, and compare its results with those of the block operations (property-based tests with random input).
2. **Verify properties of the operations**:
   - `invert(invert(m)) = m`
   - `mapped(x) ⊆ map(invert(m), map(m, x))`. Equality holds if m is injective
   - compose is associative
   - Round-tripping through parsing and output yields the same canonical form
3. **Test data collecting difficult cases**:
   - Reverse strand
   - Codons split by exon boundaries
   - codon_start of 2 or 3
   - ribosomal slippage
   - transl_except
   - Partial CDS
   - Crossing the origin of a circular genome
   - Mismatches between RefSeq and the genome (cDNA_match)
   - PDB insertion codes and unresolved residues
   - Reverse-oriented segments and bubbles in GFA
4. **Self-validation at ingest time** (§6.3) is also used as a regression test on the data.

---

## 11. Roadmap

| Phase | Description |
|---|---|
| 0 | Finalize this specification. Prepare test data. **Done (2026-09-18)**: [spec-core.md](spec-core.md), `core/test/corpus/` |
| 1 | Core library (parser, normalization, block operations, semantics) and oracle tests. **v0.1 done (2026-09-18)**: `core/` (104 tests) |
| 2 | GBFF and GFF3 adapters. To confirm generality, species with different characteristics (bacteria with circular genomes, Arabidopsis, etc.) are included from the start in addition to human and mouse. **v0.1 done (2026-09-18)**: `ingest/`, [spec-ingest.md](spec-ingest.md). Validated with real data from viruses, the human mitochondrial genome, adenovirus, plasmids, and human GRCh38 cDNA_match. Whole-genome-scale GFF3 for human and mouse, Arabidopsis, and Ensembl seqids are not yet validated |
| 3 | REST API, path search, Web UI. T1 Enrichment (SIFTS, cDNA_match, UCSC chain, MANE). **3a and 3b done (2026-09-18)**: `service/` (path search, multiple stores), [spec-service.md](spec-service.md). **3c done**: REST API, FALDO JSON-LD, paths via identical sequences (digests). **3e (SIFTS) done**: ingest of Ensembl, UniProt, and SIFTS, and validation of round trips with structures (spec-service §8). **3d (Web UI) done**: spec-service §6.1. **MANE done**: sequence tags and preference at equal cost (spec-service §2). **UCSC chain done**: human ↔ mouse liftOver (spec-service §9, spec-ingest §14). **Phase 3 done** |
| 4 | T2 (on-demand alignment and caching), annotation propagation (reproducing the poster use cases). **Part of T2 (2026-09-18)**: UniProt entries without an identical sequence are aligned at ingest time with candidates selected via ID mapping (spec-ingest §18, spec-service §13). On-demand computation and annotation propagation across assemblies are in spec-service §6 |
| 5 | T3 (user-provided GFA and chain, workspaces). Self-computation will also be considered depending on demand. **Partly started (2026-09-18)**: PAF adapter (spec-ingest §16) and minimap2 alignment of Marchantia v3.1 ↔ v7.1 (spec-service §11). GFA and workspaces not started |

---

## 12. Open Issues

- [ ] Proposal of `faldo:codonPosition` to the FALDO developers (§3.6)
- [ ] TogoCoord's domain name. Whether a prefix can be registered with identifiers.org (§3.5). The IRI form is decided to be identifiers.org style
- [ ] Where to operate the service (DBCLS / DDBJ), and the scale of computing resources
- [ ] Source of ortholog information (OrthoDB, eggNOG, Ensembl Compara)
- [ ] Handling of logs of queries for which no path was found (privacy)
