# TogoCoord Core Specification (v0.1.1)

English | [日本語](spec-core.ja.md)

2026-09-18. This document makes §3 (Location ID) and §4 (coordinate model) of [design.md](design.md) precise enough to implement. The core library (`core/`) follows this specification and is verified against the test data (`core/test/corpus/`).

---

## 1. Scope

- In scope: parsing Location IDs, converting them to semantic values, and writing them in canonical form. Building, inverting and composing mappings (block lists). Converting locations through mappings.
- Out of scope: fetching sequences, reading data sources (adapters), path search, and ID resolution (filling in the latest version, determining the canonical UniProt isoform, and so on). These are handled in the service layer. The core receives the information it needs (such as the unit of a sequence) from the caller.

---

## 2. Internal coordinates

- Everything is handled as **0-based half-open** intervals `[start, end)`. Conversion to and from 1-based closed notation happens only at parse and output time.
- The coordinate unit is determined by the sequence type.
  - Nucleotide sequences (`nt`): 1 base is 1 unit.
  - Amino acid sequences (`aa`): **1 residue is 3 units** (codon units). Character k (1..3) within the codon of residue n (1-based) corresponds to unit `3(n−1)+(k−1)`.
- With this convention, aa↔nt conversion is a one-to-one correspondence between codon units and bases, and aa↔aa conversion is a correspondence in steps of 3 units. Therefore, **the codon position is preserved even through aa↔aa mappings**.
- The caller gives the unit of a sequence through the function `unitOf(ref)`. If it is not given, the namespace default is used (§3.4).

---

## 3. Location ID

### 3.1 Overall

```
<namespace>:<accession>:<location>
```

- Everything up to the first `:` is the namespace, up to the next `:` is the accession, and the rest is the location. The accession does not contain `:`.
- The internal key of a sequence is `<namespace>:<accession>` (e.g. `refseq:NC_045512.2`, `pdb:4HHB.A`).
- All whitespace is removed before parsing.
- **Whole-sequence shorthand** (v0.1.3): `parseLocationId(text, ctx, { wholeSequence: true })` reads a `namespace:accession` with the location omitted as `1..L`, using the length obtained from `ctx.lengthOf`. This is accepted only as an input shorthand; the output is always an explicit range (e.g. `uniprot:P07203:1..203`). For sequences handled without a version (such as UniProt), the length can change with updates, so returning the length at that time explicitly preserves reproducibility. If the length is unknown, it is an error.

### 3.2 Location grammar

```ebnf
location   = [ remote ":" ] body ;
body       = "complement(" location ")"
           | "join(" location { "," location } ")"
           | "order(" location { "," location } ")"
           | span ;
span       = range | point | between | oneof ;
range      = [ "<" ] pos ".." [ ">" ] pos ;
point      = pos ;
between    = pos "^" pos ;
oneof      = INT "." INT ;
pos        = INT [ "c" ( "1" | "2" | "3" ) ] ;
INT        = decimal integer >= 1 ;
remote     = [A-Za-z0-9_.-]+ containing at least one letter ;    (* accession in the same namespace as the outer one *)
```

**Conditions that are parse errors**

| Condition | Example |
|---|---|
| 0 or a negative number | `0`, `-5` |
| The start of a range is greater than its end | `12..5` (a range across the origin is written with join) |
| `<` on the end side, or `>` on the start side | `1..<5`, `>1..5` |
| `<` or `>` on a single-point position | `<5` |
| A codon extension other than 1..3 | `5c4` |
| Between residues (`^`) where the two sides are not adjacent | `5^7` |
| `oneof` with a codon extension, `<` or `>` | `5c1.7` |
| Nested `order` | `join(order(1..2),3..4)` |

### 3.3 Conversion to semantic values (segment list)

The parsed syntax tree is converted to **a list of oriented segments**.

```ts
Segment = { ref, start, end, strand: +1 | -1, fuzzyLow?, fuzzyHigh?, uncertain? }
Location = { outer: ref, kind: "join" | "order", segments: Segment[] }
```

- **Traversal order**: segments are ordered as the feature is traversed 5'→3' (N-terminus→C-terminus).
- `complement(X)`: the segment list of X is reversed, and the strand of each segment is flipped. `complement(join(A,B))` becomes `[B⁻, A⁻]`, the same value as `join(complement(B),complement(A))`.
- `join` and `order`: the child segment lists are concatenated. Nested joins are flattened. The kind is `order` only when `order` is at the top level.
- **Fuzzy is held as a numeric end**. `<` attaches to the numerically lower side (`fuzzyLow`), and `>` to the higher side (`fuzzyHigh`). As in INSDC, the distinction is by numeric order, not by biological 5'/3'.
- `remote:` sets the ref of the segments in that part to `<outer namespace>:<remote>`.
- **Unit conversion** (for aa):
  - Position `n` → `[3(n−1), 3n)`
  - Position `nck` → `[3(n−1)+k−1, 3(n−1)+k)`
  - Range start `n[ck]` → `3(n−1) + (k−1, or 0 if omitted)`
  - Range end `n[ck]` → `3(n−1) + (k, or 3 if omitted)` (exclusive end)
- **Between residues `a^b`**: becomes a zero-length segment `[k, k)`, where k is the position of the right unit.
  - nt: `b = a+1` is required, and k = a.
  - aa: with the left unit `L = 3(a−1)+(x−1, or 2 if omitted)` and the right unit `R = 3(b−1)+(y−1, or 0 if omitted)`, `R = L+1` is required, and k = R. Therefore `60^61`, `60c3^61c1` and `60c1^60c2` are all valid.
- `oneof a.b`: becomes `[start of a, end of b)` with `uncertain` set.
- **Semantic errors**:
  - strand −1 on an aa sequence (complement cannot be used on proteins)
  - a codon extension on an nt sequence
  - on a sequence of known length, an interval beyond the sequence range

### 3.4 Namespaces

| prefix | accession pattern (approximate) | Default unit |
|---|---|---|
| `insdc` | `[A-Z]{1,6}\d{5,}(\.\d+)?` | 3 letters + digits (e.g. `AAF99721`) is aa, otherwise nt |
| `refseq` | `[A-Z]{2}_[A-Z0-9]+(\.\d+)?` | `NP_`, `XP_`, `YP_`, `WP_`, `AP_` are aa, otherwise nt |
| `ensembl` | `ENS[A-Z]*[EGTP]\d{11}(\.\d+)?` | `…P` is aa, otherwise nt |
| `uniprot` | UniProt accession + `(-\d+)?` | aa |
| `uniparc` | `UPI[0-9A-F]{10}` | aa |
| `pdb` | `[0-9][A-Za-z0-9]{3}\.<chain>` or `pdb_\d{4}[0-9][A-Za-z0-9]{3}\.<chain>` | aa |
| `refget` | `SQ\.[A-Za-z0-9_-]{32}` | No default (given by the caller) |

- The patterns in this table are a minimal check of the ID's form. Namespaces can be added with a registration function (the `test` namespace for tests is also registered this way).
- As normalization, the prefix is lowercased. A PDB ID in 4-character form is uppercased, and an extended-form ID (`pdb_0000xxxx`) that can be converted to 4-character form is converted to 4-character form. The chain is case-sensitive (auth_asym_id).
- **To be confirmed**: prefix names need to be checked against bioregistry. For example, bioregistry has `uniprot.isoform` as a separate prefix. Currently, isoforms are also accepted under `uniprot`.

### 3.5 Canonical form output

Rules for building a string from a semantic value. **If the semantic values are the same, the output strings are always the same**.

1. If there is one segment, that segment is written (wrapped in `complement(…)` if the strand is −1).
2. If there are two or more segments and **all are strand −1 with the same ref**, the form is `complement(join(…))`. The inner elements are written with the segment list reversed and the strand set to +.
3. Otherwise, it is `join(…)` (`order(…)` if the kind is order), with `complement` added to each element as needed.
4. A segment whose ref differs from the outer one is prefixed with `accession:`.
5. Adjacent intervals are not merged.
6. Writing a segment:
   - nt: `n` if the length is 1 with no fuzzy, otherwise `[<]a..[>]b`. Length 0 is `k^k+1`. If uncertain, `a.b`.
   - aa: the start is `residue[c k]` (omitted if k=1), and the end is `residue[c k]` (omitted if k=3). `n` if it covers exactly one residue, `nck` if only one unit.
   - aa of length 0 (position k): if k is a multiple of 3, `k/3 ^ k/3+1`; otherwise `LcX^RcY` (codon position).
7. **Output mode without codon extensions** (`codon: "never"`): `c` is removed, and the form covers whole residues. A `^` within a codon becomes the position `n` of that residue.

---

## 4. Mapping

### 4.1 Block

```ts
Block = { srcRef, src, tgtRef, tgt, len, rev }
```

- A src-side unit `x ∈ [src, src+len)` maps to the following position.
  - `rev = false`: `tgt + (x − src)`
  - `rev = true`: `tgt + len − 1 − (x − src)`
- Invariants: `len > 0`, `src ≥ 0`, `tgt ≥ 0`. **Overlaps are allowed on both the src side and the tgt side** (for ribosomal slippage, duplicated alignments, and closure under invert).
- A mapping (Mapping) is a set of blocks. A unit has zero or more mapped positions.

### 4.2 Operations

| Operation | Definition |
|---|---|
| `invert(m)` | Swaps src and tgt of each block (`rev` unchanged) |
| `compose(ab, bc)` | For each block a of ab and each block b of bc, if the tgt interval of a and the src interval of b intersect in `[lo, hi)`, one block is created. The src-side start is `a.src + (a.tgt + a.len − hi)` if `a.rev`, otherwise `a.src + (lo − a.tgt)`. The tgt-side start is `b.tgt + (b.src + b.len − hi)` if `b.rev`, otherwise `b.tgt + (lo − b.src)`. The length is `hi − lo`, and `rev` is `a.rev ≠ b.rev` |
| `fromLocation(F, loc)` | Mapping from feature sequence F to the reference sequence. Segments are traversed in order, and for the cumulative position `off`, `{F, off, seg.ref, seg.start, len, seg.strand = −1}` is created |
| `cdsMapping(P, cds, codonStart, aaLength)` | `compose(scale, fromLocation(P#cds, cds))`. scale is `{P, 0, P#cds, codonStart−1, 3·aaLength, false}` |

- **Leading partial codon** (`leadingPartialCodon`, v0.1.2): Ensembl places one residue `X` representing the missing codon at the start of the protein of a CDS whose 5' end is missing (INSDC `/codon_start` has no such residue). When this option is specified, the last `codonStart − 1` units of the first residue correspond to the first base of the CDS. scale is `{P, 3−(codonStart−1), P#cds, 0, 3·aaLength−(3−(codonStart−1)), false}`. When the first residue is converted to the genome, the start gets a truncation mark because part of it is missing (e.g. `<930312..930313`).
- `aaLength` is always given by the caller (the core does not infer it), because the core cannot judge the presence of a stop codon, an incomplete stop codon (TERM in `transl_except`), a CDS that is partial on the 3' side, and so on.
- If `3·aaLength + codonStart − 1` exceeds the CDS length, it is an error.

---

## 5. Conversion (`mapLocation`)

For each segment s of the input Location, in traversal order, the following is done.

1. **Pieces**: for every block that intersects s, the intersecting part is mapped. The strand of the mapped part is `s.strand × (rev ? −1 : +1)`.
2. **Piece order**: follows the traversal order of s (ascending src start if +, descending src end if −). For equal positions, the block registration order is used.
3. **Unmapped parts**: the parts of s not covered by any piece are returned, with the same strand as s, in traversal order.
4. **Merging**: consecutive pieces p, q in traversal order (possibly across input segments) are merged into one if all of the following hold.
   - The tgt ref and strand are the same
   - The tgt is contiguous in traversal order (`q.start = p.end` if +, `q.end = p.start` if −)
   - Within the same segment, src advances in traversal order (no overlap)
   - Across different segments, the input kind is `join` and neither segment is uncertain (`a.b`)

   When a genomic interval spanning an intron is mapped, it becomes one contiguous interval on the protein side. A CDS location written as a join per exon, and a codon split at an exon boundary, also become one interval on the protein side. On the other hand, when the same base maps to two positions due to slippage, the pieces are not merged (merging across segments was added in v0.1.1).
5. **Fuzzy**: for the start in traversal order (the start of the first piece within its segment) and the end (the end of the last piece within its segment) of a merged piece, whether "truncation occurred" is determined.
   - If it coincides with the end of the segment, the input fuzzy is carried over (if +, the start is `fuzzyLow`; if −, the start is `fuzzyHigh`).
   - If it is inside the segment, it is truncated if the unit immediately outside in traversal order is not covered by a piece **to the same target sequence**. Whether it maps to other sequences does not matter. Therefore, the result for one sequence does not depend on the presence of unrelated edges (changed in v0.1.1).
   - The truncation is converted to a numeric end of the tgt (if the tgt strand is +, the start becomes Low; if −, the start becomes High).
6. **Between residues (zero-length segment `[k,k)`)**: units k−1 and k are each mapped. **For each target sequence**, if all of the following hold, it maps to `[max, max)` on that sequence.
   - Within that sequence, each has exactly one mapped position
   - The mapped positions are adjacent (absolute difference of positions is 1)

   If the conditions hold for no sequence, it is unmapped. The check is per sequence because, even with overlapping genes in different reading frames (e.g. ATP8 and ATP6 in the human mitochondrial genome), the position is uniquely determined on each sequence (changed in v0.1.1; in v0.1 the condition was exactly one overall).
7. **Uncertain**: if s is uncertain and the merged result is one interval, uncertain is carried over to that interval. If it becomes multiple intervals, `uncertain: true` is set on the whole result.

**Output**
- Merged pieces are grouped by tgt ref (in order of first appearance). One Location per ref is created (kind is the same as the input).
- **When the tgt is aa**: if all strands are −1, the order is reversed, the strand is set to +, and `orientation: "reverse"` is set. If strands are mixed, all are set to + and `"mixed"` is set. Numeric fuzzy is kept as is.
- The result consists of the list of pieces (for display, before merging), the list of mapped Locations, and the unmapped Location.

---

## 6. Not handled in v0.1

- `n^1` on a circular sequence (between the end and the start). It parses, but conversion is an error.
- Checking sequence length. Done only when the length is given.
- ID resolution (filling in the version, normalizing isoforms, matching PDB chains). This is handled in the service layer.

---

## 7. Test data (`core/test/corpus/`)

| File | Source | What it verifies |
|---|---|---|
| `sars2_orf1ab_slippage.json` | NCBI NC_045512.2 (retrieved 2026-09-18). CDS `join(266..13468,13468..21555)`, YP_009724389.1 (7096 aa) | Slippage (one base belongs to two codons), the stop codon being unmapped, truncation fuzzy |
| `human_mt_nd6_minus.json` | ND6 of NC_012920.1 `complement(14149..14673)`, YP_003024037.1 (174 aa) | Minus strand, reverse output for aa, conversion of fuzzy direction |
| `human_mt_nd1_partial_stop.json` | ND1 of NC_012920.1 `3307..4262`, `transl_except=(pos:4261..4262,aa:TERM)`, 318 aa | Incomplete stop codon (aaLength given explicitly) |
| `human_mt_dloop_circular.json` | D-loop of NC_012920.1 `complement(join(16024..16569,1..576))` | Crossing the origin, canonical form of complement(join), input crossing the origin merging into one interval on the feature |
| `synthetic_plus_split_codon.json` | Synthetic data | On the + strand, a codon split at an exon boundary, merging an interval that spans an intron, `^` |
| `synthetic_minus_split_codon.json` | Synthetic data | On the − strand, a split codon |
| `synthetic_aa_to_aa.json` | Synthetic data (modeled on UniProt→PDB) | Codon position preserved in aa↔aa, unresolved regions, `^` on aa |

**Candidates for checking with real data in phase 2 (adapters)**
- Selenoproteins (with Sec specified by `transl_except`)
- Mismatches between RefSeq transcripts and the genome (`cDNA_match` and Gap in NCBI GFF)
- Partial CDSs with `codon_start` 2 or 3 (INSDC)
- PDB insertion codes (antibody structures with Kabat numbering) and unresolved residues
- CDSs crossing the origin in circular bacterial genomes
- GFA (reverse segments, bubbles)
