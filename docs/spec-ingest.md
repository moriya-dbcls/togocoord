# TogoCoord adapter specification (v0.3: GBFF / GFF3 / FASTA / SIFTS)

English | [日本語](spec-ingest.ja.md)

2026-09-18. Rules for the part of §6 of [design.md](design.md) implemented in phase 2. The implementation is in `ingest/`.

---

## 1. Output

An adapter outputs the following three kinds of records (`ingest/src/model.ts`).

| Kind | Content |
|---|---|
| `SequenceRecord` | Reference key, molecule type, unit, length, topology, taxon, refget digest (when the sequence is available) |
| `Edge` | Block list `from → to` (the unit of the core), Location ID for display, attributes, provenance, **self-validation result** |
| `Annotation` | A feature placed by Location ID (its type, and its qualifiers or attributes) |

The CLI (`togocoord-ingest`) writes these as JSON Lines (`{"record": "sequence" | "edge" | "annotation" | "warning", …}`) to standard output, and writes a summary of the validation results to standard error.

---

## 2. Reference keys

- Ensembl stable IDs (`ENST…`, `ENSP…`, `ENSMUSP…` for other species, etc.) get `ensembl:`; an accession with a prefix of the form `XX_` gets `refseq:`; anything else gets `insdc:` (e.g. `refseq:NC_012920.1`, `insdc:AAF99721.1`, `ensembl:ENSP00000407375.1`).
- When a GFF3 seqid is not an accession (such as Ensembl's `1`), `seqidToRef` supplies the mapping. In the CLI, passing an NCBI assembly report to `--seqid-map` maps sequence names (`1`, `MT`), GenBank accessions (`CM000663.2`, `KI270728.1`) and UCSC names (`chr1`) all to RefSeq accessions. When there is no mapping, a warning is issued and the features on that seqid are skipped.
- Ensembl GFF3 writes the version in a separate attribute (`transcript_id=ENST00000419783;version=3`, `protein_id=…;version=1` on CDS lines), so they are joined into `ensembl:ENST00000419783.3`. For IDs without a version, the version is filled in from the IDs of the sequences in `--fasta` (`VersionResolver`).

---

## 3. GenBank / GenPept (`ingestGenBank`)

| feature | Processing |
|---|---|
| `source` | Put taxon and organism into the SequenceRecord |
| `CDS` (nucleotide record, no `/pseudo`) | Edge from the `/protein_id` protein → the record (`cdsMapping`). Uses `/codon_start` and `/transl_table`. aaLength is the length of `/translation` |
| `CDS` (GenPept) | Taking the `/coded_by` location as the CDS, edge from the record itself (protein) → the nucleotide sequence. aaLength is the record length |
| RNA types (mRNA, ncRNA, etc.) with `/transcript_id` | Edge from the transcript → the record (`mappingFromLocation`) |
| All (except `source`) | Annotation (excluding `/translation`) |

- Features whose location cannot be parsed are skipped with a warning. Examples: `bond()`, and `n^1`, which v0.1 does not handle.

---

## 4. GFF3 (`ingestGff3`)

- **Grouping lines**: Lines with the same `ID` form one feature.
- **Ordering intervals**: Lines are sorted by start position. For the `-` strand, that order is reversed.
- **Circular sequences** (`Is_circular=true` on `region`): A line that extends past the sequence length (e.g. the D-loop `16024..17145`, sequence length 16569 nt) is wrapped at the origin and split into two intervals. On a linear sequence, a line that extends past the sequence length is skipped with a warning.
- **Partial features**: `start_range=.,N` becomes `<` (the numerically lower side), and `end_range=N,.` becomes `>` (the higher side).
- **CDS**:
  - Becomes an edge from the `protein_id` protein → the seqid sequence. Without `protein_id`, pseudogenes (`pseudo=true`) are skipped without a warning. Others (such as immunoglobulin gene segments) get a warning.
  - codon_start is `phase + 1` of the line at the 5' end.
  - GFF3 has no translated sequence, so aaLength is estimated by the rules in §5.
- **Transcripts**: Features of an RNA type, or features with `transcript_id` (such as Ensembl gene segments). However, types that represent parts of a transcript, such as exon, CDS and UTR, are not treated as transcripts even when they have `transcript_id` in NCBI data.
- **RNA types** (mRNA, ncRNA, lnc_RNA, etc.): **The location is assembled from the exon lines** linked by `Parent`, because in NCBI GFF3 the mRNA itself is only one line covering the whole gene range. Exons come after their parent, so transcripts are held pending. They are finalized when they have not been updated after a set number (1000) of features, when the input moves to another sequence, or at the end. Without exons, the transcript's own line is used. If there is a `transcript_id`, an edge from the transcript → the sequence is created.
- **Alignments** (lines with `Target`, such as `cDNA_match` and `match`): Become an alignment edge from the `Target` sequence → the seqid sequence.
  - Self-validation: Count identical bases per block; **ok if the identity counted by the adapter itself is at least 50%** (if the Gap is misread and coordinates shift, it drops to around 25%). When the count differs from NCBI's `num_mismatch`, `(counted N)` is recorded in detail (e.g. NM_001291281.3 counts 3, with `num_mismatch=4`). v0.2 judged by agreement with the reported value, so in human there were 23 cases marked as mismatching although the coordinates were correct.
  - Meaning of each Gap (CIGAR) operation:

    | Operation | Meaning |
    |---|---|
    | `M` | Aligned (creates a block) |
    | `I` | Only the Target side advances |
    | `D` | Only the genome side advances |

  - **Gap operations are written in the order of walking the genome in the direction of the line's strand** (ascending genome order for + strand lines, descending for − strand lines). The same rule applies regardless of the Target strand. This cannot be read from the GFF3 specification, so it was confirmed with real NCBI GRCh38 data.
    | Line strand / Target strand | Data used for confirmation | Matches with this rule | Matches with the reverse order |
    |---|---|---|---|
    | − / + | NM_012234.7 (`Gap=M536 I1 M1191 I3 M2314`) | 4036/4041 (5 mismatches = `num_mismatch=5`) | 1977/4041 |
    | + / − | NG_162589.1 (`M217 D5 M53`), NG_043318.1 (`M216 I3 M85`, 2 places) | 270/270, 301/301 | 150/270, 204/301 |
    | + / − (no Gap) | NG_044083.1 | 304/304 | — |
  - Not handled in v0.2: `F`/`R` (frameshifts). A warning is issued.
- **Features spanning both strands** (trans-splicing, e.g. chloroplast and mitochondrial rps12): They cannot be ordered by position, so the line order in the file is used as is, keeping each line's strand. That NCBI GFF3 writes them in transcript order was confirmed by the translation of Arabidopsis NP_051038.1 (`join(complement(69611..69724),139856..140087,140625..140650)`) matching the published sequence.

---

## 5. Estimating protein length (when there is no translated sequence)

Let `coding = CDS length − (codon_start − 1)`.

| Condition | aaLength |
|---|---|
| 3' side partial (`>`) | `floor(coding / 3)` |
| 3' side complete, and `coding` is a multiple of 3 | `coding / 3 − 1` (the last codon is the stop codon) |
| 3' side complete, with a remainder | `floor(coding / 3)` (the remainder is taken as an incomplete stop codon completed by polyadenylation) |

The values estimated by this rule matched the length of `/translation` for all 63 CDSs in the real data used for testing (SARS-CoV-2, the human mitochondrial genome, adenovirus).

---

## 6. Self-validation at ingest

| Kind | Method | Result |
|---|---|---|
| CDS | Extract the CDS sequence, translate it and compare | Compared against: `/translation` (GBFF), the protein record's sequence (GenPept), the published protein sequence given with `--fasta` (GFF3). If none is available, confirm that there is no internal stop codon. For a mismatch in a CDS with an NCBI `/exception` (e.g. "annotated by transcript or proteomic data"), `expected: /exception=…` is recorded in detail |
| Transcript model | Compare the sequence assembled from the genome with the transcript's own sequence (`--fasta rna.fna`) | ok if the lengths are equal and base substitutions are at most 5% (coordinates are preserved). If the transcript is longer and the overhang is at least 90% A, it is taken as a poly(A) tail and the comparison stops before it. Any other length difference is treated as having insertions/deletions: mismatch |
| alignment | Count identical bases per block | If `num_mismatch` is present, judge whether it agrees |

CDS validation takes the following special cases into account (confirmed that all CDSs in Ensembl GRCh38 release 116, about 370,000, are ok).

- **CDS protein length**: If the published protein sequence is available, its length is used (the estimate drops the last residue of an incomplete CDS without a stop codon).
- **Leading `X`**: If the phase is not 0 and the published sequence starts with `X`, this is taken as the Ensembl convention (a missing first codon counts as one residue), and the mapping is built with `leadingPartialCodon` (spec-core §4.2). Without this, over 7,700 cases in human Ensembl were off by one residue.
- **Initiation codon**: In a CDS complete on the 5' side, if the first codon is a start codon in the translation table, or the first residue of the published sequence is M, it is read as M (non-AUG starts such as GTG and ACG; 82 cases in Ensembl).
- **Recoded stop codons**: Positions where the published sequence has U (selenocysteine) or O (pyrrolysine) and the genome translation gives a stop codon are treated as matches and recorded (Ensembl GFF3 has no `transl_except`).
- **Inferring the translation table**: When no translation table is given explicitly (such as mitochondria in Ensembl GFF3) and a mismatch results, a translation table that matches the published sequence exactly is searched for and adopted, and recorded in the `translTable` attribute and in detail. In GBFF, following the INSDC convention, a missing `/transl_table` means translation table 1, and no inference is made.

- Start codons: In a CDS complete on the 5' side, a codon that is a start codon in the NCBI translation table is read as M (e.g. ATT of ND2).
- `/transl_except`: The position is **converted to a residue number by the core's inverse mapping**, and the amino acid is replaced (e.g. the Sec of GPX1 is residue 49). TERM lies outside the protein, so it is ignored. Amino acid names are case-insensitive (GBFF uses `OTHER`, NCBI GFF3 uses `Other`).
- Translation tables: All 27 tables generated from NCBI's `gc.prt` are used (`ingest/src/genetic-codes.ts`).
- If the sequence is not available, the result is `skipped`.
- **Basis of validation (`basis`)**: `full` when checked against the whole sequence, `partial` when only the absence of stop codons was confirmed. In path search, only edges that are ok with `full` are treated as verified even when they carry an exception mark (spec-service §3).
- For a mismatch in a feature with `/exception`, whether CDS or transcript, `expected: /exception=…` is recorded in detail.

---

## 7. Real data used for testing (`ingest/test/fixtures/`, retrieved from NCBI on 2026-09-18)

| Data | What was checked |
|---|---|
| NC_045512.2 (SARS-CoV-2, GBFF and GFF3) | slippage, 12 CDSs |
| NC_012920.1 (human mitochondrial genome, GBFF and GFF3) | translation table 2, incomplete stop codons, alternative start codons, D-loop spanning the origin, 13 CDSs |
| NC_001405.1 (human adenovirus C, GBFF and GFF3) | spliced CDSs on the minus strand, 38 CDSs |
| NC_002127.1 (circular plasmid) | translation table 11 |
| NM_000581.4 (GPX1) | selenocysteine |
| NM_014739.3 and NP_055554.1 (GenPept) | `/coded_by` matches the CDS of the mRNA record |
| 12 `cDNA_match` lines on NC_000003.12 (GRCh38) and related sequences | alignment edges including I (insertion) on the minus strand and D (deletion) on the plus strand |

**Agreement between GFF3 and GBFF**: For 63 CDSs in 3 organisms, Location ID, block list, codon_start and aaLength all matched.

**End-to-end validation**: For every residue of every protein in all fixtures (over 20,000 residues), "convert to a genome position with the core → extract the sequence → translate" matched `/translation`.

---

## 8. Ingesting large data (v0.2)

The design is in [scaling.md](scaling.md).

- **Streaming**: `ingestGff3File` and `ingestGenBankFile` read line by line (`.gz` is decompressed while reading). For GFF3, lines with the same ID are grouped within a window of up to 1000 features (in NCBI GFF3 they are always adjacent). If the same ID appears again after being finalized, a warning is issued. When streaming, the `##FASTA` section is not read (pass it with `--fasta`).
- **Output target (Sink)**: `MemorySink` (for tests and small data), `JsonlSink`, `SqliteSink`. The adapter logic is the same regardless of the output target.
- **FASTA**: Files over 64MB, or files that have a `.fai`, are accessed randomly using the `.fai` (same format as samtools; created automatically if missing).
- **Store**: SQLite (`node:sqlite`). The schema is in scaling.md §4 (v0.2 added a `basis` column to edge and set the schema version to 2). The R*Tree is built in position order after all rows are written.
- **exon**: Not stored as annotation by default (stored with `--all-annotations`).

---

## 9. FASTA adapter (sequence identity)

- `ingestFastaFile`: Outputs each sequence as a SequenceRecord with length, refget digest and MD5 (the sequence itself is not stored).
- Header interpretation: `sp|P07203|…` and `tr|…` → `uniprot:` (isoform suffixes are kept), `ENSP…` → `ensembl:`, `101m_A mol:protein` → `pdb:101M.A` (wwPDB `pdb_seqres.txt`; `mol:na` is skipped); anything else is interpreted as an accession.
- If the same sequence was registered earlier (without a digest), the SQLite store fills only the empty fields from the later record.
- CLI: Inputs `.fa`, `.faa` and `.fasta` (also `.gz`) are ingested with this adapter.

## 10. SIFTS adapter

- Input: `uniprot_segments_observed.tsv(.gz)` (EBI). One line = a UniProt interval and a SEQRES interval of a PDB chain.
- edge: alignment `uniprot:<SP_PRIMARY>` → `pdb:<PDB>.<CHAIN>`. Coordinates are SEQRES numbers (`RES_BEG..RES_END`, equivalent to label_seq_id); author numbering (`PDB_BEG-PDB_END`) is kept in the `authorNumbering` attribute. CHAIN is the author chain ID (auth_asym_id).
- Consecutive lines for the same (UniProt, chain) are merged into one edge. Lines whose lengths differ between the UniProt side and the SEQRES side cannot be made one-to-one, so they are skipped (867 lines in human, 0.06%).
- Self-validation: Count the matching residues between the UniProt and `pdb_seqres.txt` sequences. Below 50% is mismatch (engineered mutations are at most a few residues, so only coordinate shifts are detected). In human, 137 of 240,000 edges are mismatch (all short intervals).
- CLI: `.tsv(.gz)` files whose name contains `sifts` or `uniprot_segments` are ingested with this adapter. With `--sifts-known-only`, only lines for accessions present in the UniProt sequences given with `--fasta` are ingested.

## 11. Sequence tags and MANE (v0.4)

- Added `tags` (a list of strings, e.g. `MANE Select`) and `gene` (gene name) to SequenceRecord. The store schema is version 4. Within the same store, only empty fields are filled from later records. Across stores, tags are merged as a union (spec-service §1).
- MANE adapter (`ingestManeSummary`): Reads NCBI's `MANE.GRCh38.*.summary.txt(.gz)` and outputs SequenceRecords for each gene's RefSeq NM and NP and Ensembl ENST and ENSP, with the tag `MANE Select` or `MANE Plus Clinical` and the gene name. The length is 0 (unknown), and the value from other stores is used. The CLI ingests files named `MANE…summary.txt` with this adapter.
- Ingesting the MANE RNA sequences (`refseq_rna.fna`, `ensembl_rna.fna`) with the FASTA adapter makes it possible to move between the MANE NM and ENST as identical sequences via digest matching (e.g. NM_000581.4 and ENST00000419783.3 match exactly over 899 bases). The FASTA adapter sets the molecule type of transcript accessions (NM, NR, XM, XR, ENST) to RNA.

## 12. Store metadata and content summary (v0.4)

- CLI options: `--label` (display name), `--taxon`, `--organism`, `--assembly`, `--assembly-report FILE` (reads the scientific name, taxon, assembly name and assembly accession from the header of an NCBI assembly report; `--seqid-map` reads the same information). The values are recorded in the store's `meta` table.
- The organism is also picked up from the data: from GBFF `source`, NCBI GFF3 `region` (`Dbxref=taxon:N`) and UniProt FASTA headers (`OS=`, `OX=`). Proteins created from a record inherit its organism.
- When the store is closed, a content summary (`summary`) is recorded in the `meta` table. It contains the number of sequences per molecule type, the number of edges per kind, the number of blocks, the number of annotations, the number of sequences per organism, and examples to try (a residue of a verified CDS protein, that whole protein, and an alignment position). For older stores without a summary, the service computes it on the fly.

## 13. Validation when adding organisms (2026-09-18)

| Organism | Data | CDS self-validation |
|---|---|---|
| Mouse | RefSeq GRCm39 (GFF3, genome, RNA and protein sequences), RefSeq RNA, UniProt UP000000589 | 0 unexplained mismatches (37 with exception) |
| Arabidopsis | RefSeq TAIR10.1, RefSeq RNA, UniProt UP000006548 | 0 unexplained mismatches (64 with exception, such as RNA editing in chloroplasts and mitochondria) |
| Marchantia | **INSDC only** (GenBank files of GCA_003032435.1), UniProt UP000244005 | all 24,674 ok |

- In Marchantia, the genome can be reached through a path of exact matches: UniProt → identical-sequence INSDC protein → scaffold (e.g. `uniprot:A0A2R6VYA2:50` → `insdc:KZ773422.1:208..210`).
- Arabidopsis RefSeq RNA (GenBank) had 4 records whose `/translation` does not fit the CDS (e.g. NM_117687.2: the CDS is 619 bases, but the translation is 844 residues and contains `J`). The same translation appears in all 4, which seems to be an error in the source data. The adapter issues a warning and does not ingest them.
- For large gzip FASTA files (e.g. mouse RNA, which exceeds about 540 million characters when decompressed), the CLI creates a decompressed file next to it and accesses it randomly with `.fai`.

## 14. UCSC chain adapter and storage by chunk (v0.5)

- Input: UCSC chain format (`*.over.chain(.gz)`). Files whose names end in `.chain` or `.chain.gz` are ingested with this adapter. Chain sequence names (such as `chr1`) are UCSC names, so in the CLI, NCBI assembly reports are given to `--from-report` (source assembly) and `--to-report` (target) to translate them into RefSeq accessions (using the `UCSC-style-name` column). Chains on sequences without a mapping (alt, unplaced scaffolds, etc.) are skipped with a warning.
- One chain becomes one edge with `kind: "liftover"` and `directional: true`. Attributes are the chain id, score and target orientation (`qStrand`). For chains whose target is `-`, coordinates are converted to the forward direction (`qSize - q - size`) and the blocks are made `rev`.
- **Self-validation**: When both genome sequences are given with `--fasta`, up to 20 evenly spaced blocks (up to 200 bases each) are extracted per chain and compared. ok if identity is at least 0.5 (`basis: "partial"`). Aligned homologous sequences give 0.6–0.9; shifted positions give around 0.25. For hg38 ↔ mm39 it was about 70% overall.
- **Storage (schema 5)**: Blocks of directional edges go into the `chunk` table in groups of 256 blocks, instead of the `block` table. A chunk holds both sequences, the source-side range (`lo`, `hi`), the position of the first block, the number of blocks, and a byte string encoding each block as "difference from the previous block (source side, target side) and length" in zigzag LEB128. The source-side range is looked up in an R*Tree (`chunk_src`), and only the matching chunks are decoded. About 5 bytes per block: the 31.77 million blocks of hg38 → mm39 take 158MB (about 7GB was expected with the `block` table). Readers accept schemas 4 and 5.
- For the examples in the summary (§12), 30 bases are chosen from a chain's range.

## 15. Assembly sequence names, and assemblies without annotation (v0.5)

- With `--assembly-report FILE`, the store's `meta` also records the mapping from the assembly's sequence names (Sequence-Name, UCSC name, GenBank accession) to RefSeq accessions (`aliases`, JSON), and the UCSC database name (`ucsc`; for GRC assemblies it is not in the assembly report, so it is determined from the built-in table `UCSC_DATABASES`: GRCh38 → hg38, GRCh37 → hg19, GRCm39 → mm39, GRCm38 → mm10). The service uses these to translate input such as `hg19:chr7:140453136` (spec-service §2.2).
- An NCBI assembly report (`*_assembly_report.txt`) itself can also be an input. Each sequence becomes a SequenceRecord with RefSeq accession, length, organism, and DNA (circular for mitochondria and chloroplasts). An assembly without annotation (e.g. GRCh37) is put into one store together with the genome FASTA (for digests).
- A chain store also holds the sequence records (organism and length) of both assemblies from `--from-report` and `--to-report`.
- In assemblies without RefSeq numbers (INSDC only), sequences are represented by GenBank accessions (`insdc:AP031342.1`). The values of the sequence-name mapping are now namespaced sequence keys (`refseq:NC_000007.13`, `insdc:AP031342.1`; un-namespaced values in older stores are taken by the service as `refseq:`).
- The release date of the assembly report (`# Date`) is recorded in `meta.released`. It is used to choose the default assembly (spec-service §2.2).
- `--species-taxon N`: States explicitly the species to which the store's taxon belongs (when the rule that groups subspecies and strain taxa into species does not determine it).

**Human GRCh37 ↔ GRCh38 (2026-09-18)**: `grch37.sqlite` (GRCh37.p13 assembly report and genomic.fna), `chain_hg19ToHg38.sqlite`, `chain_hg38ToHg19.sqlite` (UCSC `hg19ToHg38.over.chain.gz`, `hg38ToHg19.over.chain.gz`). In addition, since the previously built `human.sqlite` has no record of sequence names, a store with only GRCh38 sequence names, `grch38_names.sqlite`, was built from the assembly report (705 sequences).

**Mouse GRCm38 ↔ GRCm39 (2026-09-18)**: `grcm38.sqlite` (GRCm38.p6 assembly report and genomic.fna), `grcm39_names.sqlite` (GRCm39 sequence names), `chain_mm10ToMm39.sqlite`, `chain_mm39ToMm10.sqlite` (UCSC `mm10ToMm39.over.chain.gz`, `mm39ToMm10.over.chain.gz`).

## 16. PAF adapter (whole-genome alignment, T3, v0.5)

- Input: PAF (`*.paf(.gz)`). Output of genome-to-genome alignment by minimap2 or similar; CIGAR (`cg:Z`; minimap2 `-c`) is required. The query is the source and the target is the destination (`minimap2 -c target.fa source.fa`). As with chain, `--from-report` (query assembly) and `--to-report` (target) are used to translate sequence names, and the store also holds the sequence records of both assemblies.
- From the CIGAR, a block list is built by merging adjacent matches (`=`, `X` and `M` advance both, `I` advances the query, `D` and `N` advance the target). On the `-` strand, the target is walked forward and the query backward.
- **Reducing to one-to-one**: PAF reports multiple matches for repeats and paralogs. As with UCSC liftOver chains, only one match is kept for each base of the source (query). Secondary alignments (`tp:A:S`) and alignments shorter than 1kb on the query are removed, and alignments are taken in descending order of score (`AS:i`, or the number of matching bases if absent). Query ranges already covered are trimmed from later alignments. Overlaps on the target side are allowed. One alignment becomes one directional `liftover` edge (with the same compressed storage as chain).
- Self-validation is the same as for chain (ok if the identity of the sampled blocks is at least 0.5).
- **Alignments not ingested** (same for chain; 2026-09-18): (1) Those whose source sequence is shared by both assemblies (the same accession: the nuclear chromosomes of TAIR10 and TAIR10.1, chrM and many unplaced scaffolds of GRCh37 and GRCh38). Their mapping is identity, and the alignments only point to other places such as repeats. (2) Those between a nuclear sequence and an organelle genome (mitochondrion, chloroplast). These are matches to organelle-derived sequences inserted into the nucleus (NUMT, NUPT), not the same position. The molecule type is determined from Assigned-Molecule-Location/Type in the assembly report. Example: aligning the TAIR10 mitochondrion to TAIR10.1 with minimap2 matched well a large mitochondrion-derived insertion on chromosome 2 of Col-0, and reducing to one-to-one selected that. In the human → mouse chain, 392 chains fall into this category.
- **Records for reproducibility** (all stores): MD5 of the input file (`meta.inputs_md5`), MD5 of the `--fasta` sequence files (`meta.sequences_md5`; the genomes the alignment was based on), the TogoCoord commit (`meta.togocoord`; `+local changes` if there are uncommitted changes), and `--method TEXT` (how the input was made: aligner version and arguments, and how it was filtered). Shown in Loaded data.

**Marchantia v3.1 ↔ v7.1 (2026-09-18)**

| | v3.1 → v7.1 | v7.1 → v3.1 |
|---|---|---|
| minimap2 2.31-r1302 `-c --eqx -x asm5 -t 8` | 27 s, 4.3GB | 56 s, 5.3GB |
| PAF lines | 6,372 | 7,150 |
| Alignments kept (one-to-one) | 4,343 (15,175 blocks) | 5,858 (28,766 blocks) |
| Query bases removed as covered by better alignments | 55,948 | 1,781,960 |
| Sampled identity | 99.53% | 98.51% |
| Store | `mp_v31_to_v71.sqlite` 2.7MB | `mp_v71_to_v31.sqlite` 3.4MB |

## 17. BED adapter, and annotations looked up by ID (fanta.bio CREs, v0.5)

- Input: BED (`*.bed(.gz)`). Standard columns 3–12, followed by custom columns. Genomic regions are stored as annotations. BED12 blocks become join, and the `-` strand becomes complement. Chromosome names (`chr1`, `1`, etc.) are translated using the sequence names in `--assembly-report`. UCSC names of patches and unplaced scaffolds not in the assembly report (`chr11_GL456060_alt`, `chr1_KI270706v1_random`) are looked up by the GenBank accession within the name (`lookupSeqid`; also used for chain and PAF sequence names).
- `--bed-type TYPE` (annotation type; default `region`), `--bed-columns NAME,...` (names of the columns after the standard columns; `attributes` splits `key:value|key:value` into attributes). The BED name column becomes the attribute `ID`.
- `--id-namespace NS`: Allows annotations to be entered as `NS:ID` (spec-service §6). Creates an `annotation_id` table (ID → annotation) in the store. `--link URL`: The URL with `{id}` replaced by the ID becomes the annotation's link.

**fanta.bio CRE v1.2.1 (2026-09-18)**: BED9+2 from `https://data.fanta.bio/cre/v1.2.1/` (column 10 is the CRE name, column 11 is `directionality:…|class:PLA/ELA`). The prefix is bioregistry's `fanta` (pattern `^FC(HS|MM)_\d+$`, `https://fanta.bio/cre/$1`).

| Store | Input | Regions | Regions not ingested | Size | Time |
|---|---|---|---|---|---|
| `fanta_human_hg38.sqlite` | `human-CREv1.2.1.hg38.cre-peaks.bed.gz`, `--assembly-report` GRCh38.p14 | 513,895 | 0 | 196MB | 4 s |
| `fanta_mouse_mm10.sqlite` | `mouse-CREv1.2.1.mm10.cre-peaks.bed.gz`, `--assembly-report` GRCm38.p6 | 307,621 | 206 (mm10 patches whose version was bumped in GRCm38.p6; coordinates are not guaranteed to be the same, so they are not ingested) | 115MB | 3 s |

```
togocoord-ingest --db fanta_mouse_mm10.sqlite --assembly-report GCF_000001635.26_GRCm38.p6_assembly_report.txt \
  --bed-type CRE --bed-columns Name,attributes --id-namespace fanta --link 'https://fanta.bio/cre/{id}' \
  mouse-CREv1.2.1.mm10.cre-peaks.bed.gz
```

The activity tables (TPM) and the JSONL annotations (associated genes, TF binding, etc., 446MB) are not ingested. They are not needed for coordinate mapping and can be referred to through the links to fanta.bio.

## 18. Alignment of proteins without an identical sequence (T2, v0.5)

A UniProt entry reaches the genome through a RefSeq, Ensembl or INSDC protein with an identical sequence (refget digest). If even one residue differs, there is no identical sequence; there may be an ID relation (UniProt cross-references), but no coordinate relation saying which residue corresponds to which, so the genome is not reached. Example: the Swiss-Prot entry `P08556` for mouse Nras differs by 2 residues from the GRCm39 translation (`NP_035067.2`) (168 L/M, 184 S/L). This adapter selects candidates through ID relations and aligns the sequences to create coordinate relations (edges).

- Input: UniProt `*_idmapping_selected.tab(.gz)` (by_organism), and the sequences of UniProt and candidate proteins in `--fasta`.
- **Targets**: UniProt entries (including isoforms) not identical to any candidate sequence. Identical ones are reached through the identical-sequence path and are not handled here.
- **Candidates** (tried in order, stopping at the first stage that finds any): (1) RefSeq (column 4), EMBL-CDS (INSDC proteins, column 18) and Ensembl_PRO (column 21) in the entry's row. (2) Proteins listed by other entries of the same gene (GeneID in column 3, Ensembl gene in column 19). (3) Proteins listed by other entries in the same UniRef90 cluster (column 9). UniProt attaches RefSeq to entries whose sequence matches, so the `P08556` row has empty RefSeq and gene, and `NP_035067.2` is in the row of TrEMBL `A0A0G2JDN6` in the same UniRef90_P08556.
- **Alignment**: 6-residue words that occur exactly once in each sequence are used as anchors, and the longest chain of anchors in the same order in both (LIS) is taken. The gaps between anchors are filled with Needleman-Wunsch (match 2, mismatch -1, gap -2). At both ends, end gaps are free. A filled segment is kept only if, between anchors, the lengths are equal (substitutions only) or at least half matches, and at the ends, at least half matches and at least 3 residues match (so that different first or last exons are not aligned by chance). Runs of gapless columns become blocks, and runs with less than half matching are discarded. Segments over 4 million cells are not aligned.
- Among the candidates, the one with the most matching residues is taken, and one `alignment` edge from the entry → the candidate is created (the provenance adapter is `protein-alignment`). Attributes: identity (`identity`), fraction of the entry aligned (`coverage`), substitutions (`substitutions`; up to 200 of `entry position/candidate position:residue>residue`), and how the candidate was chosen (`candidates`). ok if identity is at least 0.9 and coverage at least 0.5; otherwise mismatch (`approximate` in paths).

| Species | UniProt entries | With identical sequence | No candidate | Aligned (ok) | Of which same gene | Of which same UniRef90 | Not aligned |
|---|---|---|---|---|---|---|---|
| Human (RefSeq and Ensembl) | 169,651 | 162,125 | 947 | 6,556 (6,470) | 362 | 117 | 23 |
| Mouse (RefSeq) | 63,328 | 29,794 | 1,819 | 31,257 (30,009) | 24,983 | 926 | 458 |
| Arabidopsis (RefSeq) | 41,596 | 40,369 | 259 | 965 (949) | 17 | 113 | 3 |

Mouse has few identical sequences because Ensembl is not loaded, and Ensembl-derived TrEMBL entries do not match RefSeq. They are connected by aligning them to the RefSeq of the same gene. Marchantia is out of scope because by_organism has no ID mapping for it (99.2% of UniProt connects through identical sequences).

## 19. Recording mismatching bases (schema 6, 2026-09-19)

Stores do not hold the sequences themselves (because of size, synchronization with the distributors, redistribution terms, and because it goes beyond the role of a coordinate service). Instead, they hold **only the positions that differ**.

- For chain and PAF between assemblies of the same species (the taxon in both assembly reports is the same, or one name is the other's binomial name followed by a rank below species), the bases of all blocks are compared against the actual sequences at ingest, and the positions of mismatching bases (0-based position on the source), the source base, and the target base (converted to the source orientation) are recorded in the `mismatch` table (indexed on `(seq, pos)`). The count goes in the edge attribute `mismatches`. `N` is not counted as a difference.
- For chains between different species (human ↔ mouse), most bases differ, so they are not recorded (shown in the note in spec-service §14).
- T2 protein alignments hold substitutions in the edge attribute `substitutions` (§18).

| Store | Mismatching bases | Size |
|---|---|---|
| `chain_hg19ToHg38` / `chain_hg38ToHg19` | 66,673 / 185,228 | 3.5MB / 16MB |
| `chain_mm10ToMm39` / `chain_mm39ToMm10` | 26,056 / 39,541 | |
| `chain_tair10.1ToTair10` / `tair10_to_tair10.1` | 233 / 165 | |
| `mp_v31_to_v71` / `mp_v71_to_v31` | 42,516 / 336,134 | / 15MB |

Schema 6 readers accept 4–6 (a store without a `mismatch` table simply has not recorded differences).

## 20. Alignments towards the annotated assembly (a star, 2026-09-20)

As a species gains assemblies, aligning every pair takes N(N-1) alignments. Instead, **the assembly with the most annotation (the species' default assembly) is the hub, and every other assembly is aligned to it in both directions** (2(N-1) alignments). Other pairs are converted through the hub (spec-service §2.2: a path crosses at most two assemblies of one species).

- The hub is chosen like the default assembly (annotated, the newest release, then the most annotation). For Marchantia it is MpTak_v7.1 (the autosomes of Tak-1 with chrU from Tak-2, chrV from Tak-1 and the organelles, so it can take assemblies of either sex).
- A direct alignment (a chord) can be added later for a pair that needs one; that costs fewer alignments and shorter paths than a second hub.
- The minimap2 preset follows the sequence divergence (minimap2's own guidance): `asm5` (~0.1%, versions of one strain), `asm10` (~1%, another strain of the species), `asm20` (~5%, a divergent accession or subspecies). Well beyond 5% (another species) none of them fits, and a distributed chain is used instead (UCSC builds those with lastz).

**Marchantia (2026-09-20)**: hub MpTak_v7.1; spokes v3.1, MpTak2_v7.1 (Tak-2), ASM993635v2 (v5.1) and cmMarPoly1.2 (no annotation).

| Spoke | Preset | Alignments kept (to hub / from hub) | Sampled identity |
|---|---|---|---|
| v3.1 (GCA_003032435.1) | asm5 | 4,342 / 5,843 | 99.5% / 98.5% |
| MpTak2_v7.1 (GCA_037833965.1) | asm5 | 1,461 / 1,538 | 98.7% / 98.6% |
| ASM993635v2 (GCA_009936355.2, v5.1) | asm5 | 1,344 / 1,283 | 99.6% / 99.1% |
| cmMarPoly1.2 (GCA_965642975.2) | asm20 | 13,501 / 11,314 | 93.4% / 93.7% |

**Share of a spoke that reaches the hub** (`verify-genome-pair.ts`, 1,500 random windows of the spoke's genome): v5.1 97.9%, MpTak2_v7.1 91.6%, cmMarPoly1.2 63.9%. cmMarPoly1.2 is another accession, 2-3% divergent, and its genome is larger (265 Mb against the hub's 248 Mb). A higher preset lifts more of it (asm5 39.6% → asm10 54.7% → asm20 63.9%). The alignment itself covers 60.3% of it with asm10 and 71.5% with asm20, so the limit is the alignment, not the one-to-one filtering.

**Limit on recording differences**: where the differing bases exceed 1% of the aligned bases (`MAX_MISMATCH_RATE`), the positions are not recorded and only the share is kept (the `mismatchRate` attribute). The record exists to point out where two nearly identical sequences differ; between divergent ones the differences are the rule, not the exception. For cmMarPoly1.2 this took the store from 291 MB to 11.6 MB (what remains is the few regions that differ little).

**Marchantia assemblies deliberately not loaded** (2026-09-20)

| Assembly | Reason |
|---|---|
| MpTak1_v7.1 (GCA_037833805.1) | Its sequences are the standard genome's (the hub) under the same accessions (AP031342-AP031350); it would add no conversion |
| Col-CEN v1.2 (Arabidopsis) | No INSDC accession, so its sequences have no key |
| The other assemblies without annotation (Marpolrud_CA_v1, the two subspecies, ASM1997375v1) | No annotation, and each would need its own alignment; cmMarPoly1.2 stands as the example of that case |

Of MpTak2_v7.1 (Tak-2), 81.8% of the proteins (16,641 / 20,354) are identical to one of the standard genome's and the rest are its own. Its chrU is the standard genome's chrU. It is loaded so that positions written on the Tak-2 assembly can be converted.
