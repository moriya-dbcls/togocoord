// GFF3 adapter: features (Core tier T0) and cDNA_match / match alignments with Target + Gap (tier T1).
import {
  cdsMapping,
  formatLocationId,
  mappingFromLocation,
  NamespaceRegistry,
  type Block,
  type CoordContext,
  type Location,
  type Segment,
  type Unit,
} from "@togocoord/core";
import { accessionRef, ChainedSource, DEFAULT_EXCLUDED_ANNOTATIONS, extent, ingestContext, ownString, RNA_TYPES, type AdapterOptions } from "./common.ts";
import { attr, parseGff3, type Gff3Document, type GffFeature, type GffRow } from "./gff3.ts";
import { MemorySink, type Edge, type IngestResult, type Provenance, type SequenceRecord, type Sink } from "./model.ts";
import { checksums, MemorySequenceSource, type SequenceSource } from "./sequence.ts";
import { alignmentIdentity, cdsLength, expectMismatch, inferAaLength, validateCds, validateTranscript } from "./validate.ts";

/** Ingest a whole GFF3 text (including a `##FASTA` section) into memory. For large files use `ingestGff3File`. */
export function ingestGff3(text: string, options: Gff3Options = {}): IngestResult {
  const sink = new MemorySink();
  ingestGff3Document(parseGff3(text), sink, options);
  return sink.result;
}

export function ingestGff3Document(doc: Gff3Document, sink: Sink, options: Gff3Options = {}): void {
  const ingestor = new Gff3Ingestor(sink, options);
  for (const [seqid, length] of doc.regions) ingestor.sequenceRegion(seqid, length);
  for (const [seqid, residues] of doc.sequences) ingestor.residues(seqid, residues);
  for (const f of doc.features) if (f.type === "region" && attr(f, "Is_circular") === "true") ingestor.markCircular(f.seqid);
  for (const f of doc.features) ingestor.feature(f);
  ingestor.finish();
}

export interface Gff3Options extends AdapterOptions {
  /** seqid -> internal key; defaults to INSDC/RefSeq accession detection (see assemblyReportSeqids for `1`, `chr1`). */
  seqidToRef?: (seqid: string) => string | undefined;
  /** Completes version-less keys (e.g. Ensembl protein IDs); default: identity. */
  resolveRef?: (ref: string) => string;
}

/**
 * Sequence key of a transcript or protein named in GFF3 attributes. Ensembl writes the version separately
 * (`transcript_id=ENST00000456328;version=2`); NCBI writes `transcript_id=NM_000581.4`.
 */
function attributeRef(id: string | undefined, version: string | undefined, registry: NamespaceRegistry, resolve: (r: string) => string): string | undefined {
  if (!id) return undefined;
  const versioned = version !== undefined && !/\.\d+$/.test(id) ? `${id}.${version}` : id;
  const ref = accessionRef(versioned, registry);
  return ref && resolve(ref);
}

const TRANSCRIPT_PARTS = new Set([
  "exon", "CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon", "intron",
  "polyA_signal_sequence", "polyA_site", "transcription_start_site", "mature_protein_region_of_CDS", "signal_peptide_region_of_CDS",
]);

/** Identity below which an alignment edge is reported as a mismatch (misaligned coordinates). */
const MIN_ALIGNMENT_IDENTITY = 0.5;

const ALIGNMENT_TYPES = new Set(["cDNA_match", "EST_match", "match", "nucleotide_match", "translated_nucleotide_match", "protein_match"]);

interface Molecule {
  ref: string;
  length?: number;
  circular: boolean;
}

/**
 * Feature-at-a-time GFF3 adapter. Call `sequenceRegion` for `##sequence-region` directives, `feature` for each
 * grouped feature (a molecule's `region` feature must precede its other features for topology to apply), then
 * `finish`.
 */
export class Gff3Ingestor {
  readonly #sink: Sink;
  readonly #registry: NamespaceRegistry;
  readonly #units = new Map<string, Unit>();
  readonly #ctx: CoordContext;
  readonly #toRef: (seqid: string) => string | undefined;
  readonly #residues = new MemorySequenceSource();
  readonly #source: SequenceSource;
  readonly #base: Provenance = { adapter: "gff3" };
  readonly #molecules = new Map<string, Molecule>();
  readonly #regions = new Map<string, number>();
  readonly #fasta = new Map<string, string>();
  readonly #emitted = new Set<string>();
  readonly #unknown = new Set<string>();
  readonly #exclude: ReadonlySet<string>;
  readonly #resolve: (ref: string) => string;
  /** Transcripts waiting for their exons (NCBI writes mRNA as one row spanning the gene; exons are children). */
  readonly #pending = new Map<string, { f: GffFeature; m: Molecule; provenance: Provenance; exons: GffRow[]; touched: number }>();
  #counter = 0;

  constructor(sink: Sink, options: Gff3Options = {}) {
    this.#sink = sink;
    this.#registry = options.registry ?? new NamespaceRegistry();
    this.#ctx = ingestContext(this.#registry, this.#units);
    this.#toRef = options.seqidToRef ?? ((seqid) => accessionRef(seqid, this.#registry));
    this.#source = new ChainedSource(this.#residues, options.source);
    this.#exclude = options.excludeAnnotations ?? DEFAULT_EXCLUDED_ANNOTATIONS;
    this.#resolve = options.resolveRef ?? ((r) => r);
    if (options.file) this.#base.file = options.file;
  }

  sequenceRegion(seqid: string, length: number): void {
    this.#regions.set(ownString(seqid), length);
    const m = this.#molecules.get(seqid);
    if (m && m.length === undefined) m.length = length;
  }

  /** Residues of a molecule (e.g. from a `##FASTA` section); must be given before its features. */
  residues(seqid: string, residues: string): void {
    this.#fasta.set(seqid, residues);
  }

  markCircular(seqid: string): void {
    const m = this.#molecule(seqid);
    if (m) m.circular = true;
  }

  #molecule(seqid: string): Molecule | undefined {
    let m = this.#molecules.get(seqid);
    if (m) return m;
    const found = this.#toRef(seqid);
    if (!found) return undefined;
    const ref = ownString(found);
    m = { ref, circular: false };
    const length = this.#regions.get(seqid);
    if (length !== undefined) m.length = length;
    this.#molecules.set(ownString(seqid), m);
    this.#units.set(ref, "nt");
    const residues = this.#fasta.get(seqid);
    if (residues) this.#residues.add(ref, residues);
    return m;
  }

  #sequence(record: SequenceRecord): void {
    if (this.#emitted.has(record.ref)) return;
    this.#emitted.add(record.ref);
    this.#sink.sequence(record);
  }

  feature(f: GffFeature): void {
    const m = this.#molecule(f.seqid);
    if (!m) {
      if (!this.#unknown.has(f.seqid)) this.#sink.warning(`seqid '${f.seqid}' is not a known accession; supply seqidToRef`);
      this.#unknown.add(ownString(f.seqid));
      return;
    }
    const provenance: Provenance = { ...this.#base, record: f.seqid, feature: `${f.type} ${f.id ?? `line ${f.rows[0]!.line}`}` };
    if (f.type === "region") {
      if (attr(f, "Is_circular") === "true") m.circular = true;
      this.#sequence(moleculeRecord(f, m, this.#fasta.get(f.seqid), provenance));
      return;
    }
    if (ALIGNMENT_TYPES.has(f.type) || attr(f, "Target") !== undefined) {
      this.#alignment(f, m.ref, provenance);
      return;
    }

    this.#counter++;
    // Transcripts: RNA feature types, or any other feature naming a transcript (Ensembl gene segments, ...);
    // parts of a transcript also carry transcript_id in NCBI GFF3 and are not transcripts themselves.
    const namesTranscript = attr(f, "transcript_id") !== undefined && !TRANSCRIPT_PARTS.has(f.type);
    if ((RNA_TYPES.has(f.type) || namesTranscript) && f.id !== undefined) {
      this.#pending.set(`${f.seqid}\t${f.id}`, { f, m, provenance, exons: [], touched: this.#counter });
      this.#flushTranscripts(f.seqid);
      return;
    }
    if (f.type === "exon") {
      for (const parent of f.attributes.Parent ?? []) {
        const p = this.#pending.get(`${f.seqid}\t${parent}`);
        if (!p) continue;
        p.exons.push(...f.rows);
        p.touched = this.#counter;
      }
    }

    const loc = featureLocation(f, m, this.#sink);
    if (loc) {
      if (!this.#exclude.has(f.type)) {
        this.#sink.annotation({ location: formatLocationId(loc, this.#ctx), type: f.type, attributes: f.attributes, extent: extent(loc), provenance });
      }
      if (f.type === "CDS") this.#cds(f, loc, m.ref, provenance);
    }
    this.#flushTranscripts(f.seqid);
  }

  /** Emit transcripts whose exons are complete: not touched for a while, or on another molecule. */
  #flushTranscripts(seqid: string | undefined, all = false): void {
    for (const [key, p] of this.#pending) {
      if (!all && p.f.seqid === seqid && this.#counter - p.touched < 1000) continue;
      this.#pending.delete(key);
      this.#transcript(p);
    }
  }

  /** A transcript located by its exons (or by its own rows when it has none). */
  #transcript(p: { f: GffFeature; m: Molecule; provenance: Provenance; exons: GffRow[] }): void {
    const { f, m, provenance } = p;
    const modelled: GffFeature = p.exons.length ? { ...f, rows: p.exons } : f;
    const loc = featureLocation(modelled, m, this.#sink);
    if (!loc) return;
    if (!this.#exclude.has(f.type)) {
      this.#sink.annotation({ location: formatLocationId(loc, this.#ctx), type: f.type, attributes: f.attributes, extent: extent(loc), provenance });
    }
    const transcript = attributeRef(attr(f, "transcript_id"), attr(f, "version"), this.#registry, this.#resolve);
    if (!transcript) return;
    this.#units.set(transcript, "nt");
    this.#sink.edge({
      kind: "annotation",
      from: transcript,
      to: m.ref,
      blocks: [...mappingFromLocation(transcript, loc).blocks],
      location: formatLocationId(loc, this.#ctx),
      attributes: { type: f.type, ...(attr(f, "exception") !== undefined && { exception: attr(f, "exception")! }) },
      provenance,
      validation: expectMismatch(validateTranscript(loc, transcript, this.#source, this.#ctx), attr(f, "exception")),
    });
  }

  finish(): void {
    this.#flushTranscripts(undefined, true);
    for (const [seqid, m] of this.#molecules) {
      this.#sequence(moleculeRecord(undefined, m, this.#fasta.get(seqid), { ...this.#base, record: seqid }));
    }
  }

  #cds(f: GffFeature, loc: Location, ref: string, provenance: Provenance): void {
    const st = { registry: this.#registry, units: this.#units, ctx: this.#ctx, source: this.#source };
    // Ensembl writes the protein version as `version=` on CDS rows; NCBI writes it in protein_id.
    const protein = attributeRef(attr(f, "protein_id"), attr(f, "version"), st.registry, this.#resolve);
    if (!protein) {
      // Pseudogene CDSs have no product; gene segments (e.g. IGKC) have no protein record of their own.
      if (attr(f, "pseudo") !== "true") this.#sink.warning(`${f.seqid}: CDS ${f.id ?? ""}: no protein_id; no protein edge`);
      return;
    }
    st.units.set(protein, "aa");
    // Phase of the 5'-most row gives the codon start.
    const first = loc.segments[0]!;
    const firstRow = f.rows.find((r) => (first.strand === 1 ? r.start - 1 === first.start : r.end === first.end)) ?? f.rows[0]!;
    const codonStart = (firstRow.phase ?? 0) + 1;
    const table = Number(attr(f, "transl_table") ?? 1);
    // Published protein residues (e.g. --fasta protein.faa) make validation exact; GFF3 has no translation. Their
    // length also replaces the inferred one (Ensembl CDSs without a stop codon would otherwise lose a residue).
    const protLength = st.source.length?.(protein);
    const inferred = inferAaLength(loc, codonStart);
    // Ensembl proteins of 5'-incomplete CDSs start with an X for the incomplete first codon (one residue longer).
    const leadingPartialCodon = codonStart > 1 && st.source.get(protein, 0, 1) === "X";
    const fits = protLength !== undefined && (leadingPartialCodon ? 3 * protLength - 3 + codonStart - 1 : 3 * protLength + codonStart - 1) <= cdsLength(loc);
    const aaLength = fits ? protLength! : inferred;
    let mapping;
    try {
      mapping = cdsMapping({ protein, cds: loc, codonStart, aaLength, leadingPartialCodon });
    } catch (e) {
      this.#sink.warning(`${f.seqid}: CDS ${f.id ?? ""}: ${(e as Error).message}`);
      return;
    }
    const attributes: Record<string, string> = {
      codonStart: String(codonStart),
      translTable: String(table),
      aaLength: String(aaLength),
      aaLengthSource: aaLength === protLength ? "protein sequence" : "inferred",
      ...(leadingPartialCodon && { leadingPartialCodon: "true" }),
    };
    for (const name of ["gene", "product", "exception"]) {
      const v = attr(f, name);
      if (v !== undefined) attributes[name] = v;
    }
    const validation = validateCds({
      cds: loc,
      mapping,
      codonStart,
      table,
      aaLength,
      ...(protLength !== undefined && { translation: st.source.get(protein, 0, protLength)! }),
      translExcept: f.attributes.transl_except ?? [],
      leadingPartialCodon,
      tableGiven: attr(f, "transl_table") !== undefined,
      outer: ref,
      ...(attr(f, "exception") !== undefined && { exception: attr(f, "exception")! }),
      ctx: st.ctx,
      source: st.source,
    });
    if (validation.table !== undefined) attributes.translTable = String(validation.table);
    delete validation.table;
    this.#sink.edge({
      kind: "annotation",
      from: protein,
      to: ref,
      blocks: [...mapping.blocks],
      location: formatLocationId(loc, st.ctx),
      attributes,
      provenance,
      validation,
    });
    // Published residues (--fasta protein.faa) give the checksums used to identify the protein across databases.
    const residues = protLength !== undefined ? st.source.get(protein, 0, protLength) : undefined;
    this.#sequence({
      ref: protein,
      moltype: "protein",
      unit: "aa",
      length: protLength ?? aaLength,
      provenance,
      ...(residues !== undefined && checksums(residues)),
    });
  }

  /**
   * Target + Gap alignment -> blocks from the target (e.g. a RefSeq transcript) to the genome.
   * Gap operations follow the genome in the direction of the row's strand (spec-ingest §4).
   */
  #alignment(f: GffFeature, genome: string, provenance: Provenance): void {
      const st = { registry: this.#registry, units: this.#units, source: this.#source };
    const blocks: Block[] = [];
    let target: string | undefined;
    for (const row of f.rows) {
      const [tid, tstart, tend, tstrand = "+"] = (row.attributes.Target?.[0] ?? "").split(/\s+/);
      const ref = tid && accessionRef(tid, st.registry);
      if (!ref || tstart === undefined || tend === undefined) {
        this.#sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: unusable Target`);
        return;
      }
      if (target !== undefined && ref !== target) {
        this.#sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: rows target different sequences`);
        return;
      }
      if (row.strand === "." || row.strand === "?" || (tstrand !== "+" && tstrand !== "-")) {
        this.#sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: unstranded alignment rows are not supported`);
        return;
      }
      target = ref;
      // Gap operations follow the genome in the direction of the row's strand (ascending on '+' rows, descending
      // on '-' rows), whichever strand the Target is on. Verified on NCBI GRCh38: '-' row / '+' Target NM_012234.7
      // (4036/4041 identical in the gapped row, num_mismatch=5); '+' row / '-' Target NG_162589.1 and NG_043318.1
      // (all bases identical; the opposite order gives 150/270 and 204/301).
      const gAsc = row.strand === "+";
      const tAsc = tstrand === "+";
      let g = gAsc ? row.start - 1 : row.end;
      let t = tAsc ? Number(tstart) - 1 : Number(tend);
      const ops = (row.attributes.Gap?.[0] ?? `M${row.end - row.start + 1}`).trim().split(/\s+/);
      for (const op of ops) {
        const n = Number(op.slice(1));
        switch (op[0]) {
          case "M":
            blocks.push({ srcRef: ref, src: tAsc ? t : t - n, tgtRef: genome, tgt: gAsc ? g : g - n, len: n, rev: gAsc !== tAsc });
            g += gAsc ? n : -n;
            t += tAsc ? n : -n;
            break;
          case "I":
            t += tAsc ? n : -n;
            break;
          case "D":
            g += gAsc ? n : -n;
            break;
          default:
            this.#sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: Gap operation '${op}' is not supported`);
            return;
        }
      }
      if (g !== (gAsc ? row.end : row.start - 1)) {
        this.#sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: Gap does not span the row (line ${row.line})`);
        return;
      }
      if (t !== (tAsc ? Number(tend) : Number(tstart) - 1)) {
        this.#sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: Gap does not span the Target range (line ${row.line})`);
        return;
      }
    }
    if (!target) return;
    st.units.set(target, "nt");

    const attributes: Record<string, string> = { type: f.type };
    for (const name of ["identity", "pct_coverage", "pct_identity_gap", "num_mismatch", "gap_count"]) {
      const v = attr(f, name);
      if (v !== undefined) attributes[name] = v;
    }
    const edge: Edge = { kind: "alignment", from: target, to: genome, blocks, attributes, provenance, validation: { status: "skipped" } };
    const identity = alignmentIdentity(blocks, st.source);
    if (!identity) {
      edge.validation = { status: "skipped", detail: "sequences not available" };
    } else {
      // Coordinates are judged by our own count: a misread Gap drops identity to ~25%, while real sequence differences
      // (alternate haplotypes, RefSeq corrections) keep it high. NCBI's num_mismatch is recorded, not trusted blindly:
      // it can differ from the count on the released sequences (e.g. NM_001291281.3: 3 counted, num_mismatch=4).
      const mismatches = identity.aligned - identity.identical;
      const reported = attributes.num_mismatch === undefined ? undefined : Number(attributes.num_mismatch);
      edge.validation = {
        status: identity.identical >= identity.aligned * MIN_ALIGNMENT_IDENTITY ? "ok" : "mismatch",
        basis: "full",
        detail:
          `${identity.identical}/${identity.aligned} aligned bases identical` +
          (reported === undefined ? "" : `; num_mismatch=${reported}`) +
          (reported !== undefined && reported !== mismatches ? ` (counted ${mismatches})` : ""),
      };
    }
    this.#sink.edge(edge);
  }
}

function moleculeRecord(
  region: GffFeature | undefined,
  m: Molecule,
  residues: string | undefined,
  provenance: Provenance,
): SequenceRecord {
  const mol = region ? attr(region, "mol_type") ?? "" : "";
  const seq: SequenceRecord = {
    ref: m.ref,
    moltype: /RNA/i.test(mol) ? "RNA" : "DNA",
    unit: "nt",
    length: m.length ?? residues?.length ?? region?.rows[0]?.end ?? 0,
    topology: m.circular ? "circular" : "linear",
    provenance,
  };
  const taxon = region?.attributes.Dbxref?.find((x) => x.startsWith("taxon:"));
  if (taxon) seq.taxon = Number(taxon.slice(6));
  if (residues) Object.assign(seq, checksums(residues));
  return seq;
}

/**
 * Rows -> Location in traversal order. Rows extending past the end of a circular molecule
 * (NCBI writes origin-spanning features as e.g. 16024..17145 on a 16569 nt genome) are wrapped.
 */
function featureLocation(f: GffFeature, m: Molecule, sink: Sink): Location | undefined {
  const strands = new Set(f.rows.map((r) => r.strand));
  if (strands.size > 1) {
    sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: rows on different strands; feature skipped`);
    return undefined;
  }
  const strand: 1 | -1 = f.rows[0]!.strand === "-" ? -1 : 1;
  const rows = [...f.rows].sort((a, b) => a.start - b.start);
  let segments: Segment[] = rows.map((r) => rowSegment(r, m.ref, strand));
  if (strand === -1) segments.reverse();
  if (m.circular && m.length !== undefined) {
    const len = m.length;
    segments = segments.flatMap((s) => wrap(s, len));
  } else if (m.length !== undefined && segments.some((s) => s.end > m.length!)) {
    sink.warning(`${f.seqid}: ${f.type} ${f.id ?? ""}: extends past the end of a linear sequence; feature skipped`);
    return undefined;
  }
  return { outer: m.ref, kind: "join", segments };
}

function rowSegment(r: GffRow, ref: string, strand: 1 | -1): Segment {
  const seg: Segment = { ref, start: r.start - 1, end: r.end, strand };
  // NCBI marks partial ends with start_range=.,N / end_range=N,.
  if (r.attributes.start_range?.[0] === ".") seg.fuzzyLow = true;
  if (r.attributes.end_range?.[1] === ".") seg.fuzzyHigh = true;
  return seg;
}

function wrap(s: Segment, len: number): Segment[] {
  if (s.end <= len) return [s];
  if (s.start >= len) return [{ ...s, start: s.start - len, end: s.end - len }];
  const { fuzzyLow, fuzzyHigh, ...plain } = s;
  const head: Segment = { ...plain, end: len };
  const tail: Segment = { ...plain, start: 0, end: s.end - len };
  if (fuzzyLow) head.fuzzyLow = true;
  if (fuzzyHigh) tail.fuzzyHigh = true;
  return s.strand === 1 ? [head, tail] : [tail, head];
}

