// GenBank / GenPept flat file adapter (Core tier T0).
import {
  cdsMapping,
  formatLocationId,
  mappingFromLocation,
  NamespaceRegistry,
  parseLocationId,
  type CoordContext,
  type Location,
  type Unit,
} from "@togocoord/core";
import { accessionRef, ChainedSource, DEFAULT_EXCLUDED_ANNOTATIONS, DROPPED_ATTRIBUTES, extent, ingestContext, RNA_TYPES, type AdapterOptions } from "./common.ts";
import { parseGenBank, qualifier, qualifiers, type GbFeature, type GbRecord } from "./gbff.ts";
import { MemorySink, type Annotation, type Edge, type IngestResult, type Provenance, type SequenceRecord, type Sink } from "./model.ts";
import { MemorySequenceSource, refgetDigest } from "./sequence.ts";
import { expectMismatch, inferAaLength, validateCds, validateTranscript } from "./validate.ts";

export function ingestGenBank(text: string, options: AdapterOptions = {}): IngestResult {
  return ingestGenBankRecords(parseGenBank(text), options);
}

export function ingestGenBankRecords(records: GbRecord[], options: AdapterOptions = {}): IngestResult {
  const sink = new MemorySink();
  const ingestor = new GenBankIngestor(sink, options);
  for (const r of records) ingestor.record(r);
  return sink.result;
}

/** Record-at-a-time GenBank/GenPept adapter; residues are held only for the record being processed. */
export class GenBankIngestor {
  readonly #sink: Sink;
  readonly #options: AdapterOptions;
  readonly #registry: NamespaceRegistry;
  readonly #units = new Map<string, Unit>();
  readonly #ctx: CoordContext;
  readonly #emitted = new Set<string>();

  constructor(sink: Sink, options: AdapterOptions = {}) {
    this.#sink = sink;
    this.#options = options;
    this.#registry = options.registry ?? new NamespaceRegistry();
    this.#ctx = ingestContext(this.#registry, this.#units);
  }

  record(record: GbRecord): void {
    const acc = record.version ?? record.accession ?? record.name;
    const ref = accessionRef(acc, this.#registry);
    if (!ref) {
      this.#sink.warning(`${acc}: not an INSDC/RefSeq accession; record skipped`);
      return;
    }
    this.#units.set(ref, record.unit === "aa" ? "aa" : "nt");
    const residues = new MemorySequenceSource();
    if (record.sequence) residues.add(ref, record.sequence);
    const provenance: Provenance = { adapter: "gbff", record: ref.slice(ref.indexOf(":") + 1) };
    if (this.#options.file) provenance.file = this.#options.file;
    const st: State = {
      ctx: this.#ctx,
      registry: this.#registry,
      source: new ChainedSource(residues, this.#options.source),
      residues,
      sink: this.#sink,
      provenance,
      units: this.#units,
      exclude: this.#options.excludeAnnotations ?? DEFAULT_EXCLUDED_ANNOTATIONS,
      sequence: (s) => {
        if (this.#emitted.has(s.ref)) return;
        this.#emitted.add(s.ref);
        this.#sink.sequence(s);
      },
    };
    st.sequence(sequenceRecord(record, ref, provenance));
    for (const feature of record.features) {
      if (feature.key === "source") continue;
      ingestFeature(feature, record, ref, st);
    }
  }
}

function sequenceRecord(record: GbRecord, ref: string, provenance: Provenance): SequenceRecord {
  const moltype = record.unit === "aa" ? "protein" : /RNA/i.test(record.moltype ?? "") ? "RNA" : "DNA";
  const seq: SequenceRecord = { ref, moltype, unit: record.unit === "aa" ? "aa" : "nt", length: record.length, provenance };
  if (record.topology) seq.topology = record.topology;
  const src = record.features.find((f) => f.key === "source");
  const taxon = src && qualifiers(src, "db_xref").find((x) => x.startsWith("taxon:"));
  if (taxon) seq.taxon = Number(taxon.slice(6));
  const organism = src && qualifier(src, "organism");
  if (organism) seq.organism = organism;
  if (record.sequence) seq.digest = refgetDigest(record.sequence);
  return seq;
}

interface State {
  ctx: CoordContext;
  registry: NamespaceRegistry;
  source: ChainedSource;
  residues: MemorySequenceSource;
  sink: Sink;
  provenance: Provenance;
  units: Map<string, Unit>;
  exclude: ReadonlySet<string>;
  sequence: (s: SequenceRecord) => void;
}

function ingestFeature(f: GbFeature, record: GbRecord, ref: string, st: State): void {
  const label = `${f.key} ${f.location}`;
  const provenance = { ...st.provenance, feature: label };
  let loc: Location;
  try {
    loc = parseLocationId(`${ref}:${f.location}`, st.ctx);
  } catch (e) {
    st.sink.warning(`${provenance.record}: ${label}: ${(e as Error).message}`);
    return;
  }
  if (!st.exclude.has(f.key)) st.sink.annotation(annotation(f, loc, st.ctx, provenance));

  if (f.key === "CDS" && qualifier(f, "pseudo") === undefined) {
    if (record.unit === "aa") ingestCodedBy(f, record, ref, provenance, st);
    else ingestCds(f, loc, ref, provenance, st);
  } else if (RNA_TYPES.has(f.key)) {
    const tid = qualifier(f, "transcript_id");
    const transcript = tid && accessionRef(tid, st.registry);
    if (!transcript) return;
    st.units.set(transcript, "nt");
    st.sink.edge({
      kind: "annotation",
      from: transcript,
      to: ref,
      blocks: [...mappingFromLocation(transcript, loc).blocks],
      location: formatLocationId(loc, st.ctx),
      attributes: { type: f.key, ...(qualifier(f, "exception") !== undefined && { exception: qualifier(f, "exception") || "true" }) },
      provenance,
      validation: expectMismatch(validateTranscript(loc, transcript, st.source, st.ctx), qualifier(f, "exception")),
    });
  }
}

function annotation(f: GbFeature, loc: Location, ctx: CoordContext, provenance: Provenance): Annotation {
  const attributes: Record<string, string[]> = {};
  for (const [name, value] of f.qualifiers) {
    if (DROPPED_ATTRIBUTES.has(name)) continue;
    (attributes[name] ??= []).push(value === true ? "" : value);
  }
  return { location: formatLocationId(loc, ctx), type: f.key, attributes, extent: extent(loc), provenance };
}

function ingestCds(f: GbFeature, loc: Location, ref: string, provenance: Provenance, st: State): void {
  const pid = qualifier(f, "protein_id");
  const protein = pid && accessionRef(pid, st.registry);
  if (!protein) {
    st.sink.warning(`${provenance.record}: ${provenance.feature}: CDS without a usable /protein_id`);
    return;
  }
  st.units.set(protein, "aa");
  const codonStart = Number(qualifier(f, "codon_start") ?? 1);
  const table = Number(qualifier(f, "transl_table") ?? 1);
  const translation = qualifier(f, "translation");
  const aaLength = translation ? translation.length : inferAaLength(loc, codonStart);
  if (translation) st.residues.add(protein, translation);
  pushCdsEdge({ f, cds: loc, outer: ref, protein, codonStart, table, aaLength, translation, provenance, st });

  const seq: SequenceRecord = { ref: protein, moltype: "protein", unit: "aa", length: aaLength, provenance };
  if (translation) seq.digest = refgetDigest(translation);
  st.sequence(seq);
}

/** GenPept CDS: /coded_by gives the nucleotide location of this protein. */
function ingestCodedBy(f: GbFeature, record: GbRecord, protein: string, provenance: Provenance, st: State): void {
  const codedBy = qualifier(f, "coded_by");
  const acc = codedBy && /([A-Za-z0-9_]+\.\d+):/.exec(codedBy)?.[1];
  const outer = acc && accessionRef(acc, st.registry);
  if (!codedBy || !outer) {
    st.sink.warning(`${provenance.record}: CDS without a usable /coded_by`);
    return;
  }
  st.units.set(outer, "nt");
  let cds: Location;
  try {
    cds = parseLocationId(`${outer}:${codedBy.replace(/\s+/g, "")}`, st.ctx);
  } catch (e) {
    st.sink.warning(`${provenance.record}: coded_by ${codedBy}: ${(e as Error).message}`);
    return;
  }
  pushCdsEdge({
    f,
    cds,
    outer,
    protein,
    codonStart: Number(qualifier(f, "codon_start") ?? 1),
    table: Number(qualifier(f, "transl_table") ?? 1),
    aaLength: record.length,
    translation: record.sequence,
    provenance: { ...provenance, feature: `CDS coded_by ${codedBy}` },
    st,
  });
}

function pushCdsEdge(a: {
  f: GbFeature;
  cds: Location;
  outer: string;
  protein: string;
  codonStart: number;
  table: number;
  aaLength: number;
  translation: string | undefined;
  provenance: Provenance;
  st: State;
}): void {
  const { st } = a;
  let mapping;
  try {
    mapping = cdsMapping({ protein: a.protein, cds: a.cds, codonStart: a.codonStart, aaLength: a.aaLength });
  } catch (e) {
    st.sink.warning(`${a.provenance.record}: ${a.provenance.feature}: ${(e as Error).message}`);
    return;
  }
  const attributes: Record<string, string> = {
    codonStart: String(a.codonStart),
    translTable: String(a.table),
    aaLength: String(a.aaLength),
    aaLengthSource: a.translation ? "translation" : "inferred",
  };
  for (const name of ["gene", "product", "exception"]) {
    const v = qualifier(a.f, name);
    if (v !== undefined) attributes[name] = v || "true";
  }
  if (qualifier(a.f, "ribosomal_slippage") !== undefined) attributes.ribosomalSlippage = "true";
  const edge: Edge = {
    kind: "annotation",
    from: a.protein,
    to: a.outer,
    blocks: [...mapping.blocks],
    location: formatLocationId(a.cds, st.ctx),
    attributes,
    provenance: a.provenance,
    validation: validateCds({
      cds: a.cds,
      mapping,
      codonStart: a.codonStart,
      table: a.table,
      aaLength: a.aaLength,
      ...(a.translation !== undefined && { translation: a.translation }),
      translExcept: qualifiers(a.f, "transl_except"),
      outer: a.outer,
      ...(qualifier(a.f, "exception") !== undefined && { exception: qualifier(a.f, "exception")! }),
      ctx: st.ctx,
      source: st.source,
    }),
  };
  st.sink.edge(edge);
}
