// Output model of adapters (design §5, §6.1) and the sinks they write to (scaling §3).
import { Mapping, type Block, type Unit } from "@togocoord/core";

export interface Provenance {
  adapter: "gbff" | "gff3" | "fasta" | "sifts";
  /** Source file or stream label. */
  file?: string;
  /** Record accession.version (GBFF) or seqid (GFF3). */
  record?: string;
  /** Feature key/type and its original location or ID. */
  feature?: string;
}

export interface Validation {
  status: "ok" | "mismatch" | "skipped";
  detail?: string;
  /** `full`: compared residue by residue with the actual sequence; `partial`: weaker checks only (no internal stop). */
  basis?: "full" | "partial";
}

export interface SequenceRecord {
  /** Internal key `namespace:accession.version`. */
  ref: string;
  moltype: "DNA" | "RNA" | "protein";
  unit: Unit;
  /** Length in residues. */
  length: number;
  topology?: "linear" | "circular";
  taxon?: number;
  organism?: string;
  /** GA4GH refget digest (`SQ.` + sha512t24u), when the residues are known. The key for sequence identity. */
  digest?: string;
  /** Lower-case hex MD5 of the upper-case residues (as used by UniParc and refget v1), for joining external data. */
  md5?: string;
  provenance: Provenance;
}

export interface Edge {
  kind: "annotation" | "alignment";
  from: string;
  to: string;
  /** Blocks in internal units (core §4.1); from -> to. */
  blocks: Block[];
  /** Canonical Location ID of the source feature on `to`, for display. */
  location?: string;
  attributes: Record<string, string>;
  provenance: Provenance;
  validation: Validation;
}

export interface Annotation {
  /** Canonical Location ID. */
  location: string;
  /** INSDC feature key or GFF3 type. */
  type: string;
  attributes: Record<string, string[]>;
  /** Bounding interval on the outer sequence, in internal units (for spatial indexes). */
  extent: { ref: string; start: number; end: number };
  provenance: Provenance;
}

/** Receiver of adapter output; implementations may stream to memory, JSON Lines or a database. */
export interface Sink {
  sequence(record: SequenceRecord): void;
  edge(edge: Edge): void;
  annotation(annotation: Annotation): void;
  warning(message: string): void;
}

export interface IngestResult {
  sequences: SequenceRecord[];
  edges: Edge[];
  annotations: Annotation[];
  warnings: string[];
}

/** Collects everything in memory (tests and small inputs). */
export class MemorySink implements Sink {
  readonly result: IngestResult = emptyResult();
  sequence(record: SequenceRecord): void {
    addSequence(this.result, record);
  }
  edge(edge: Edge): void {
    this.result.edges.push(edge);
  }
  annotation(annotation: Annotation): void {
    this.result.annotations.push(annotation);
  }
  warning(message: string): void {
    this.result.warnings.push(message);
  }
}

/** Core mapping of one or more edges (from -> to). */
export function edgeMapping(...edges: Edge[]): Mapping {
  return new Mapping(edges.flatMap((e) => e.blocks));
}

export function emptyResult(): IngestResult {
  return { sequences: [], edges: [], annotations: [], warnings: [] };
}

const sequenceKeys = new WeakMap<IngestResult, Set<string>>();

/** Append a sequence record unless one with the same key exists (O(1); results can hold 10^5+ proteins). */
export function addSequence(result: IngestResult, seq: SequenceRecord): boolean {
  let keys = sequenceKeys.get(result);
  if (!keys) sequenceKeys.set(result, (keys = new Set(result.sequences.map((s) => s.ref))));
  if (keys.has(seq.ref)) return false;
  keys.add(seq.ref);
  result.sequences.push(seq);
  return true;
}

/** Merge results, keeping the first record per sequence key. */
export function mergeResults(...results: IngestResult[]): IngestResult {
  const out = emptyResult();
  for (const r of results) {
    for (const s of r.sequences) addSequence(out, s);
    out.edges.push(...r.edges);
    out.annotations.push(...r.annotations);
    out.warnings.push(...r.warnings);
  }
  return out;
}
