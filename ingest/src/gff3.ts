// GFF3 parser (https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).
// Line-level pieces (parseGffLine, FeatureGrouper) are shared by the in-memory and streaming readers.
import { parseFasta } from "./fasta.ts";

export interface GffRow {
  seqid: string;
  source: string;
  type: string;
  /** 1-based closed, as written (may exceed the sequence length on circular molecules). */
  start: number;
  end: number;
  score?: number;
  strand: "+" | "-" | "." | "?";
  phase?: 0 | 1 | 2;
  attributes: Record<string, string[]>;
  line: number;
}

/** Rows sharing an ID form one (possibly discontinuous) feature. */
export interface GffFeature {
  id?: string;
  type: string;
  seqid: string;
  rows: GffRow[];
  attributes: Record<string, string[]>;
}

export interface Gff3Document {
  /** `##sequence-region` lengths (end coordinate). */
  regions: Map<string, number>;
  features: GffFeature[];
  /** `##FASTA` section, if any. */
  sequences: Map<string, string>;
}

/** Parse a whole GFF3 text. Rows sharing an ID are grouped wherever they occur. */
export function parseGff3(text: string): Gff3Document {
  const regions = new Map<string, number>();
  const features: GffFeature[] = [];
  const grouper = new FeatureGrouper(Infinity);
  const lines = text.split(/\r?\n/);
  let fastaFrom = -1;
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i]!;
    if (line.startsWith("##FASTA")) {
      fastaFrom = i + 1;
      break;
    }
    const region = sequenceRegion(line);
    if (region) regions.set(region[0], region[1]);
    const row = parseGffLine(line, i + 1);
    if (row) features.push(...grouper.push(row));
  }
  features.push(...grouper.end());
  const sequences = fastaFrom < 0 ? new Map<string, string>() : parseFasta(lines.slice(fastaFrom).join("\n"));
  return { regions, features, sequences };
}

/** `##sequence-region seqid start end` -> [seqid, end]. */
export function sequenceRegion(line: string): [string, number] | undefined {
  if (!line.startsWith("##sequence-region")) return undefined;
  const [, seqid, , end] = line.trim().split(/\s+/);
  return seqid && end ? [decode(seqid), Number(end)] : undefined;
}

/** One feature line, or undefined for blank lines, comments and directives. */
export function parseGffLine(line: string, lineNo: number): GffRow | undefined {
  if (line === "" || line.startsWith("#")) return undefined;
  const cols = line.split("\t");
  if (cols.length !== 9) throw new Error(`GFF3 line ${lineNo}: expected 9 columns, found ${cols.length}`);
  const [seqid, source, type, start, end, score, strand, phase, attrs] = cols as [string, string, string, string, string, string, string, string, string];
  const row: GffRow = {
    seqid: decode(seqid),
    source,
    type,
    start: Number(start),
    end: Number(end),
    strand: (["+", "-", ".", "?"].includes(strand) ? strand : ".") as GffRow["strand"],
    attributes: parseAttributes(attrs),
    line: lineNo,
  };
  if (score !== ".") row.score = Number(score);
  if (phase === "0" || phase === "1" || phase === "2") row.phase = Number(phase) as 0 | 1 | 2;
  return row;
}

/**
 * Groups rows into features by (seqid, ID), emitting features in order of first appearance.
 * With a finite `window`, at most that many features are held open; this is enough for files where rows of one
 * feature are adjacent (as in NCBI RefSeq GFF3). A row whose ID was already emitted is reported through
 * `violations` and starts a new, separate feature.
 */
export class FeatureGrouper {
  readonly #window: number;
  readonly #open = new Map<string, GffFeature>();
  /** Recently emitted keys (ring buffer), to detect non-adjacent rows cheaply. */
  readonly #recent = new Set<string>();
  readonly #ring: Array<string | undefined> = new Array(100_000);
  #ringPos = 0;
  #anonymous = 0;
  violations: string[] = [];

  constructor(window = 1000) {
    this.#window = window;
  }

  push(row: GffRow): GffFeature[] {
    const id = row.attributes.ID?.[0];
    const key = id === undefined ? `\0${this.#anonymous++}` : `${row.seqid}\t${id}`;
    const open = this.#open.get(key);
    if (open) {
      open.rows.push(row);
      return [];
    }
    if (id !== undefined && this.#recent.has(key)) {
      this.violations.push(`line ${row.line}: rows of ${row.seqid} ${id} are not adjacent; grouped separately`);
    }
    const feature: GffFeature = { type: row.type, seqid: row.seqid, rows: [row], attributes: row.attributes };
    if (id !== undefined) feature.id = id;
    this.#open.set(key, feature);
    return this.#open.size > this.#window ? this.#flushOldest() : [];
  }

  end(): GffFeature[] {
    const out = [...this.#open.values()];
    this.#open.clear();
    return out;
  }

  #flushOldest(): GffFeature[] {
    const [openKey, feature] = this.#open.entries().next().value!;
    this.#open.delete(openKey);
    if (!openKey.startsWith("\0")) {
      // Stored beyond the life of the input chunk: copy (see ownString in common.ts).
      const key = JSON.parse(JSON.stringify(openKey)) as string;
      const evicted = this.#ring[this.#ringPos];
      if (evicted !== undefined) this.#recent.delete(evicted);
      this.#ring[this.#ringPos] = key;
      this.#ringPos = (this.#ringPos + 1) % this.#ring.length;
      this.#recent.add(key);
    }
    return [feature];
  }
}

function decode(s: string): string {
  return s.includes("%") ? decodeURIComponent(s) : s;
}

function parseAttributes(text: string): Record<string, string[]> {
  const attrs: Record<string, string[]> = {};
  if (text === "." || text === "") return attrs;
  for (const pair of text.split(";")) {
    if (!pair) continue;
    const eq = pair.indexOf("=");
    if (eq < 0) continue;
    attrs[decode(pair.slice(0, eq))] = pair.slice(eq + 1).split(",").map(decode);
  }
  return attrs;
}

export function attr(f: GffFeature, name: string): string | undefined {
  return f.attributes[name]?.[0];
}
