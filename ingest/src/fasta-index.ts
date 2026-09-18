// Random access to large FASTA files through a samtools-compatible .fai index (scaling §3, P4).
import { closeSync, existsSync, openSync, readFileSync, readSync, writeFileSync } from "node:fs";
import type { SequenceSource } from "./sequence.ts";

export interface FaiEntry {
  name: string;
  length: number;
  /** Byte offset of the first residue. */
  offset: number;
  lineBases: number;
  lineBytes: number;
}

/** Scan an uncompressed FASTA and compute its index. Lines of a record must share one width (except the last). */
export function buildFai(path: string): FaiEntry[] {
  const fd = openSync(path, "r");
  const entries: FaiEntry[] = [];
  const chunk = Buffer.allocUnsafe(1 << 24);
  let carry = Buffer.alloc(0);
  let filePos = 0; // byte offset of carry[0]
  let current: FaiEntry | undefined;
  let lastLineShort = false; // a line shorter than lineBases was seen in the current record

  const line = (buf: Buffer, start: number, end: number, absolute: number) => {
    // [start, end) excludes the '\n'; absolute = file offset of buf[start]
    const hasCr = end > start && buf[end - 1] === 13;
    const bases = end - start - (hasCr ? 1 : 0);
    if (buf[start] === 62 /* '>' */) {
      if (current) entries.push(current);
      const header = buf.toString("latin1", start + 1, end - (hasCr ? 1 : 0));
      current = { name: header.trim().split(/\s+/)[0] ?? "", length: 0, offset: absolute + (end - start) + 1, lineBases: 0, lineBytes: 0 };
      lastLineShort = false;
      return;
    }
    if (!current || bases === 0) return;
    if (current.lineBases === 0) {
      current.lineBases = bases;
      current.lineBytes = end - start + 1;
    } else if (lastLineShort || bases > current.lineBases) {
      throw new Error(`${path}: irregular line length in ${current.name}; cannot index`);
    }
    if (bases < current.lineBases) lastLineShort = true;
    current.length += bases;
  };

  for (;;) {
    const n = readSync(fd, chunk, 0, chunk.length, null);
    if (n === 0) break;
    const buf = carry.length ? Buffer.concat([carry, chunk.subarray(0, n)]) : chunk.subarray(0, n);
    let start = 0;
    for (let nl = buf.indexOf(10, start); nl >= 0; nl = buf.indexOf(10, start)) {
      line(buf, start, nl, filePos + start);
      start = nl + 1;
    }
    carry = Buffer.from(buf.subarray(start));
    filePos += start;
  }
  if (carry.length) line(carry, 0, carry.length, filePos);
  if (current) entries.push(current);
  closeSync(fd);
  return entries;
}

export function writeFai(path: string, entries: FaiEntry[]): void {
  writeFileSync(path, entries.map((e) => `${e.name}\t${e.length}\t${e.offset}\t${e.lineBases}\t${e.lineBytes}\n`).join(""));
}

export function readFai(path: string): FaiEntry[] {
  return readFileSync(path, "utf8")
    .split("\n")
    .filter(Boolean)
    .map((l) => {
      const [name, length, offset, lineBases, lineBytes] = l.split("\t");
      return { name: name!, length: Number(length), offset: Number(offset), lineBases: Number(lineBases), lineBytes: Number(lineBytes) };
    });
}

/** Load `<fasta>.fai`, building (and writing) it when missing. */
export function loadOrBuildFai(fasta: string): FaiEntry[] {
  const fai = `${fasta}.fai`;
  if (existsSync(fai)) return readFai(fai);
  const entries = buildFai(fasta);
  writeFai(fai, entries);
  return entries;
}

/** Residues read on demand from an indexed FASTA; `nameToRef` maps FASTA ids to sequence keys. */
export class FaiSequenceSource implements SequenceSource {
  readonly #fd: number;
  readonly #entries = new Map<string, FaiEntry>();

  constructor(fasta: string, nameToRef: (name: string) => string | undefined, entries = loadOrBuildFai(fasta)) {
    this.#fd = openSync(fasta, "r");
    for (const e of entries) {
      const ref = nameToRef(e.name);
      if (ref) this.#entries.set(ref, e);
    }
  }

  refs(): string[] {
    return [...this.#entries.keys()];
  }

  length(ref: string): number | undefined {
    return this.#entries.get(ref)?.length;
  }

  get(ref: string, start: number, end: number): string | undefined {
    const e = this.#entries.get(ref);
    if (!e || start < 0 || end > e.length || start > end) return undefined;
    if (start === end) return "";
    const byte = (p: number) => e.offset + Math.floor(p / e.lineBases) * e.lineBytes + (p % e.lineBases);
    const from = byte(start);
    const to = byte(end - 1) + 1;
    const buf = Buffer.allocUnsafe(to - from);
    readSync(this.#fd, buf, 0, buf.length, from);
    let out = "";
    let runStart = 0;
    for (let i = 0; i < buf.length; i++) {
      const c = buf[i];
      if (c === 10 || c === 13) {
        if (i > runStart) out += buf.toString("latin1", runStart, i);
        runStart = i + 1;
      }
    }
    if (runStart < buf.length) out += buf.toString("latin1", runStart);
    return out.toUpperCase();
  }

  close(): void {
    closeSync(this.#fd);
  }
}
