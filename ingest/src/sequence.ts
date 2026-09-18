// Residue-level helpers: sequence access, reverse complement, translation, refget digests.
import { createHash } from "node:crypto";
import type { Location } from "@togocoord/core";
import { GENETIC_CODES } from "./genetic-codes.ts";

/** Random access to residues; coordinates are 0-based half-open in residues. */
export interface SequenceSource {
  get(ref: string, start: number, end: number): string | undefined;
  /** Length in residues, when the whole sequence is available. */
  length?(ref: string): number | undefined;
}

/** In-memory source. A molecule may be given as several slices, each starting at `offset` (0-based). */
export class MemorySequenceSource implements SequenceSource {
  readonly #seqs = new Map<string, Array<{ offset: number; residues: string }>>();

  add(ref: string, residues: string, offset = 0): this {
    const slices = this.#seqs.get(ref) ?? [];
    slices.push({ offset, residues: residues.toUpperCase() });
    this.#seqs.set(ref, slices);
    return this;
  }

  has(ref: string): boolean {
    return this.#seqs.has(ref);
  }

  /** Length of a sequence given in full (a single slice at offset 0). */
  length(ref: string): number | undefined {
    const slices = this.#seqs.get(ref);
    return slices?.length === 1 && slices[0]!.offset === 0 ? slices[0]!.residues.length : undefined;
  }

  get(ref: string, start: number, end: number): string | undefined {
    const s = this.#seqs.get(ref)?.find((x) => start >= x.offset && end <= x.offset + x.residues.length);
    return s?.residues.slice(start - s.offset, end - s.offset);
  }
}

const COMPLEMENT: Record<string, string> = {
  A: "T", C: "G", G: "C", T: "A", U: "A", R: "Y", Y: "R", K: "M", M: "K", S: "S", W: "W", B: "V", V: "B", D: "H", H: "D", N: "N",
};

export function reverseComplement(seq: string): string {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) out += COMPLEMENT[seq[i]!] ?? "N";
  return out;
}

/** Residues of a nucleotide location in traversal order, or undefined if any part is unavailable. */
export function extract(loc: Location, source: SequenceSource): string | undefined {
  let out = "";
  for (const seg of loc.segments) {
    const part = source.get(seg.ref, seg.start, seg.end);
    if (part === undefined) return undefined;
    out += seg.strand === 1 ? part : reverseComplement(part);
  }
  return out;
}

const BASE: Record<string, number> = { T: 0, U: 0, C: 1, A: 2, G: 3 };

function codonIndex(codon: string): number | undefined {
  const [a, b, c] = [BASE[codon[0]!], BASE[codon[1]!], BASE[codon[2]!]];
  return a === undefined || b === undefined || c === undefined ? undefined : 16 * a + 4 * b + c;
}

export function hasGeneticCode(table: number): boolean {
  return table in GENETIC_CODES;
}

/** Translate complete codons; ambiguous codons give 'X', stops give '*'. */
export function translate(nt: string, table = 1): string {
  const code = GENETIC_CODES[table];
  if (!code) throw new Error(`unknown genetic code ${table}`);
  let aa = "";
  for (let i = 0; i + 3 <= nt.length; i += 3) {
    const idx = codonIndex(nt.slice(i, i + 3).toUpperCase());
    aa += idx === undefined ? "X" : code.aa[idx];
  }
  return aa;
}

/** Whether `codon` is an initiation codon in `table` (NCBI sncbieaa 'M'). */
export function isStartCodon(codon: string, table = 1): boolean {
  const idx = codonIndex(codon.toUpperCase());
  return idx !== undefined && GENETIC_CODES[table]?.starts[idx] === "M";
}

/** GA4GH refget sequence digest: `SQ.` + base64url(sha512(residues.upper())[0:24]). */
export function refgetDigest(residues: string): string {
  const hash = createHash("sha512").update(residues.toUpperCase()).digest();
  return `SQ.${hash.subarray(0, 24).toString("base64url")}`;
}

/** INSDC /transl_except amino acid names to one-letter codes ('*' for TERM). */
export const AA_ABBREVIATIONS: Readonly<Record<string, string>> = {
  Ala: "A", Arg: "R", Asn: "N", Asp: "D", Cys: "C", Gln: "Q", Glu: "E", Gly: "G", His: "H", Ile: "I",
  Leu: "L", Lys: "K", Met: "M", Phe: "F", Pro: "P", Ser: "S", Thr: "T", Trp: "W", Tyr: "Y", Val: "V",
  Sec: "U", Pyl: "O", Asx: "B", Glx: "Z", Xle: "J", Xaa: "X", OTHER: "X", TERM: "*",
};
