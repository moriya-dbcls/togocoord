// Self-validation at ingest time (design §6.3).
import { formatLocationId, mapLocation, parseLocationId, type CoordContext, type Location, type Mapping } from "@togocoord/core";
import type { Validation } from "./model.ts";
import { AA_ABBREVIATIONS, extract, GENETIC_CODE_IDS, hasGeneticCode, isStartCodon, reverseComplement, translate, type SequenceSource } from "./sequence.ts";

export function fivePrimePartial(cds: Location): boolean {
  const s = cds.segments[0];
  return !!s && (s.strand === 1 ? !!s.fuzzyLow : !!s.fuzzyHigh);
}

export function threePrimePartial(cds: Location): boolean {
  const s = cds.segments.at(-1);
  return !!s && (s.strand === 1 ? !!s.fuzzyHigh : !!s.fuzzyLow);
}

export function cdsLength(cds: Location): number {
  return cds.segments.reduce((n, s) => n + s.end - s.start, 0);
}

/**
 * Protein length when no translation is given: all complete codons, minus the stop codon when the 3' end is
 * complete and ends on a codon boundary. A trailing partial codon on a 3'-complete CDS is taken to be a stop
 * completed by polyadenylation (e.g. vertebrate mitochondria).
 */
export function inferAaLength(cds: Location, codonStart: number): number {
  const coding = cdsLength(cds) - (codonStart - 1);
  const codons = Math.floor(coding / 3);
  if (threePrimePartial(cds)) return codons;
  return coding % 3 === 0 ? codons - 1 : codons;
}

export interface CdsCheck {
  cds: Location;
  /** Protein (codon units) -> nucleotide mapping built from `cds`. */
  mapping: Mapping;
  codonStart: number;
  table: number;
  aaLength: number;
  /** /translation (GBFF) or the protein record's residues; absent in GFF3. */
  translation?: string;
  /** Raw /transl_except values, e.g. `(pos:220..222,aa:Sec)`. */
  translExcept: string[];
  /** Sequence key that transl_except positions refer to. */
  outer: string;
  /** INSDC /exception (e.g. "annotated by transcript or proteomic data"): a mismatch is then expected. */
  exception?: string;
  /** The first residue stands for the incomplete first codon (Ensembl `X`); see cdsMapping. */
  leadingPartialCodon?: boolean;
  /** `table` came from the data (transl_table); otherwise other codes may be tried against the published protein. */
  tableGiven?: boolean;
  ctx: CoordContext;
  source: SequenceSource;
}

/** Translate the CDS from sequence and compare with the expected protein. */
export function validateCds(c: CdsCheck): Validation & { table?: number } {
  let v: Validation & { table?: number } = translateAndCompare(c);
  // Without an explicit genetic code (e.g. Ensembl GFF3 on mitochondria), take the code that reproduces the protein.
  if (v.status === "mismatch" && !c.tableGiven && c.translation !== undefined) {
    for (const table of GENETIC_CODE_IDS) {
      if (table === c.table) continue;
      const alt = translateAndCompare({ ...c, table });
      if (alt.status === "ok") {
        v = { ...alt, table, detail: `${alt.detail}; genetic code ${table} inferred (none given)` };
        break;
      }
    }
  }
  if (v.status === "mismatch" && c.exception) v.detail = `${v.detail}; expected: /exception="${c.exception}"`;
  return v;
}

function translateAndCompare(c: CdsCheck): Validation {
  if (!hasGeneticCode(c.table)) return { status: "skipped", detail: `unknown genetic code ${c.table}` };
  const nt = extract(c.cds, c.source);
  if (nt === undefined) return { status: "skipped", detail: "nucleotide sequence not available" };

  const coding = nt.slice(c.codonStart - 1);
  const aa = [...(c.leadingPartialCodon ? "X" : ""), ...translate(coding, c.table)];
  const notes: string[] = c.leadingPartialCodon ? ["first residue is the incomplete first codon (X)"] : [];

  // Initiation codons are read as Met: alternative starts of the genetic code, and any codon when the published
  // protein says M (non-AUG initiation, e.g. Ensembl GTG/ACG starts).
  const startCodon = coding.slice(0, 3);
  const completeStart = !fivePrimePartial(c.cds) && c.codonStart === 1 && aa.length > 0 && aa[0] !== "M";
  if (completeStart && (isStartCodon(startCodon, c.table) || c.translation?.[0]?.toUpperCase() === "M")) {
    aa[0] = "M";
    notes.push(isStartCodon(startCodon, c.table) ? `alternative start codon ${startCodon}` : `non-AUG initiation codon ${startCodon}`);
  }

  const inverse = c.mapping.inverse();
  for (const raw of c.translExcept) {
    const m = /^\(pos:(.+),aa:([A-Za-z]+)\)$/.exec(raw.replace(/\s+/g, ""));
    const code = m ? aaCode(m[2]!) : undefined;
    if (!m || code === undefined) {
      notes.push(`unparsed transl_except ${raw}`);
      continue;
    }
    if (code === "*") continue; // stop codons are outside the protein
    const hit = mapLocation(parseLocationId(`${c.outer}:${m[1]}`, c.ctx), inverse, c.ctx).pieces[0];
    if (!hit) {
      notes.push(`transl_except ${raw} is outside the CDS`);
      continue;
    }
    const residue = Math.floor(Math.min(hit.target.start, hit.target.end - 1) / 3);
    if (residue < aa.length) aa[residue] = code;
    notes.push(`transl_except ${m[2]} at residue ${residue + 1}`);
  }

  const product = aa.slice(0, c.aaLength).join("");
  if (product.length < c.aaLength) {
    return { status: "mismatch", detail: `CDS encodes ${product.length} complete codons, expected ${c.aaLength} aa` };
  }
  if (!threePrimePartial(c.cds)) {
    const next = aa[c.aaLength];
    if (next !== undefined && next !== "*") notes.push(`codon after the last residue translates to '${next}', not a stop`);
  }

  if (c.translation === undefined) {
    const stop = product.indexOf("*");
    return stop < 0
      ? { status: "ok", basis: "partial", detail: ["no internal stop codon", ...notes].join("; ") }
      : { status: "mismatch", basis: "partial", detail: [`internal stop codon at residue ${stop + 1}`, ...notes].join("; ") };
  }

  const expected = c.translation.toUpperCase();
  const diffs: number[] = [];
  for (let i = 0; i < Math.max(product.length, expected.length); i++) {
    const [a, b] = [product[i], expected[i]];
    if (a === "*" && (b === "U" || b === "O")) {
      // Stop codon recoded as selenocysteine / pyrrolysine without /transl_except (e.g. Ensembl GFF3).
      notes.push(`${b === "U" ? "selenocysteine" : "pyrrolysine"} at residue ${i + 1} (recoded stop codon)`);
      continue;
    }
    if (a !== b && a !== "X" && b !== "X") diffs.push(i + 1);
  }
  if (diffs.length === 0) return { status: "ok", basis: "full", detail: ["translation matches", ...notes].join("; ") };
  const shown = diffs.slice(0, 5).map((i) => `${i}:${product[i - 1] ?? "-"}/${expected[i - 1] ?? "-"}`);
  return { status: "mismatch", basis: "full", detail: [`${diffs.length} residue(s) differ (${shown.join(", ")})`, ...notes].join("; ") };
}

/** Record that a mismatch is expected when NCBI flags the feature with /exception (same wording as for CDS). */
export function expectMismatch(v: Validation, exception: string | undefined): Validation {
  return v.status === "mismatch" && exception !== undefined ? { ...v, detail: `${v.detail}; expected: /exception="${exception}"` } : v;
}

/**
 * Compare a genome-derived transcript model with the transcript's own sequence. Equal length with a few
 * substitutions keeps coordinates colinear (ok); a length difference means indels, so positions shift (mismatch).
 */
export function validateTranscript(loc: Location, transcript: string, source: SequenceSource, ctx: CoordContext): Validation {
  const modelLength = cdsLength(loc);
  const ownLength = source.length?.(transcript);
  const own = source.get(transcript, 0, ownLength ?? modelLength);
  if (own === undefined) return { status: "skipped", detail: "transcript sequence not available" };
  const model = extract(loc, source);
  if (model === undefined) return { status: "skipped", detail: "genome sequence not available" };
  // RefSeq transcripts often end in a poly(A) tail that the genome does not contain; positions before it align.
  const tail = ownLength !== undefined && ownLength > modelLength ? own.slice(modelLength) : "";
  const polyA = tail.length > 0 && [...tail].filter((c) => c === "A").length >= tail.length * 0.9;
  if (ownLength !== undefined && ownLength !== modelLength && !polyA) {
    return { status: "mismatch", basis: "full", detail: `genome model is ${modelLength} nt, ${transcript} is ${ownLength} nt (indels); use an alignment` };
  }
  const note = polyA ? `; ${tail.length} nt poly(A) tail not in the genome` : "";
  let substitutions = 0;
  for (let i = 0; i < model.length; i++) if (model[i] !== own[i]) substitutions++;
  if (substitutions === 0) return { status: "ok", basis: "full", detail: `identical to the transcript sequence${note}` };
  if (substitutions <= model.length * 0.05) {
    return { status: "ok", basis: "full", detail: `${substitutions} substitution(s) against ${transcript}; coordinates are colinear${note}` };
  }
  return { status: "mismatch", basis: "full", detail: `genome model differs from ${transcript} at ${substitutions} of ${model.length} positions (${formatLocationId(loc, ctx)})` };
}

/** Count identities over the aligned blocks of an alignment edge. */
export function alignmentIdentity(
  blocks: ReadonlyArray<{ srcRef: string; src: number; tgtRef: string; tgt: number; len: number; rev: boolean }>,
  source: SequenceSource,
): { identical: number; aligned: number } | undefined {
  let identical = 0;
  let aligned = 0;
  for (const b of blocks) {
    const s = source.get(b.srcRef, b.src, b.src + b.len);
    const t0 = source.get(b.tgtRef, b.tgt, b.tgt + b.len);
    if (s === undefined || t0 === undefined) return undefined;
    const t = b.rev ? reverseComplement(t0) : t0;
    for (let i = 0; i < b.len; i++) if (s[i] === t[i]) identical++;
    aligned += b.len;
  }
  return { identical, aligned };
}

const AA_BY_LOWER = new Map(Object.entries(AA_ABBREVIATIONS).map(([k, v]) => [k.toLowerCase(), v]));

/** transl_except amino acid, case-insensitively (GBFF writes `OTHER`, NCBI GFF3 writes `Other`). */
function aaCode(name: string): string | undefined {
  return AA_BY_LOWER.get(name.toLowerCase());
}

