// Shared adapter plumbing.
import { createContext, NamespaceRegistry, type CoordContext, type Location, type Unit } from "@togocoord/core";
import type { SequenceSource } from "./sequence.ts";

export interface AdapterOptions {
  /** Label recorded in provenance. */
  file?: string;
  registry?: NamespaceRegistry;
  /** Extra residues (e.g. transcript or genome sequences from FASTA) used for validation. */
  source?: SequenceSource;
  /** Feature types not stored as annotations (default: exon, which is implied by transcript edges). */
  excludeAnnotations?: ReadonlySet<string>;
}

export const DEFAULT_EXCLUDED_ANNOTATIONS: ReadonlySet<string> = new Set(["exon"]);

/**
 * A copy of `s` that does not share storage with a larger string. V8 substrings keep their parent alive, so a key
 * sliced from a 1 MB input chunk would pin the whole chunk for as long as the key is stored in a map.
 */
export function ownString(s: string): string {
  return JSON.parse(JSON.stringify(s)) as string;
}

/**
 * Accession -> internal key: Ensembl stable IDs (`ENST…`, `ENSP…`, also other species' `ENSMUSP…`) -> `ensembl`,
 * `XX_` prefixes -> `refseq`, others -> `insdc`.
 */
export function accessionRef(accession: string, registry: NamespaceRegistry): string | undefined {
  const namespace = /^ENS[A-Z]*[GTPER]\d{11}/.test(accession) ? "ensembl" : /^[A-Z]{2}_/.test(accession) ? "refseq" : "insdc";
  try {
    return ownString(registry.refKey(namespace, accession));
  } catch {
    return undefined;
  }
}

/** Context whose units come from the records being ingested, then from namespace defaults. */
export function ingestContext(registry: NamespaceRegistry, units: Map<string, Unit>): CoordContext {
  return createContext({ registry, units: (ref) => units.get(ref) });
}

/** First source that has the requested residues. */
export class ChainedSource implements SequenceSource {
  readonly #sources: SequenceSource[];
  constructor(...sources: Array<SequenceSource | undefined>) {
    this.#sources = sources.filter((s): s is SequenceSource => s !== undefined);
  }
  get(ref: string, start: number, end: number): string | undefined {
    for (const s of this.#sources) {
      const r = s.get(ref, start, end);
      if (r !== undefined) return r;
    }
    return undefined;
  }
  length(ref: string): number | undefined {
    for (const s of this.#sources) {
      const n = s.length?.(ref);
      if (n !== undefined) return n;
    }
    return undefined;
  }
}

/** Qualifiers/attributes that are copied to annotations (large or redundant ones are dropped). */
export const DROPPED_ATTRIBUTES = new Set(["translation"]);

export const RNA_TYPES = new Set(["mRNA", "ncRNA", "rRNA", "tRNA", "misc_RNA", "precursor_RNA", "tmRNA", "lnc_RNA", "primary_transcript", "transcript"]);

/** Bounding interval of a location on its outer sequence. */
export function extent(loc: Location): { ref: string; start: number; end: number } {
  const own = loc.segments.filter((s) => s.ref === loc.outer);
  const segs = own.length ? own : loc.segments;
  return { ref: loc.outer, start: Math.min(...segs.map((s) => s.start)), end: Math.max(...segs.map((s) => s.end)) };
}

/**
 * Completes version-less accessions (Ensembl GFF3 writes `protein_id=ENSP00000334393` without the version that its
 * protein FASTA carries) from the sequence keys known to the sources. Unknown or ambiguous keys are returned unchanged.
 */
export class VersionResolver {
  readonly #latest = new Map<string, string | null>();

  add(ref: string): void {
    const m = /^(.*)\.(\d+)$/.exec(ref);
    if (!m) return;
    const bare = m[1]!;
    const seen = this.#latest.get(bare);
    this.#latest.set(bare, seen === undefined ? ref : seen === ref ? ref : null);
  }

  resolve(ref: string): string {
    return this.#latest.get(ref) ?? ref;
  }
}

/**
 * Sequence-name -> RefSeq key from an NCBI assembly report (`*_assembly_report.txt`): the sequence name (`1`, `MT`),
 * GenBank accession (`CM000663.2`, `KI270728.1`), UCSC name (`chr1`) and RefSeq accession all map to `refseq:<RefSeq-Accn>`.
 */
export function assemblyReportSeqids(text: string): Map<string, string> {
  const out = new Map<string, string>();
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith("#")) continue;
    const cols = line.split("\t");
    const [name, , , , genbank, , refseq, , , ucsc] = cols;
    if (!refseq || refseq === "na") continue;
    const key = `refseq:${refseq}`;
    for (const alias of [name, genbank, refseq, ucsc]) if (alias && alias !== "na") out.set(alias, key);
  }
  return out;
}
