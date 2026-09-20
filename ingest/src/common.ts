// Shared adapter plumbing.
import { createContext, NamespaceRegistry, type CoordContext, type Location, type Unit } from "@togocoord/core";
import type { SequenceRecord } from "./model.ts";
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
    const key = reportRef(refseq, genbank);
    if (!key) continue;
    for (const alias of [name, genbank, refseq, ucsc]) if (alias && alias !== "na") out.set(alias, key);
  }
  return out;
}

/**
 * UCSC database names of GRC assemblies (the assembly report does not carry them), so that input such as
 * `hg19:chr7:140453136` can name the assembly.
 */
export const UCSC_DATABASES: Record<string, string> = {
  GRCh38: "hg38",
  GRCh37: "hg19",
  GRCm39: "mm39",
  GRCm38: "mm10",
  "T2T-CHM13v2.0": "hs1",
};

/** Assembly name without its patch level (`GRCh37.p13` -> `GRCh37`). */
export function assemblyBaseName(name: string): string {
  return name.replace(/\.p\d+$/i, "");
}

/** Sequence key of an assembly report row: the RefSeq accession, else (INSDC-only assemblies) the GenBank one. */
function reportRef(refseq: string | undefined, genbank: string | undefined): string | undefined {
  if (refseq && refseq !== "na") return `refseq:${refseq}`;
  if (genbank && genbank !== "na") return `insdc:${genbank}`;
  return undefined;
}

/**
 * Sequence names of an assembly (Sequence-Name `7`, UCSC `chr7`, GenBank `CM000669.1`) -> sequence key
 * (`refseq:NC_000007.13`, or `insdc:AP031342.1` for an INSDC-only assembly), for input written with an assembly's own
 * names.
 */
export function assemblyReportAliases(text: string): Record<string, string> {
  const out: Record<string, string> = {};
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith("#")) continue;
    const [name, , , , genbank, , refseq, , , ucsc] = line.split("\t");
    const key = reportRef(refseq, genbank);
    if (!key) continue;
    for (const alias of [name, ucsc, genbank]) if (alias && alias !== "na" && out[alias] === undefined) out[alias] = key;
  }
  return out;
}

/**
 * Kind of molecule of each sequence of an assembly: `nuclear`, or the organelle (`mitochondrion`, `chloroplast`, ...),
 * from the Assigned-Molecule-Location/Type column. Unplaced scaffolds (na) are nuclear.
 */
export function assemblyReportMolecules(text: string): Map<string, string> {
  const out = new Map<string, string>();
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith("#")) continue;
    const [, , , type, genbank, , refseq] = line.split("\t");
    const ref = reportRef(refseq, genbank);
    if (!ref) continue;
    const t = (type ?? "").toLowerCase();
    out.set(ref, /mitochondri/.test(t) ? "mitochondrion" : /chloroplast|plastid|apicoplast|cyanelle/.test(t) ? "plastid" : "nuclear");
  }
  return out;
}

/** The sequences of an assembly as records (RefSeq or GenBank accession, length, species), from its NCBI assembly report. */
export function assemblyReportSequences(text: string, file?: string): SequenceRecord[] {
  const info = assemblyReportInfo(text);
  const out: SequenceRecord[] = [];
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith("#")) continue;
    const [name, , , type, genbank, , refseq, , length] = line.split("\t");
    const ref = reportRef(refseq, genbank);
    if (!ref || !length || !/^\d+$/.test(length)) continue;
    out.push({
      ref,
      moltype: "DNA",
      unit: "nt",
      length: Number(length),
      ...(/mitochondri|chloroplast|plastid/i.test(type ?? "") && { topology: "circular" as const }),
      ...(info.taxon && { taxon: Number(info.taxon) }),
      ...(info.organism && { organism: info.organism }),
      provenance: { adapter: "assembly-report", ...(file && { file }), record: name! },
    });
  }
  return out;
}

/**
 * Look up a sequence name, falling back to the GenBank accession inside UCSC names of patches and unplaced scaffolds
 * that assembly reports do not list (`chr11_GL456060_alt`, `chrUn_JH584304`, `chr1_KI270706v1_random` ->
 * GL456060.1, JH584304.1, KI270706.1).
 */
export function lookupSeqid(names: Map<string, string>, name: string): string | undefined {
  const direct = names.get(name);
  if (direct) return direct;
  const m = /^chr[0-9A-Za-z]+_([A-Z]{2}\d+)(?:v(\d+))?(?:_(?:alt|fix|random))?$/.exec(name);
  return m ? names.get(`${m[1]}.${m[2] ?? "1"}`) : undefined;
}

/** Header of an NCBI assembly report: organism, taxon, assembly name and accession, release date. */
export function assemblyReportInfo(text: string): { organism?: string; taxon?: string; assembly?: string; accession?: string; released?: string } {
  const field = (name: string) => {
    const v = new RegExp(`^# ${name}:\\s*(.+?)\\s*$`, "m").exec(text)?.[1];
    return v && v !== "n/a" ? v : undefined;
  };
  const out: { organism?: string; taxon?: string; assembly?: string; accession?: string; released?: string } = {};
  const organism = field("Organism name");
  const taxon = field("Taxid");
  const assembly = field("Assembly name");
  const accession = field("RefSeq assembly accession") ?? field("GenBank assembly accession");
  const released = field("Date");
  if (organism) out.organism = organism;
  if (taxon) out.taxon = taxon;
  if (assembly) out.assembly = assembly;
  if (accession) out.accession = accession.split(/\s/)[0]!;
  if (released && /^\d{4}-\d{2}-\d{2}$/.test(released)) out.released = released;
  return out;
}
