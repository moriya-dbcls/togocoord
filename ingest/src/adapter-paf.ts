// PAF adapter (enrichment, tier T3): pairwise genome alignments from minimap2 and other aligners, computed by us or
// brought by users. Format: https://github.com/lh3/miniasm/blob/master/PAF.md
//   qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq [tags: cg:Z:<CIGAR> tp:A:<P|S|I> AS:i:<score>]
// The query is the assembly converted from and the target the one converted to (`minimap2 -c TO.fa FROM.fa`).
// Coordinates are 0-based half-open; on a '-' strand the CIGAR runs along the target forward and the query backward.
import type { Block } from "@togocoord/core";
import { NamespaceRegistry } from "@togocoord/core";
import { validateAlignedBlocks } from "./adapter-chain.ts";
import { ownString } from "./common.ts";
import type { Edge, Provenance, Sink } from "./model.ts";
import type { SequenceSource } from "./sequence.ts";
import { readLines } from "./stream.ts";

export interface PafRecord {
  qname: string;
  qlen: number;
  qstart: number;
  qend: number;
  strand: "+" | "-";
  tname: string;
  tlen: number;
  tstart: number;
  tend: number;
  nmatch: number;
  alen: number;
  mapq: number;
  /** P primary, S secondary, I inversion (minimap2 tp:A). */
  type?: string;
  score?: number;
  cigar?: string;
}

export function parsePafLine(line: string): PafRecord | undefined {
  const f = line.split("\t");
  if (f.length < 12) return undefined;
  const tag = (name: string) => f.slice(12).find((x) => x.startsWith(`${name}:`))?.split(":").slice(2).join(":");
  const score = tag("AS");
  return {
    qname: f[0]!,
    qlen: Number(f[1]),
    qstart: Number(f[2]),
    qend: Number(f[3]),
    strand: f[4] === "-" ? "-" : "+",
    tname: f[5]!,
    tlen: Number(f[6]),
    tstart: Number(f[7]),
    tend: Number(f[8]),
    nmatch: Number(f[9]),
    alen: Number(f[10]),
    mapq: Number(f[11]),
    ...(tag("tp") && { type: tag("tp") }),
    ...(score !== undefined && { score: Number(score) }),
    ...(tag("cg") && { cigar: tag("cg") }),
  };
}

/**
 * Ungapped blocks of an alignment from its CIGAR, from the query (`fromRef`) to the target (`toRef`), adjacent
 * matches merged. `=`/`X`/`M` consume both, `I` the query, `D`/`N` the target.
 */
export function pafBlocks(r: PafRecord, fromRef: string, toRef: string): Block[] {
  if (!r.cigar) throw new Error(`PAF record ${r.qname}:${r.qstart}-${r.qend} has no cg:Z CIGAR (run minimap2 with -c)`);
  const rev = r.strand === "-";
  const out: Array<{ -readonly [K in keyof Block]: Block[K] }> = [];
  let q = rev ? r.qend : r.qstart;
  let t = r.tstart;
  for (const [, n, op] of r.cigar.matchAll(/(\d+)([MIDNSHP=X])/g)) {
    const len = Number(n);
    if (op === "M" || op === "=" || op === "X") {
      const src = rev ? q - len : q;
      const last = out.at(-1);
      if (last && last.tgt + last.len === t && (rev ? last.src === src + len : last.src + last.len === src)) {
        if (rev) last.src = src;
        last.len += len;
      } else out.push({ srcRef: fromRef, src, tgtRef: toRef, tgt: t, len, rev });
      q += rev ? -len : len;
      t += len;
    } else if (op === "I") q += rev ? -len : len;
    else if (op === "D" || op === "N") t += len;
  }
  return out;
}

/** The part of a block whose source lies in [a, b) (a and b within the block). */
function clip(b: Block, a: number, e: number): Block {
  const offset = b.rev ? b.src + b.len - e : a - b.src;
  return { ...b, src: a, tgt: b.tgt + offset, len: e - a };
}

/** Parts of [start, end) not covered by `covered` (sorted, disjoint). */
function uncovered(covered: Array<[number, number]>, start: number, end: number): Array<[number, number]> {
  const out: Array<[number, number]> = [];
  let at = start;
  for (const [a, b] of covered) {
    if (b <= at) continue;
    if (a >= end) break;
    if (a > at) out.push([at, a]);
    at = Math.max(at, b);
  }
  if (at < end) out.push([at, end]);
  return out;
}

export interface PafOptions {
  file?: string;
  registry?: NamespaceRegistry;
  /** PAF sequence name of the "from" (query) assembly -> sequence key (see assemblyReportSeqids). */
  fromRef: (name: string) => string | undefined;
  /** PAF sequence name of the "to" (target) assembly -> sequence key. */
  toRef: (name: string) => string | undefined;
  source?: SequenceSource;
  /** Blocks sampled per alignment for validation (default 20). */
  sample?: number;
  minIdentity?: number;
  /** Alignments shorter than this on the query are ignored (default 1000). */
  minLength?: number;
  /** Source sequences that belong to both assemblies: identity, so their alignments are skipped (see ChainOptions). */
  shared?: (ref: string) => boolean;
  /** Whether two sequences may correspond (see ChainOptions.compatible). */
  compatible?: (from: string, to: string) => boolean;
}

export interface PafStats {
  records: number;
  alignments: number;
  blocks: number;
  skipped: number;
  /** Skipped because the source sequence belongs to both assemblies. */
  shared?: number;
  /** Skipped between a nuclear sequence and an organelle genome. */
  crossMolecule?: number;
  /** Query bases dropped because a better alignment already covered them (one-to-one on the source side). */
  overlapBases: number;
  sampledBases: number;
  identicalBases: number;
}

/**
 * One directional `liftover` edge per alignment, one-to-one on the source side like UCSC liftOver chains: alignments
 * are taken best first (AS, else matching bases) and each keeps only the source bases no better alignment covers.
 * Secondary alignments are ignored. The opposite direction comes from aligning the other way round.
 */
export async function ingestPafFile(path: string, sink: Sink, options: PafOptions): Promise<PafStats> {
  const provenance: Provenance = { adapter: "paf", ...(options.file && { file: options.file }) };
  const stats: PafStats = { records: 0, alignments: 0, blocks: 0, skipped: 0, overlapBases: 0, sampledBases: 0, identicalBases: 0 };
  const minLength = options.minLength ?? 1000;
  const records: PafRecord[] = [];
  const missing = new Set<string>();
  for await (const line of readLines(path)) {
    if (!line || line.startsWith("#")) continue;
    const r = parsePafLine(line);
    if (!r) continue;
    stats.records++;
    if (r.type === "S" || r.qend - r.qstart < minLength) {
      stats.skipped++;
      continue;
    }
    const names: Array<[string, string | undefined]> = [
      [`from ${r.qname}`, options.fromRef(r.qname)],
      [`to ${r.tname}`, options.toRef(r.tname)],
    ];
    if (names.some(([, ref]) => !ref)) {
      stats.skipped++;
      for (const [n, ref] of names) {
        if (!ref && !missing.has(n)) {
          missing.add(n);
          sink.warning(`${options.file ?? path}: no sequence for ${n}; its alignments are skipped`);
        }
      }
      continue;
    }
    if (options.compatible && !options.compatible(options.fromRef(r.qname)!, options.toRef(r.tname)!)) {
      stats.skipped++;
      stats.crossMolecule = (stats.crossMolecule ?? 0) + 1;
      continue;
    }
    if (options.shared?.(options.fromRef(r.qname)!)) {
      stats.skipped++;
      stats.shared = (stats.shared ?? 0) + 1;
      continue;
    }
    records.push({ ...r, qname: ownString(r.qname), tname: ownString(r.tname) });
  }
  records.sort((a, b) => (b.score ?? b.nmatch) - (a.score ?? a.nmatch) || b.nmatch - a.nmatch);
  const covered = new Map<string, Array<[number, number]>>();
  for (const r of records) {
    const cov = covered.get(r.qname) ?? [];
    const free = uncovered(cov, r.qstart, r.qend);
    stats.overlapBases += r.qend - r.qstart - free.reduce((n, [a, b]) => n + b - a, 0);
    if (free.length === 0) continue;
    const blocks = pafBlocks(r, options.fromRef(r.qname)!, options.toRef(r.tname)!).flatMap((b) =>
      free.filter(([a, e]) => a < b.src + b.len && e > b.src).map(([a, e]) => clip(b, Math.max(a, b.src), Math.min(e, b.src + b.len))),
    );
    for (const [a, e] of free) cov.push([a, e]);
    cov.sort((x, y) => x[0] - y[0]);
    covered.set(r.qname, cov);
    if (blocks.length === 0) continue;
    blocks.sort((x, y) => x.src - y.src);
    const edge: Edge = {
      kind: "liftover",
      directional: true,
      from: blocks[0]!.srcRef,
      to: blocks[0]!.tgtRef,
      blocks,
      attributes: {
        strand: r.strand,
        mapq: String(r.mapq),
        matches: String(r.nmatch),
        length: String(r.alen),
        ...(r.score !== undefined && { score: String(r.score) }),
        query: `${r.qname}:${r.qstart}-${r.qend}`,
        target: `${r.tname}:${r.tstart}-${r.tend}`,
      },
      provenance: { ...provenance, record: `${r.qname}:${r.qstart}-${r.qend}` },
      validation: validateAlignedBlocks(blocks, options.source, options.sample ?? 20, options.minIdentity ?? 0.5, stats),
    };
    sink.edge(edge);
    stats.alignments++;
    stats.blocks += blocks.length;
  }
  return stats;
}
