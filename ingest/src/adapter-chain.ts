// UCSC chain adapter (enrichment, tier T1/T3): whole-genome liftOver alignments between assemblies or species.
// Format: https://genome.ucsc.edu/goldenPath/help/chain.html
//   chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
//   size dt dq        (blocks separated by gaps dt on the target, dq on the query)
//   size              (last block)
// In liftOver files "t" is the assembly being lifted from and "q" the one lifted to. Coordinates are 0-based; on a
// '-' query strand they count from the end of the reverse complement.
import { NamespaceRegistry, type Block } from "@togocoord/core";
import { ownString } from "./common.ts";
import type { Edge, Provenance, Sink } from "./model.ts";
import { reverseComplement, type SequenceSource } from "./sequence.ts";
import { readLines } from "./stream.ts";

export interface ChainHeader {
  score: number;
  tName: string;
  tSize: number;
  tStart: number;
  tEnd: number;
  qName: string;
  qSize: number;
  qStrand: "+" | "-";
  qStart: number;
  qEnd: number;
  id: string;
}

export interface Chain extends ChainHeader {
  /** Ungapped blocks as [tStart, qStart (on the chain's query strand), size]. */
  blocks: Array<[number, number, number]>;
}

export async function* readChains(path: string): AsyncGenerator<Chain> {
  let chain: Chain | undefined;
  let t = 0;
  let q = 0;
  for await (const line of readLines(path)) {
    if (line.startsWith("chain")) {
      if (chain) yield chain;
      const f = line.trim().split(/\s+/);
      if (f.length < 13 || f[4] !== "+") throw new Error(`${path}: unsupported chain header: ${line}`);
      chain = {
        score: Number(f[1]),
        tName: f[2]!,
        tSize: Number(f[3]),
        tStart: Number(f[5]),
        tEnd: Number(f[6]),
        qName: f[7]!,
        qSize: Number(f[8]),
        qStrand: f[9] === "-" ? "-" : "+",
        qStart: Number(f[10]),
        qEnd: Number(f[11]),
        id: f[12]!,
        blocks: [],
      };
      t = chain.tStart;
      q = chain.qStart;
      continue;
    }
    if (!chain || line.trim() === "") continue;
    const [size, dt, dq] = line.trim().split(/\s+/).map(Number) as [number, number | undefined, number | undefined];
    chain.blocks.push([t, q, size]);
    t += size + (dt ?? 0);
    q += size + (dq ?? 0);
  }
  if (chain) yield chain;
}

/** Blocks of a chain in forward coordinates of both sequences (`rev` when the query is on the '-' strand). */
export function chainBlocks(c: Chain, fromRef: string, toRef: string): Block[] {
  const rev = c.qStrand === "-";
  return c.blocks.map(([t, q, size]) => ({
    srcRef: fromRef,
    src: t,
    tgtRef: toRef,
    tgt: rev ? c.qSize - q - size : q,
    len: size,
    rev,
  }));
}

export interface ChainOptions {
  file?: string;
  registry?: NamespaceRegistry;
  /** UCSC sequence name of the "from" (t) assembly -> sequence key (see assemblyReportSeqids). */
  fromRef: (name: string) => string | undefined;
  /** UCSC sequence name of the "to" (q) assembly -> sequence key. */
  toRef: (name: string) => string | undefined;
  /** Genome residues of both assemblies, for sampled identity checks. */
  source?: SequenceSource;
  /** Blocks sampled per chain for validation (default 20). */
  sample?: number;
  /** Identity below which a chain is reported as a mismatch (coordinates wrong); default 0.5. */
  minIdentity?: number;
}

export interface ChainStats {
  chains: number;
  blocks: number;
  skipped: number;
  sampledBases: number;
  identicalBases: number;
}

/**
 * One directional `liftover` edge per chain (from -> to). LiftOver chains are filtered on the "from" side only, so they
 * are not used in reverse; the opposite direction comes from the opposite file (e.g. mm39ToHg38).
 */
export async function ingestChainFile(path: string, sink: Sink, options: ChainOptions): Promise<ChainStats> {
  const provenance: Provenance = { adapter: "chain", ...(options.file && { file: options.file }) };
  const sample = options.sample ?? 20;
  const minIdentity = options.minIdentity ?? 0.5;
  const stats: ChainStats = { chains: 0, blocks: 0, skipped: 0, sampledBases: 0, identicalBases: 0 };
  const missing = new Set<string>();
  for await (const c of readChains(path)) {
    const from = options.fromRef(c.tName);
    const to = options.toRef(c.qName);
    if (!from || !to) {
      stats.skipped++;
      for (const n of [!from && `from ${c.tName}`, !to && `to ${c.qName}`]) {
        if (n && !missing.has(n)) {
          missing.add(n);
          sink.warning(`${options.file ?? path}: no sequence for ${n}; its chains are skipped`);
        }
      }
      continue;
    }
    const blocks = chainBlocks(c, ownString(from), ownString(to));
    const edge: Edge = {
      kind: "liftover",
      directional: true,
      from: blocks[0]!.srcRef,
      to: blocks[0]!.tgtRef,
      blocks,
      attributes: { chain: c.id, score: String(c.score), qStrand: c.qStrand },
      provenance: { ...provenance, record: `chain ${c.id}` },
      validation: validateAlignedBlocks(blocks, options.source, sample, minIdentity, stats),
    };
    sink.edge(edge);
    stats.chains++;
    stats.blocks += blocks.length;
  }
  return stats;
}

/**
 * Identity over up to `sample` evenly spaced blocks of a genome alignment (chain, PAF): aligned orthologous DNA
 * ~0.6-0.9, the same genome ~1.0, misplaced coordinates ~0.25.
 */
export function validateAlignedBlocks(
  blocks: Block[],
  source: SequenceSource | undefined,
  sample: number,
  minIdentity: number,
  stats: Pick<ChainStats, "sampledBases" | "identicalBases">,
): Edge["validation"] {
  if (!source) return { status: "skipped", detail: "genome sequences not available" };
  const step = Math.max(1, Math.floor(blocks.length / sample));
  let same = 0;
  let total = 0;
  for (let i = 0; i < blocks.length; i += step) {
    const b = blocks[i]!;
    const len = Math.min(b.len, 200);
    const s = source.get(b.srcRef, b.src, b.src + len);
    const t0 = b.rev ? source.get(b.tgtRef, b.tgt + b.len - len, b.tgt + b.len) : source.get(b.tgtRef, b.tgt, b.tgt + len);
    if (s === undefined || t0 === undefined) return { status: "skipped", detail: "genome sequences not available" };
    const t = b.rev ? reverseComplement(t0) : t0;
    for (let k = 0; k < len; k++) if (s[k] === t[k] && s[k] !== "N") same++;
    total += len;
  }
  stats.sampledBases += total;
  stats.identicalBases += same;
  return {
    status: total === 0 || same >= total * minIdentity ? "ok" : "mismatch",
    basis: "partial",
    detail: `${same}/${total} sampled aligned bases identical`,
  };
}
