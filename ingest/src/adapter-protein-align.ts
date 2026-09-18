// Protein alignment adapter (tier T2, spec-ingest §18): residue correspondences between a UniProt entry and the
// annotated proteins its ID mapping names (RefSeq, Ensembl, INSDC), for entries with no identical protein.
// An ID relation alone (UniProt P08556 <-> RefSeq NP_035067) does not say which residue is which; the alignment does.
// Candidates come from UniProt's idmapping_selected.tab (columns: 1 UniProtKB-AC, 4 RefSeq, 18 EMBL-CDS,
// 21 Ensembl_PRO); residues from the FASTA sources.
import { residueBlock, type Block } from "@togocoord/core";
import { ownString } from "./common.ts";
import type { Edge, Provenance, Sink } from "./model.ts";
import { refgetDigest, type SequenceSource } from "./sequence.ts";
import { readLines } from "./stream.ts";

// ---- alignment ----------------------------------------------------------------------------------------------------

export interface ProteinAlignment {
  /** Aligned runs (no gaps) as [start in a, start in b, length], 0-based residues. */
  blocks: Array<[number, number, number]>;
  identical: number;
  aligned: number;
  /** Aligned but different residues: [position in a, position in b (1-based), residue in a, residue in b]. */
  substitutions: Array<[number, number, string, string]>;
}

const MATCH = 2;
const MISMATCH = -1;
const GAP = -2;
/** Largest gap region filled by dynamic programming (cells); larger ones stay unaligned. */
const MAX_CELLS = 4_000_000;

/**
 * Align two similar proteins: k-mers occurring once in each are anchors; the longest colinear chain of anchors is kept
 * and the regions between anchors are filled by Needleman-Wunsch (free end gaps before the first and after the last
 * anchor). Runs of aligned columns less than half identical are dropped. Undefined when there is no anchor.
 */
export function alignProteins(a: string, b: string, k = 6): ProteinAlignment | undefined {
  const unique = (s: string) => {
    const pos = new Map<string, number>();
    for (let i = 0; i + k <= s.length; i++) {
      const w = s.slice(i, i + k);
      pos.set(w, pos.has(w) ? -1 : i);
    }
    return pos;
  };
  const ua = unique(a);
  const ub = unique(b);
  const pairs: Array<[number, number]> = [];
  for (const [w, i] of ua) {
    const j = ub.get(w);
    if (i >= 0 && j !== undefined && j >= 0) pairs.push([i, j]);
  }
  if (pairs.length === 0) return undefined;
  pairs.sort((x, y) => x[0] - y[0]);
  // Longest chain increasing in both (patience LIS on b positions).
  const tails: number[] = [];
  const prev = new Int32Array(pairs.length).fill(-1);
  for (let p = 0; p < pairs.length; p++) {
    let lo = 0;
    let hi = tails.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (pairs[tails[mid]!]![1] < pairs[p]![1]) lo = mid + 1;
      else hi = mid;
    }
    if (lo > 0) prev[p] = tails[lo - 1]!;
    tails[lo] = p;
  }
  const chain: Array<[number, number]> = [];
  for (let p = tails.at(-1)!; p >= 0; p = prev[p]!) chain.push(pairs[p]!);
  chain.reverse();

  // Columns: [i, j] aligned pairs, in order. A region filled by dynamic programming is kept when it looks aligned:
  // between anchors, equal lengths (substitutions only) or at least half identical; at an end (a different first or
  // last exon is common), at least half identical with 3 identical residues, so that chance pairs are not kept.
  const columns: Array<[number, number]> = [];
  let ai = 0;
  let bj = 0;
  const fill = (i1: number, j1: number, freeStart: boolean, freeEnd: boolean) => {
    const region = needlemanWunsch(a, b, ai, i1, bj, j1, freeStart, freeEnd);
    const same = region.filter(([i, j]) => a[i] === b[j]).length;
    const end = freeStart || freeEnd;
    const keep = end ? same * 2 >= region.length && same >= Math.min(3, region.length) : i1 - ai === j1 - bj || same * 2 >= region.length;
    if (keep) for (const c of region) columns.push(c);
    ai = i1;
    bj = j1;
  };
  for (const [i, j] of chain) {
    if (i < ai || j < bj) continue; // overlaps the previous anchor's k-mer
    fill(i, j, columns.length === 0, false);
    for (let d = 0; d < k; d++) columns.push([i + d, j + d]);
    ai = i + k;
    bj = j + k;
  }
  fill(a.length, b.length, columns.length === 0, true);

  // Runs of consecutive columns; drop runs less than half identical.
  const out: ProteinAlignment = { blocks: [], identical: 0, aligned: 0, substitutions: [] };
  let run: Array<[number, number]> = [];
  const flush = () => {
    if (run.length === 0) return;
    const same = run.filter(([i, j]) => a[i] === b[j]).length;
    if (same * 2 >= run.length) {
      out.blocks.push([run[0]![0], run[0]![1], run.length]);
      out.identical += same;
      out.aligned += run.length;
      for (const [i, j] of run) if (a[i] !== b[j]) out.substitutions.push([i + 1, j + 1, a[i]!, b[j]!]);
    }
    run = [];
  };
  for (const c of columns) {
    const last = run.at(-1);
    if (last && (c[0] !== last[0] + 1 || c[1] !== last[1] + 1)) flush();
    run.push(c);
  }
  flush();
  return out;
}

/** Aligned columns of a[i0, i1) and b[j0, j1); gaps are omitted. Regions over MAX_CELLS stay unaligned. */
function needlemanWunsch(a: string, b: string, i0: number, i1: number, j0: number, j1: number, freeStart: boolean, freeEnd: boolean): Array<[number, number]> {
  const n = i1 - i0;
  const m = j1 - j0;
  if (n === 0 || m === 0 || n * m > MAX_CELLS) return [];
  const W = m + 1;
  const S = new Int32Array((n + 1) * W);
  const T = new Uint8Array((n + 1) * W); // 0 diagonal, 1 up (gap in b), 2 left (gap in a)
  for (let i = 1; i <= n; i++) {
    S[i * W] = freeStart ? 0 : GAP * i;
    T[i * W] = 1;
  }
  for (let j = 1; j <= m; j++) {
    S[j] = freeStart ? 0 : GAP * j;
    T[j] = 2;
  }
  for (let i = 1; i <= n; i++) {
    const ca = a.charCodeAt(i0 + i - 1);
    for (let j = 1; j <= m; j++) {
      const d = S[(i - 1) * W + j - 1]! + (ca === b.charCodeAt(j0 + j - 1) ? MATCH : MISMATCH);
      const u = S[(i - 1) * W + j]! + GAP;
      const l = S[i * W + j - 1]! + GAP;
      const at = i * W + j;
      if (d >= u && d >= l) (S[at] = d), (T[at] = 0);
      else if (u >= l) (S[at] = u), (T[at] = 1);
      else (S[at] = l), (T[at] = 2);
    }
  }
  // With free end gaps, start the traceback from the best cell of the last row or column.
  let i = n;
  let j = m;
  if (freeEnd) {
    let best = S[n * W + m]!;
    for (let x = 0; x <= n; x++) if (S[x * W + m]! > best) (best = S[x * W + m]!), (i = x), (j = m);
    for (let y = 0; y <= m; y++) if (S[n * W + y]! > best) (best = S[n * W + y]!), (i = n), (j = y);
  }
  const out: Array<[number, number]> = [];
  while (i > 0 && j > 0) {
    const t = T[i * W + j];
    if (t === 0) {
      out.push([i0 + i - 1, j0 + j - 1]);
      i--;
      j--;
    } else if (t === 1) i--;
    else j--;
  }
  return out.reverse();
}

// ---- ID mapping ---------------------------------------------------------------------------------------------------

export interface IdMapping {
  /** UniProt accession -> the protein keys its row names (RefSeq, EMBL-CDS, Ensembl_PRO). */
  proteins: Map<string, string[]>;
  /** UniProt accession -> its genes (Entrez Gene ID, Ensembl gene ID). */
  genes: Map<string, string[]>;
  /** Gene -> the protein keys named by any entry of the gene. */
  geneProteins: Map<string, Set<string>>;
  /** UniProt accession -> its UniRef90 cluster (sequences at least 90% identical). */
  cluster: Map<string, string>;
  /** UniRef90 cluster -> the protein keys named by any of its entries. */
  clusterProteins: Map<string, Set<string>>;
}

/**
 * Candidates from idmapping_selected.tab: the proteins an entry's row names; else, through its genes (columns 3 GeneID
 * and 19 Ensembl), those named by the other entries of the gene; else those named by the other entries of its UniRef90
 * cluster (column 9). UniProt names a RefSeq protein on the entry with the identical sequence: NP_035067 on TrEMBL
 * A0A0G2JDN6, not on Swiss-Prot P08556, whose row has no gene either; both are in UniRef90_P08556.
 */
export async function readIdMapping(path: string): Promise<IdMapping> {
  const out: IdMapping = { proteins: new Map(), genes: new Map(), geneProteins: new Map(), cluster: new Map(), clusterProteins: new Map() };
  const split = (v: string | undefined) => (v ?? "").split(";").map((x) => x.trim()).filter((x) => x && x !== "-");
  for await (const line of readLines(path)) {
    if (!line) continue;
    const f = line.split("\t");
    const accession = ownString(f[0]!);
    const refs = [
      ...split(f[3]).map((x) => `refseq:${x}`),
      ...split(f[17]).map((x) => `insdc:${x}`),
      ...split(f[20]).map((x) => `ensembl:${x}`),
    ].map(ownString);
    const genes = [...split(f[2]).map((x) => `GeneID:${x}`), ...split(f[18]).map((x) => x.replace(/\.\d+$/, ""))].map(ownString);
    if (refs.length) out.proteins.set(accession, refs);
    if (genes.length) out.genes.set(accession, genes);
    for (const g of genes) {
      const set = out.geneProteins.get(g) ?? new Set<string>();
      for (const r of refs) set.add(r);
      out.geneProteins.set(g, set);
    }
    const cluster = f[8]?.trim();
    if (cluster) {
      out.cluster.set(accession, ownString(cluster));
      const set = out.clusterProteins.get(cluster) ?? new Set<string>();
      for (const r of refs) set.add(r);
      out.clusterProteins.set(ownString(cluster), set);
    }
  }
  return out;
}

// ---- adapter ------------------------------------------------------------------------------------------------------

export interface ProteinAlignOptions {
  file?: string;
  /** Residues of the UniProt entries and of the candidate proteins. */
  source: SequenceSource;
  /** Every sequence key the sources hold (UniProt entries and candidates). */
  refs: Iterable<string>;
  /** Completes version-less candidate keys (e.g. an Ensembl protein without its version). */
  resolveRef?: (ref: string) => string;
  /** Alignments below this identity (over aligned residues) are recorded as mismatches (default 0.9). */
  minIdentity?: number;
  /** ... or covering less than this share of the UniProt entry (default 0.5). */
  minCoverage?: number;
}

export interface ProteinAlignStats {
  entries: number;
  identical: number;
  noCandidate: number;
  /** Aligned to a protein of the same gene (no protein named on the entry's own row). */
  viaGene: number;
  /** Aligned to a protein of the same UniRef90 cluster (neither on the row nor through a gene). */
  viaCluster: number;
  aligned: number;
  ok: number;
  unaligned: number;
}

/**
 * One `alignment` edge (UniProt entry -> best candidate) per UniProt entry that has no identical protein among the
 * sources but has candidates in the ID mapping. Entries with an identical protein are reached by identity already.
 */
export async function ingestProteinAlignments(idmapping: string, sink: Sink, options: ProteinAlignOptions): Promise<ProteinAlignStats> {
  const provenance: Provenance = { adapter: "protein-alignment", ...(options.file && { file: options.file }) };
  const minIdentity = options.minIdentity ?? 0.9;
  const minCoverage = options.minCoverage ?? 0.5;
  const whole = (ref: string) => {
    const n = options.source.length?.(ref);
    return n === undefined ? undefined : options.source.get(ref, 0, n);
  };
  const all = [...options.refs];
  const entries = all.filter((r) => r.startsWith("uniprot:"));
  // Digests of the annotated proteins: a UniProt entry identical to one of them needs no alignment.
  const annotated = new Set<string>();
  for (const r of all) {
    if (r.startsWith("uniprot:")) continue;
    const s = whole(r);
    if (s) annotated.add(refgetDigest(s));
  }
  const mapping = await readIdMapping(idmapping);
  const stats: ProteinAlignStats = { entries: entries.length, identical: 0, noCandidate: 0, viaGene: 0, viaCluster: 0, aligned: 0, ok: 0, unaligned: 0 };
  for (const entry of entries) {
    const residues = whole(entry);
    if (!residues) continue;
    if (annotated.has(refgetDigest(residues))) {
      stats.identical++;
      continue;
    }
    const accession = entry.slice("uniprot:".length).replace(/-\d+$/, "");
    const usable = (refs: Iterable<string>) => [...new Set([...refs].map((r) => options.resolveRef?.(r) ?? r))].filter((r) => whole(r) !== undefined);
    let candidates = usable(mapping.proteins.get(accession) ?? []);
    let source = "UniProt ID mapping";
    if (candidates.length === 0) {
      candidates = usable((mapping.genes.get(accession) ?? []).flatMap((g) => [...(mapping.geneProteins.get(g) ?? [])]));
      source = "UniProt ID mapping, same gene";
    }
    const cluster = mapping.cluster.get(accession);
    if (candidates.length === 0 && cluster) {
      candidates = usable(mapping.clusterProteins.get(cluster) ?? []);
      source = `UniProt ID mapping, same ${cluster.split("_")[0]} cluster`;
    }
    if (candidates.length === 0) {
      stats.noCandidate++;
      continue;
    }
    let best: { ref: string; al: ProteinAlignment } | undefined;
    for (const c of candidates) {
      const al = alignProteins(residues, whole(c)!);
      if (al && (!best || al.identical > best.al.identical)) best = { ref: c, al };
    }
    if (!best || best.al.blocks.length === 0) {
      stats.unaligned++;
      continue;
    }
    const { al, ref } = best;
    const identity = al.identical / al.aligned;
    const coverage = al.aligned / residues.length;
    const ok = identity >= minIdentity && coverage >= minCoverage;
    const blocks: Block[] = al.blocks.map(([i, j, len]) => residueBlock({ srcRef: entry, srcBegin: i + 1, tgtRef: ref, tgtBegin: j + 1, length: len }));
    const edge: Edge = {
      kind: "alignment",
      from: entry,
      to: ref,
      blocks,
      attributes: {
        method: "anchored global protein alignment (TogoCoord)",
        candidates: source,
        identity: identity.toFixed(4),
        coverage: coverage.toFixed(4),
        // Aligned residues that differ, as from-position/to-position:from-residue>to-residue (at most 200).
        ...(al.substitutions.length && {
          substitutions: al.substitutions
            .slice(0, 200)
            .map(([i, j, x, y]) => `${i}/${j}:${x}>${y}`)
            .join(";"),
        }),
        ...(al.substitutions.length > 200 && { substitutionsTotal: String(al.substitutions.length) }),
      },
      provenance: { ...provenance, record: entry },
      validation: {
        status: ok ? "ok" : "mismatch",
        basis: "full",
        detail: `${al.identical}/${al.aligned} aligned residues identical; ${al.aligned}/${residues.length} of the entry aligned`,
      },
    };
    sink.edge(edge);
    stats.aligned++;
    if (source.endsWith("same gene")) stats.viaGene++;
    if (source.endsWith("cluster")) stats.viaCluster++;
    if (ok) stats.ok++;
  }
  return stats;
}
