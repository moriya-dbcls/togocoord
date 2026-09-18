// Location-aware path search over the edges of a StoreSet (design §8).
import { formatLocationId, Mapping, mapLocation, splitRef, type CoordContext, type Location } from "@togocoord/core";
import type { StoredEdge } from "@togocoord/ingest";
import { LAYER, type Category } from "./category.ts";
import type { SetBlock, StoreSet } from "./stores.ts";

export type Target = { ref: string } | { category: Category } | { namespace: string };

export interface ConvertOptions {
  /** Where to go; without it every sequence reached within `maxHops` is returned. */
  to?: Target | Target[];
  /** Default: 4 with a target, 1 without. */
  maxHops?: number;
  /** Stop after this many results (default: unlimited). */
  maxResults?: number;
}

export interface Step {
  /** `<store>:<edge id>`, or `digest:<refget digest>` for an identity step. */
  edge: string;
  /** `identity`: same residues (refget digest), same coordinates. */
  kind: StoredEdge["kind"] | "identity";
  from: string;
  to: string;
  /** `forward` when the step follows the edge from its `from` to its `to`. */
  direction: "forward" | "inverse";
  /** Canonical ID of the location entering this step. */
  input: string;
  /** Part of `input` this step could not map, or null. */
  unmapped: string | null;
  cost: number;
  edgeLocation?: string;
  attributes: Record<string, string>;
  validation: StoredEdge["validation"];
  /** Where the edge came from (absent for identity steps). */
  provenance?: StoredEdge["provenance"];
}

export interface Conversion {
  location: Location;
  /** Canonical Location ID of `location`. */
  id: string;
  category: Category;
  cost: number;
  /**
   * The path uses an edge that failed self-validation or that NCBI flags with /exception without full verification:
   * positions may be shifted by indels between the sequences.
   */
  approximate: boolean;
  orientation: "forward" | "reverse" | "mixed";
  path: Step[];
}

/**
 * Cost order (spec-service §3): identity (same digest) 0 < edge verified against or derived from primary data 1
 * < last resort 11. Exact paths therefore win over any path with an unverified or approximate step.
 */
export const IDENTITY_COST = 0;

/**
 * Extra cost of edges taken only as a last resort: those that failed self-validation, and those NCBI marks with
 * /exception unless validation against the actual sequence succeeded (e.g. a transcript that differs from the genome
 * only by substitutions keeps its coordinates).
 */
export const EXCEPTION_PENALTY = 10;
export const MISMATCH_PENALTY = 10;

export function edgeCost(e: { attributes: Record<string, string>; validation: StoredEdge["validation"] }): number {
  return isApproximate(e) ? 1 + (e.validation.status === "mismatch" ? MISMATCH_PENALTY : EXCEPTION_PENALTY) : 1;
}

export function isApproximate(e: { attributes: Record<string, string>; validation: StoredEdge["validation"] }): boolean {
  if (e.validation.status === "mismatch") return true;
  const verified = e.validation.status === "ok" && e.validation.basis === "full";
  return e.attributes.exception !== undefined && !verified;
}

interface State {
  location: Location;
  cost: number;
  path: Step[];
  orientation: Conversion["orientation"];
  /** Last non-lateral direction through the layers: -1 towards the genome, +1 towards structures, 0 none yet. */
  trend: -1 | 0 | 1;
  /** Direction changes so far (at most one U-turn is allowed). */
  turns: number;
  /** Insertion order, for deterministic tie-breaking. */
  seq: number;
}

/**
 * Convert `input` to the requested targets along the cheapest chain of edges (Dijkstra over sequences).
 * Only edges whose blocks overlap the current location are followed, so the search stays local even in whole-genome
 * stores. Each sequence is reached once, by its cheapest path; a matching target is not expanded further.
 */
export function convert(stores: StoreSet, input: Location, options: ConvertOptions = {}, ctx: CoordContext = stores.context()): Conversion[] {
  const targets = options.to === undefined ? undefined : Array.isArray(options.to) ? options.to : [options.to];
  const maxHops = options.maxHops ?? (targets ? 4 : 1);
  const matches = (ref: string) =>
    !targets ||
    targets.some((t) =>
      "ref" in t ? t.ref === ref : "category" in t ? stores.category(ref) === t.category : splitRef(ref).namespace === t.namespace,
    );

  const layerOf = (ref: string) => LAYER[stores.category(ref)];
  const cap = depthCap(stores, input.outer, targets, layerOf);
  // A sequence may be worth reaching in several layer states (trend, turns); keep the cheapest per state.
  const key = (s: Pick<State, "location" | "trend" | "turns">) => `${s.location.outer}|${s.trend}|${s.turns}`;
  const start: State = { location: input, cost: 0, path: [], orientation: "forward", trend: 0, turns: 0, seq: 0 };
  const best = new Map<string, number>([[key(start), 0]]);
  const queue = new MinHeap<State>((a, b) => a.cost - b.cost || a.seq - b.seq);
  let seq = 1;
  queue.push(start);
  const results: Conversion[] = [];
  const reached = new Set<string>([input.outer]);

  while (queue.size) {
    const state = queue.pop()!;
    const ref = state.location.outer;
    if (state.cost > (best.get(key(state)) ?? Infinity)) continue;
    if (state.path.length > 0 && matches(ref)) {
      if (reached.has(ref)) continue; // already returned through a cheaper state
      reached.add(ref);
      results.push({
        location: state.location,
        id: formatLocationId(state.location, ctx),
        category: stores.category(ref),
        cost: state.cost,
        approximate: state.path.some(isApproximate),
        orientation: state.orientation,
        path: state.path,
      });
      if (options.maxResults !== undefined && results.length >= options.maxResults) break;
      if (targets) continue; // a reached target is terminal
    }
    if (state.path.length >= maxHops) continue;

    const here = layerOf(ref);
    const withinCap = (r: string) => (layerOf(r) ?? -Infinity) <= cap;
    for (const next of expand(stores, state, ctx, withinCap)) {
      const there = layerOf(next.location.outer);
      if (there !== undefined && there > cap) continue; // deeper than both ends (spec-service §2.1, rule 1)
      const dir = here === undefined || there === undefined ? 0 : Math.sign(there - here);
      let { trend, turns } = state;
      if (dir !== 0) {
        if (trend !== 0 && dir !== trend) turns++;
        trend = dir as -1 | 1;
      }
      if (turns > MAX_TURNS) continue; // rule 2
      const candidate = { ...next, trend, turns, seq: seq++ };
      const k = key(candidate);
      if (candidate.cost < (best.get(k) ?? Infinity)) {
        best.set(k, candidate.cost);
        queue.push(candidate);
      }
    }
  }
  return results;
}

/** At most one U-turn through the layers (up towards the genome and back down, or the reverse). */
export const MAX_TURNS = 1;

/**
 * Deepest layer a path may enter: the deeper of the source and the targets (rule 1). Targets given only by namespace
 * have no known layer, so they do not bound the search.
 */
function depthCap(stores: StoreSet, source: string, targets: Target[] | undefined, layerOf: (ref: string) => number | undefined): number {
  if (!targets) return Infinity;
  let cap = layerOf(source) ?? Infinity;
  for (const t of targets) {
    const layer = "ref" in t ? layerOf(t.ref) : "category" in t ? LAYER[t.category] : undefined;
    cap = Math.max(cap, layer ?? Infinity);
  }
  return cap;
}

/** Binary min-heap; ties are broken by insertion order so results are deterministic. */
class MinHeap<T> {
  readonly #items: T[] = [];
  readonly #less: (a: T, b: T) => number;
  constructor(compare: (a: T, b: T) => number) {
    this.#less = compare;
  }
  get size(): number {
    return this.#items.length;
  }
  push(item: T): void {
    const a = this.#items;
    a.push(item);
    for (let i = a.length - 1; i > 0; ) {
      const p = (i - 1) >> 1;
      if (this.#less(a[i]!, a[p]!) >= 0) break;
      [a[i], a[p]] = [a[p]!, a[i]!];
      i = p;
    }
  }
  pop(): T | undefined {
    const a = this.#items;
    const top = a[0];
    const last = a.pop();
    if (a.length && last !== undefined) {
      a[0] = last;
      for (let i = 0; ; ) {
        const l = 2 * i + 1;
        const r = l + 1;
        let m = i;
        if (l < a.length && this.#less(a[l]!, a[m]!) < 0) m = l;
        if (r < a.length && this.#less(a[r]!, a[m]!) < 0) m = r;
        if (m === i) break;
        [a[i], a[m]] = [a[m]!, a[i]!];
        i = m;
      }
    }
    return top;
  }
}

/** One hop from `state` through every edge that overlaps its location (`seq` is assigned by the caller). */
function expand(
  stores: StoreSet,
  state: State,
  ctx: CoordContext,
  allowed: (ref: string) => boolean = () => true,
): Array<Omit<State, "seq" | "trend" | "turns">> {
  const loc = state.location;
  const byEdge = new Map<string, SetBlock[]>();
  for (const seg of loc.segments) {
    const [a, b] = seg.start === seg.end ? [seg.start - 1, seg.start + 1] : [seg.start, seg.end];
    for (const blk of stores.blocksAt(seg.ref, a, b)) {
      if (!allowed(blk.tgtRef)) continue; // blocks are oriented away from seg.ref
      const list = byEdge.get(blk.key);
      if (list) {
        if (!list.some((x) => x.src === blk.src && x.tgt === blk.tgt && x.len === blk.len && x.srcRef === blk.srcRef)) list.push(blk);
      } else byEdge.set(blk.key, [blk]);
    }
  }
  const out: Array<Omit<State, "seq" | "trend" | "turns">> = [];
  const inputId = formatLocationId(loc, ctx);

  // Identity: the same coordinates on every sequence with identical residues.
  if (loc.segments.every((s) => s.ref === loc.outer)) {
    const digest = stores.sequence(loc.outer)?.digest;
    for (const other of stores.identical(loc.outer)) {
      if (!allowed(other)) continue;
      const location: Location = { ...loc, outer: other, segments: loc.segments.map((s) => ({ ...s, ref: other })) };
      const step: Step = {
        edge: `digest:${digest}`,
        kind: "identity",
        from: loc.outer,
        to: other,
        direction: "forward",
        input: inputId,
        unmapped: null,
        cost: IDENTITY_COST,
        attributes: { digest: digest! },
        validation: { status: "ok", basis: "full", detail: "identical residues (refget digest)" },
      };
      out.push({ location, cost: state.cost + IDENTITY_COST, path: [...state.path, step], orientation: state.orientation });
    }
  }

  const used = new Set(state.path.map((s) => s.edge));
  const edges = [...byEdge.keys()].filter((k) => !used.has(k)).map((k) => ({ key: k, edge: stores.edge(k)! }));

  // Where an alignment connects the same two sequences (e.g. RefSeq transcript and genome via cDNA_match),
  // it describes the transcript's own sequence; the genome-model annotation edge does not.
  const aligned = new Set(edges.filter((x) => x.edge.kind === "alignment").map((x) => pairKey(x.edge)));
  for (const { key, edge } of edges) {
    if (edge.kind === "annotation" && aligned.has(pairKey(edge))) continue;
    const result = mapLocation(loc, new Mapping(byEdge.get(key)!), ctx);
    const cost = edgeCost(edge);
    for (const t of result.targets) {
      if (t.location.outer === loc.outer) continue;
      const step: Step = {
        edge: key,
        kind: edge.kind,
        from: edge.from,
        to: edge.to,
        direction: edge.from === loc.outer ? "forward" : "inverse",
        input: inputId,
        unmapped: result.unmapped && formatLocationId(result.unmapped, ctx),
        cost,
        attributes: edge.attributes,
        validation: edge.validation,
        provenance: edge.provenance,
      };
      if (edge.location !== undefined) step.edgeLocation = edge.location;
      out.push({
        location: t.location,
        cost: state.cost + cost,
        path: [...state.path, step],
        orientation: combine(state.orientation, t.orientation),
      });
    }
  }
  return out;
}

function pairKey(e: StoredEdge): string {
  return e.from < e.to ? `${e.from}\t${e.to}` : `${e.to}\t${e.from}`;
}

function combine(a: Conversion["orientation"], b: Conversion["orientation"]): Conversion["orientation"] {
  if (a === "mixed" || b === "mixed") return "mixed";
  return a === b ? "forward" : "reverse";
}
