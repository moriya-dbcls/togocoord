// Location-aware path search over the edges of a StoreSet (design §8).
import { formatLocationId, Mapping, mapLocation, splitRef, type CoordContext, type Location } from "@togocoord/core";
import type { StoredEdge } from "@togocoord/ingest";
import type { Category } from "./category.ts";
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
  /** `<store>:<edge id>` */
  edge: string;
  kind: StoredEdge["kind"];
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
  provenance: StoredEdge["provenance"];
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
 * Extra cost of edges taken only as a last resort: those that failed self-validation, and those NCBI marks with
 * /exception unless validation against the actual sequence succeeded (e.g. a transcript that differs from the genome
 * only by substitutions keeps its coordinates).
 */
export const EXCEPTION_PENALTY = 10;
export const MISMATCH_PENALTY = 10;

export function edgeCost(e: Pick<StoredEdge, "attributes" | "validation">): number {
  return isApproximate(e) ? 1 + (e.validation.status === "mismatch" ? MISMATCH_PENALTY : EXCEPTION_PENALTY) : 1;
}

export function isApproximate(e: Pick<StoredEdge, "attributes" | "validation">): boolean {
  if (e.validation.status === "mismatch") return true;
  const verified = e.validation.status === "ok" && e.validation.basis === "full";
  return e.attributes.exception !== undefined && !verified;
}

interface State {
  location: Location;
  cost: number;
  path: Step[];
  orientation: Conversion["orientation"];
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

  const best = new Map<string, number>([[input.outer, 0]]);
  const queue: State[] = [{ location: input, cost: 0, path: [], orientation: "forward" }];
  const results: Conversion[] = [];

  while (queue.length) {
    const state = popCheapest(queue);
    const ref = state.location.outer;
    if (state.cost > (best.get(ref) ?? Infinity)) continue;
    if (state.path.length > 0 && matches(ref)) {
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

    for (const next of expand(stores, state, ctx)) {
      const r = next.location.outer;
      if (next.cost < (best.get(r) ?? Infinity)) {
        best.set(r, next.cost);
        queue.push(next);
      }
    }
  }
  return results;
}

function popCheapest(queue: State[]): State {
  let k = 0;
  for (let i = 1; i < queue.length; i++) if (queue[i]!.cost < queue[k]!.cost) k = i;
  return queue.splice(k, 1)[0]!;
}

/** One hop from `state` through every edge that overlaps its location. */
function expand(stores: StoreSet, state: State, ctx: CoordContext): State[] {
  const loc = state.location;
  const byEdge = new Map<string, SetBlock[]>();
  for (const seg of loc.segments) {
    const [a, b] = seg.start === seg.end ? [seg.start - 1, seg.start + 1] : [seg.start, seg.end];
    for (const blk of stores.blocksAt(seg.ref, a, b)) {
      const list = byEdge.get(blk.key);
      if (list) {
        if (!list.some((x) => x.src === blk.src && x.tgt === blk.tgt && x.len === blk.len && x.srcRef === blk.srcRef)) list.push(blk);
      } else byEdge.set(blk.key, [blk]);
    }
  }
  const used = new Set(state.path.map((s) => s.edge));
  const edges = [...byEdge.keys()].filter((k) => !used.has(k)).map((k) => ({ key: k, edge: stores.edge(k)! }));

  // Where an alignment connects the same two sequences (e.g. RefSeq transcript and genome via cDNA_match),
  // it describes the transcript's own sequence; the genome-model annotation edge does not.
  const aligned = new Set(edges.filter((x) => x.edge.kind === "alignment").map((x) => pairKey(x.edge)));
  const out: State[] = [];
  const inputId = formatLocationId(loc, ctx);
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
