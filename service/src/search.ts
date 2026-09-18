// Location-aware path search over the edges of a StoreSet (design §8).
import { formatLocationId, Mapping, mapLocation, splitRef, type CoordContext, type Location } from "@togocoord/core";
import type { StoredEdge } from "@togocoord/ingest";
import { LAYER, type Category } from "./category.ts";
import type { SetBlock, StoreSet } from "./stores.ts";

export type Target = { ref: string } | { category: Category } | { namespace: string };

export interface ConvertOptions {
  /** Where to go; without it every sequence reached within `maxHops` is returned. */
  to?: Target | Target[];
  /** Default: 4 with a target, 1 without; plus CROSSING_HOPS when crossing species or assemblies. */
  maxHops?: number;
  /** Stop after this many results (default: unlimited). */
  maxResults?: number;
  /**
   * Tags of preferred sequences (e.g. "MANE Select"). Among paths of equal cost, those through preferred
   * intermediate sequences win; among results of equal cost, preferred targets come first. Costs are unchanged.
   */
  prefer?: string[];
  /**
   * Species of the targets (NCBI taxon). Default: the species of the input. Steps into another species (liftOver
   * chains, identical sequences of another species) are taken only when this names a species other than the input's
   * (spec-service §2.2).
   */
  taxon?: number;
  /** Assembly of genome targets (e.g. GRCm39). Default: the input's; another assembly is reached only when named. */
  assembly?: string;
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
  /** Tags of the target sequence (e.g. "MANE Select"). */
  tags: string[];
  /** Species of the target: tracked along the path (its own record may not carry one, e.g. Ensembl). */
  taxon?: number;
  cost: number;
  /**
   * The path uses an edge that failed self-validation or that NCBI flags with /exception without full verification:
   * positions may be shifted by indels between the sequences.
   */
  approximate: boolean;
  /**
   * Strand of the target: for nucleotides, the strand of its segments; for proteins (always written N->C), `reverse`
   * when the input corresponds to the antisense strand of the coding sequence.
   */
  orientation: "forward" | "reverse" | "mixed";
  path: Step[];
}

/**
 * Cost order (spec-service §3): identity (same digest) 0 < edge verified against or derived from primary data 1
 * < last resort 11. Exact paths therefore win over any path with an unverified or approximate step.
 */
export const IDENTITY_COST = 0;

/** Cross-assembly / cross-species liftOver: dearer than any edge within one assembly, so those paths win. */
export const LIFTOVER_COST = 2;

/**
 * Extra cost of edges taken only as a last resort: those that failed self-validation, and those NCBI marks with
 * /exception unless validation against the actual sequence succeeded (e.g. a transcript that differs from the genome
 * only by substitutions keeps its coordinates).
 */
export const EXCEPTION_PENALTY = 10;
export const MISMATCH_PENALTY = 10;

export function edgeCost(e: { kind?: StoredEdge["kind"]; attributes: Record<string, string>; validation: StoredEdge["validation"] }): number {
  const base = e.kind === "liftover" ? LIFTOVER_COST : 1;
  return isApproximate(e) ? base + (e.validation.status === "mismatch" ? MISMATCH_PENALTY : EXCEPTION_PENALTY) : base;
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
  /** Intermediate sequences without a preferred tag (secondary key after cost). */
  detours: number;
  /** Species and genome assembly the path is in (undefined until known). */
  taxon: number | undefined;
  assembly: string | undefined;
  /** The path has stepped into another species or assembly (at most once). */
  crossed: boolean;
  /** Insertion order, for deterministic tie-breaking. */
  seq: number;
}

/**
 * Extra hops allowed when crossing species or assemblies: the crossing step and the way back to the target's layer.
 * Without a target the search stops where it lands in the other species: the result is what the input corresponds to
 * there (spec-service §2.2).
 */
export const CROSSING_HOPS = 2;

/** Edge kinds whose ends are in different species or assemblies. */
const CROSSING_KINDS = new Set<Step["kind"]>(["liftover"]);

/**
 * Convert `input` to the requested targets along the cheapest chain of edges (Dijkstra over sequences).
 * Only edges whose blocks overlap the current location are followed, so the search stays local even in whole-genome
 * stores. Each sequence is reached once, by its cheapest path; a matching target is not expanded further.
 */
export function convert(stores: StoreSet, input: Location, options: ConvertOptions = {}, ctx: CoordContext = stores.context()): Conversion[] {
  const targets = options.to === undefined ? undefined : Array.isArray(options.to) ? options.to : [options.to];
  // Scope (spec-service §2.2): stay in the input's species and assembly unless another one is requested.
  const inputTaxon = stores.taxonOf(input.outer);
  const inputAssembly = stores.assemblyOf(input.outer);
  const crossing =
    (options.taxon !== undefined && options.taxon !== inputTaxon) || (options.assembly !== undefined && options.assembly !== inputAssembly);
  const maxHops = options.maxHops ?? (targets ? 4 : 1) + (crossing ? CROSSING_HOPS : 0);
  const inScope = (s: State) =>
    (options.taxon === undefined || s.taxon === options.taxon || (s.taxon === undefined && !crossing)) &&
    (options.assembly === undefined || stores.category(s.location.outer) !== "genome" || stores.assemblyOf(s.location.outer) === options.assembly);
  const matches = (ref: string) =>
    !targets ||
    targets.some((t) =>
      "ref" in t ? t.ref === ref : "category" in t ? stores.category(ref) === t.category : splitRef(ref).namespace === t.namespace,
    );

  const layerOf = (ref: string) => LAYER[stores.category(ref)];
  const cap = depthCap(stores, input.outer, targets, layerOf);
  // A sequence may be worth reaching in several layer states (trend, turns); keep the cheapest per state.
  const key = (s: Pick<State, "location" | "trend" | "turns" | "crossed">) => `${s.location.outer}|${s.trend}|${s.turns}|${s.crossed}`;
  const prefer = new Set(options.prefer ?? []);
  const preferred = (ref: string) => prefer.size > 0 && stores.tags(ref).some((t) => prefer.has(t));
  const start: State = {
    location: input,
    cost: 0,
    path: [],
    orientation: ctx.unitOf(input.outer) === "aa" ? "forward" : strandOrientation(input),
    trend: 0,
    turns: 0,
    detours: 0,
    taxon: inputTaxon,
    assembly: inputAssembly,
    crossed: false,
    seq: 0,
  };
  const best = new Map<string, [number, number]>([[key(start), [0, 0]]]);
  const better = (s: State, b: [number, number] | undefined) => !b || s.cost < b[0] || (s.cost === b[0] && s.detours < b[1]);
  const queue = new MinHeap<State>((a, b) => a.cost - b.cost || a.detours - b.detours || a.seq - b.seq);
  let seq = 1;
  queue.push(start);
  const results: Conversion[] = [];
  const reached = new Set<string>([input.outer]);

  while (queue.size) {
    const state = queue.pop()!;
    const ref = state.location.outer;
    const recorded = best.get(key(state));
    if (recorded && (state.cost > recorded[0] || (state.cost === recorded[0] && state.detours > recorded[1]))) continue;
    const target = state.path.length > 0 && matches(ref);
    if (target && inScope(state)) {
      if (reached.has(ref)) continue; // already returned through a cheaper state
      reached.add(ref);
      results.push({
        location: state.location,
        id: formatLocationId(state.location, ctx),
        category: stores.category(ref),
        tags: stores.tags(ref),
        ...(state.taxon !== undefined && { taxon: state.taxon }),
        cost: state.cost,
        approximate: state.path.some(isApproximate),
        orientation: state.orientation,
        path: state.path,
      });
      if (options.maxResults !== undefined && results.length >= options.maxResults) break;
      if (targets) continue; // a reached target is terminal
    }
    if (state.path.length >= maxHops) continue;
    if (!targets && state.crossed) continue; // landed in the other species (see CROSSING_HOPS)
    // A target outside the requested scope (e.g. the human protein when mouse is asked for) only leads on by crossing.
    const crossingOnly = target;

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
      const scope = nextScope(state, next);
      if (!scope || (crossingOnly && !(scope.crossed && !state.crossed))) continue;
      // The sequence being left becomes an intermediate node of the path (the source is not counted).
      const detours = state.detours + (state.path.length > 0 && prefer.size > 0 && !preferred(ref) ? 1 : 0);
      const candidate = { ...next, ...scope, trend, turns, detours, seq: seq++ };
      const k = key(candidate);
      if (better(candidate, best.get(k))) {
        best.set(k, [candidate.cost, candidate.detours]);
        queue.push(candidate);
      }
    }
  }
  // Where a step leads in species and assembly; undefined when the step leaves the requested scope.
  function nextScope(state: State, next: { location: Location; path: Step[] }): Pick<State, "taxon" | "assembly" | "crossed"> | undefined {
    const ref = next.location.outer;
    const taxon = stores.taxonOf(ref);
    const assembly = stores.assemblyOf(ref);
    const across =
      CROSSING_KINDS.has(next.path.at(-1)!.kind) ||
      (taxon !== undefined && state.taxon !== undefined && taxon !== state.taxon) ||
      (assembly !== undefined && state.assembly !== undefined && assembly !== state.assembly);
    if (!across) return { taxon: state.taxon ?? taxon, assembly: assembly ?? state.assembly, crossed: state.crossed };
    if (!crossing || state.crossed) return undefined;
    // Only into the requested species / assembly.
    const toTaxon = taxon ?? options.taxon ?? state.taxon;
    if (options.taxon !== undefined && toTaxon !== options.taxon) return undefined;
    if (options.taxon === undefined && toTaxon !== state.taxon) return undefined;
    if (options.assembly !== undefined && assembly !== undefined && assembly !== options.assembly) return undefined;
    return { taxon: toTaxon, assembly, crossed: true };
  }

  // Equal-cost results: preferred targets first (stable otherwise).
  if (prefer.size > 0) results.sort((a, b) => a.cost - b.cost || Number(preferred(b.location.outer)) - Number(preferred(a.location.outer)));
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
): Array<Omit<State, "seq" | "trend" | "turns" | "detours" | "taxon" | "assembly" | "crossed">> {
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
  const out: Array<Omit<State, "seq" | "trend" | "turns" | "detours" | "taxon" | "assembly" | "crossed">> = [];
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
        // Nucleotide segments carry their strand, so the step reports the true strand of its output; protein segments
        // are kept N->C (strand +), so a protein reached antisense keeps that in the state.
        orientation: ctx.unitOf(loc.outer) === "aa" ? combine(state.orientation, t.orientation) : t.orientation,
      });
    }
  }
  return out;
}

function pairKey(e: StoredEdge): string {
  return e.from < e.to ? `${e.from}\t${e.to}` : `${e.to}\t${e.from}`;
}

function strandOrientation(loc: Location): Conversion["orientation"] {
  const minus = loc.segments.filter((s) => s.strand === -1).length;
  return minus === 0 ? "forward" : minus === loc.segments.length ? "reverse" : "mixed";
}

function combine(a: Conversion["orientation"], b: Conversion["orientation"]): Conversion["orientation"] {
  if (a === "mixed" || b === "mixed") return "mixed";
  return a === b ? "forward" : "reverse";
}
