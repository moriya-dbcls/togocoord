// A set of stores (typically one per assembly/species) queried together.
import { createContext, NamespaceRegistry, type CoordContext, type Unit } from "@togocoord/core";
import { Lru, TogoCoordStore, type StoreOptions, type StoredBlock, type StoredEdge } from "@togocoord/ingest";
import { categoryOf, type Category } from "./category.ts";

/** A block with a globally unique edge key `<store>:<edge id>`. */
export interface SetBlock extends StoredBlock {
  key: string;
}

export class StoreSet {
  readonly stores: TogoCoordStore[] = [];
  readonly registry: NamespaceRegistry;
  /** Merged sequence records and identity lists are looked up for every node of every search: keep them. */
  readonly #sequences = new Lru<string, ReturnType<TogoCoordStore["sequence"]> | null>(100_000);
  readonly #identical = new Lru<string, string[]>(100_000);
  #meta: Array<Record<string, unknown>> | undefined;

  constructor(paths: string[] = [], options: StoreOptions & { registry?: NamespaceRegistry } = {}) {
    this.registry = options.registry ?? new NamespaceRegistry();
    for (const p of paths) this.stores.push(new TogoCoordStore(p, options));
  }

  add(store: TogoCoordStore): this {
    this.stores.push(store);
    this.#sequences.clear();
    this.#identical.clear();
    this.#meta = undefined;
    return this;
  }

  unitOf(ref: string): Unit | undefined {
    for (const s of this.stores) {
      const u = s.unitOf(ref);
      if (u) return u;
    }
    return undefined;
  }

  /** Units and lengths from the stores (a recorded length of 0 means unknown), then namespace defaults. */
  context(): CoordContext {
    return createContext({
      registry: this.registry,
      units: (ref) => this.unitOf(ref),
      lengths: (ref) => {
        const n = this.sequence(ref)?.length;
        return n && n > 0 ? n : undefined;
      },
    });
  }

  /**
   * The sequence record merged over all stores (first value per field). A store may only know a sequence as the end
   * of an edge (e.g. a protein in a genome GFF3 without its residues), another may hold its length and digest.
   */
  sequence(ref: string): ReturnType<TogoCoordStore["sequence"]> {
    const cached = this.#sequences.get(ref);
    if (cached !== undefined) return cached ?? undefined;
    let merged: ReturnType<TogoCoordStore["sequence"]>;
    for (const s of this.stores) {
      const r = s.sequence(ref);
      if (!r) continue;
      if (!merged) {
        merged = { ...r };
        continue;
      }
      // First known value per field; a length of 0 means unknown (e.g. MANE records); tags are united.
      const m = merged as Record<string, unknown>;
      for (const [k, v] of Object.entries(r)) {
        if (k === "tags") continue;
        if (m[k] === undefined || m[k] === null || (k === "length" && m[k] === 0)) m[k] = v;
      }
      if (r.tags?.length) merged.tags = [...new Set([...(merged.tags ?? []), ...r.tags])];
    }
    this.#sequences.set(ref, merged ?? null);
    return merged;
  }

  tags(ref: string): string[] {
    return this.sequence(ref)?.tags ?? [];
  }

  category(ref: string): Category {
    const seq = this.sequence(ref);
    const unit = seq?.unit ?? this.unitOf(ref) ?? this.registry.defaultUnit(ref);
    return categoryOf(ref, { ...(unit && { unit }), ...(seq?.moltype && { moltype: seq.moltype }) });
  }

  /** Blocks touching `[start, end)` of `ref` in every store, oriented away from `ref`. */
  blocksAt(ref: string, start: number, end: number): SetBlock[] {
    const out: SetBlock[] = [];
    this.stores.forEach((s, i) => {
      for (const b of s.blocksAt(ref, start, end)) out.push({ ...b, key: `${i}:${b.edge}` });
    });
    return out;
  }

  /** Other sequences with identical residues (same refget digest), in any store. */
  identical(ref: string): string[] {
    let out = this.#identical.get(ref);
    if (out === undefined) {
      const digest = this.sequence(ref)?.digest;
      out = digest ? this.refsByDigest(digest).filter((r) => r !== ref) : [];
      this.#identical.set(ref, out);
    }
    return out;
  }

  /** Sequences with this refget digest in any store. */
  refsByDigest(digest: string): string[] {
    const out = new Set<string>();
    for (const s of this.stores) for (const r of s.refsByDigest(digest)) out.add(r);
    return [...out];
  }

  /** Edges from or to `ref` in every store, keyed like blocks (`<store>:<edge id>`). */
  edges(ref: string): Array<StoredEdge & { key: string }> {
    return this.stores.flatMap((s, i) => s.edges(ref).map((e) => ({ ...e, key: `${i}:${e.id}` })));
  }

  /** Annotations overlapping `[start, end)` (internal units) of `ref` in every store. */
  annotations(ref: string, start: number, end: number): ReturnType<TogoCoordStore["annotations"]> {
    return this.stores.flatMap((s) => s.annotations(ref, start, end));
  }

  /** Per store: file name, recorded metadata (label, organism, assembly, inputs, ...) and content summary. */
  meta(): Array<Record<string, unknown>> {
    this.#meta ??= this.stores.map((s) => {
      const { summary: _recorded, ...meta } = s.meta();
      return { file: s.path.split(/[\\/]/).pop()!, ...meta, summary: s.summary() };
    });
    return this.#meta;
  }

  edge(key: string): StoredEdge | undefined {
    const [store, id] = key.split(":");
    return this.stores[Number(store)]?.edge(Number(id));
  }

  close(): void {
    for (const s of this.stores) s.close();
  }
}
