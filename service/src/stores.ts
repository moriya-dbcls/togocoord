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
  readonly #scope = new Lru<string, { taxon: number | null; assembly: string | null }>(100_000);
  #meta: Array<Record<string, unknown>> | undefined;
  #species: Array<{ taxon: number; organism?: string; names: string[]; assemblies: string[] }> | undefined;

  constructor(paths: string[] = [], options: StoreOptions & { registry?: NamespaceRegistry } = {}) {
    this.registry = options.registry ?? new NamespaceRegistry();
    for (const p of paths) this.stores.push(new TogoCoordStore(p, options));
  }

  add(store: TogoCoordStore): this {
    this.stores.push(store);
    this.#sequences.clear();
    this.#identical.clear();
    this.#scope.clear();
    this.#meta = undefined;
    this.#species = undefined;
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

  /**
   * Species of a sequence: its own record, else the recorded taxon of a store holding it (e.g. an Ensembl store built
   * with --taxon), else an identical sequence of a single species.
   */
  taxonOf(ref: string): number | undefined {
    return this.#scopeOf(ref).taxon ?? undefined;
  }

  /** Assembly of a genome sequence: the recorded assembly of the first store holding it (spec-service §2.2). */
  assemblyOf(ref: string): string | undefined {
    return this.#scopeOf(ref).assembly ?? undefined;
  }

  #scopeOf(ref: string): { taxon: number | null; assembly: string | null } {
    let out = this.#scope.get(ref);
    if (!out) {
      out = { taxon: this.#taxon(ref) ?? null, assembly: this.#assembly(ref) ?? null };
      this.#scope.set(ref, out);
    }
    return out;
  }

  #taxon(ref: string): number | undefined {
    const own = this.sequence(ref)?.taxon;
    if (own !== undefined && own !== null) return Number(own);
    const meta = this.meta();
    for (const [i, s] of this.stores.entries()) {
      const t = meta[i]!.taxon;
      if (t && s.sequence(ref)) return Number(t);
    }
    const taxa = new Set(this.identical(ref).map((r) => this.sequence(r)?.taxon).filter((t) => t !== undefined && t !== null));
    return taxa.size === 1 ? Number([...taxa][0]) : undefined;
  }

  #assembly(ref: string): string | undefined {
    if (this.category(ref) !== "genome") return undefined;
    const meta = this.meta();
    for (const [i, s] of this.stores.entries()) {
      const a = meta[i]!.assembly;
      if (typeof a === "string" && s.sequence(ref)) return a;
    }
    return undefined;
  }

  /** Loaded species with their names and genome assemblies (from store metadata and content summaries). */
  species(): Array<{ taxon: number; organism?: string; names: string[]; assemblies: string[] }> {
    if (this.#species) return this.#species;
    const out = new Map<number, { taxon: number; organism?: string; names: string[]; assemblies: string[] }>();
    // Names as recorded ("Mus musculus", "Mus musculus (house mouse)"); the shortest is shown.
    const add = (taxon: number, organism?: string, assembly?: string) => {
      const e = out.get(taxon) ?? { taxon, names: [], assemblies: [] };
      if (organism && !e.names.includes(organism)) e.names.push(organism);
      if (organism && (!e.organism || organism.length < e.organism.length)) e.organism = organism;
      if (assembly && !e.assemblies.includes(assembly)) e.assemblies.push(assembly);
      out.set(taxon, e);
    };
    for (const m of this.meta()) {
      if (m.taxon) add(Number(m.taxon), m.organism as string | undefined, m.assembly as string | undefined);
      const taxa = (m.summary as { taxa?: Array<{ taxon: number; organism?: string }> } | undefined)?.taxa ?? [];
      for (const t of taxa) add(t.taxon, t.organism);
    }
    this.#species = [...out.values()];
    return this.#species;
  }

  organismName(taxon: number): string | undefined {
    return this.species().find((s) => s.taxon === taxon)?.organism;
  }

  edge(key: string): StoredEdge | undefined {
    const [store, id] = key.split(":");
    return this.stores[Number(store)]?.edge(Number(id));
  }

  close(): void {
    for (const s of this.stores) s.close();
  }
}
