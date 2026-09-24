// A set of stores (typically one per assembly/species) queried together.
import { createContext, NamespaceRegistry, type CoordContext, type Unit } from "@togocoord/core";
import { Lru, TogoCoordStore, type StoreOptions, type StoredBlock, type StoredEdge } from "@togocoord/ingest";
import { categoryOf, type Category } from "./category.ts";

export interface Species {
  taxon: number;
  organism?: string;
  /** Names as recorded, e.g. "Mus musculus" and "Mus musculus (house mouse)". */
  names: string[];
  /** Taxa folded into this species (subspecies, strains; see speciesTaxon). */
  taxa?: number[];
  assemblies: string[];
  /** Assembly of genome targets when none is requested and the input is not on a genome of the species. */
  defaultAssembly?: string;
}

export type StoredAnnotation = ReturnType<TogoCoordStore["annotations"]>[number] & { id?: string; link?: string };

export interface Assembly {
  name: string;
  taxon?: number;
  /** UCSC database name, e.g. hg19. */
  ucsc?: string;
  accession?: string;
  /** Sequence names (chr7, 7, CM000669.1) -> sequence key (refseq:NC_000007.13, insdc:AP031342.1). */
  aliases: Record<string, string>;
  /** Sequence keys of its sequences. */
  refs: Set<string>;
  /** Release date (YYYY-MM-DD, from the assembly report). */
  released?: string;
  /** Annotation edges in the stores of this assembly. */
  annotated: number;
}

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
  #species: Species[] | undefined;
  #assemblies: Assembly[] | undefined;
  #crossings: ReturnType<StoreSet["crossings"]> | undefined;
  #tagSpecies: ReturnType<StoreSet["tagSpecies"]> | undefined;
  #speciesOf: Map<number, number> | undefined;

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
    this.#assemblies = undefined;
    this.#crossings = undefined;
    this.#tagSpecies = undefined;
    this.#speciesOf = undefined;
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

  /**
   * Annotations overlapping `[start, end)` (internal units) of `ref` in every store. Those of a store built with
   * --id-namespace carry their ID in that namespace (`fanta:FCHS_1`) and, with --link, a URL.
   */
  annotations(ref: string, start: number, end: number): Array<StoredAnnotation> {
    const meta = this.meta();
    return this.stores.flatMap((s, i) => s.annotations(ref, start, end).map((a) => this.#withId(a, meta[i]!)));
  }

  /** An annotation by its ID in a namespace of annotation IDs (`fanta`, `FCHS_1`). */
  annotationById(namespace: string, id: string): StoredAnnotation | undefined {
    const meta = this.meta();
    for (const [i, s] of this.stores.entries()) {
      if (String(meta[i]!.id_namespace ?? "").toLowerCase() !== namespace.toLowerCase()) continue;
      const a = s.annotationById(id);
      if (a) return this.#withId(a, meta[i]!);
    }
    return undefined;
  }

  /** Namespaces of annotation IDs (e.g. fanta). */
  annotationNamespaces(): string[] {
    return [...new Set(this.meta().map((m) => m.id_namespace).filter((n): n is string => typeof n === "string"))];
  }

  #withId(a: ReturnType<TogoCoordStore["annotations"]>[number], meta: Record<string, unknown>): StoredAnnotation {
    const ns = meta.id_namespace as string | undefined;
    const id = a.attributes.ID?.[0];
    if (!ns || !id) return a;
    const link = typeof meta.link === "string" ? meta.link.replace("{id}", encodeURIComponent(id)) : undefined;
    return { ...a, id: `${ns}:${id}`, ...(link && { link }) };
  }

  /** Per store: file name, recorded metadata (label, organism, assembly, inputs, ...) and content summary. */
  meta(): Array<Record<string, unknown>> {
    this.#meta ??= this.stores.map((s) => {
      const { summary: _recorded, aliases: _aliases, ...meta } = s.meta();
      return { file: s.path.split(/[\\/]/).pop()!, ...meta, summary: s.summary() };
    });
    return this.#meta;
  }

  /**
   * Species a store's genome alignments convert from and to (a sample). A chain store holds the sequences of both
   * assemblies, so its content alone would put a mouse-to-human chain under human; between assemblies of one species
   * both are the same species. Not part of meta(), which the species lookup itself uses.
   */
  alignmentSpecies(store: number): { alignsFrom?: number; alignsTo?: number } {
    if (!this.#holds(store, "liftover")) return {};
    const ends = this.stores[store]?.edgeEnds("liftover", 1)[0];
    if (!ends) return {};
    const from = this.taxonOf(ends.from);
    const to = this.taxonOf(ends.to);
    return { ...(from !== undefined && { alignsFrom: from }), ...(to !== undefined && { alignsTo: to }) };
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
    if (own !== undefined && own !== null) return this.speciesTaxon(Number(own));
    const meta = this.meta();
    for (const [i, s] of this.stores.entries()) {
      const t = meta[i]!.taxon;
      if (t && s.sequence(ref)) return this.speciesTaxon(Number(t));
    }
    const taxa = new Set(
      this.identical(ref)
        .map((r) => this.sequence(r)?.taxon)
        .filter((t) => t !== undefined && t !== null)
        .map((t) => this.speciesTaxon(Number(t))),
    );
    return taxa.size === 1 ? [...taxa][0] : undefined;
  }

  /**
   * The species a taxon belongs to (spec-service §2.2): assemblies of subspecies, varieties and strains are assemblies
   * of their species, as are several individuals of one species. A taxon is folded into a species only when that is
   * explicit: the store says so (`--species-taxon`), or the NCBI name marks an infraspecific rank after the binomial
   * ("Marchantia polymorpha subsp. ruderalis" -> the loaded "Marchantia polymorpha"). Names alone are not enough
   * ("Human immunodeficiency virus 1" and "... 2" are two species).
   */
  speciesTaxon(taxon: number): number {
    if (!this.#speciesOf) {
      const map = new Map<number, number>();
      const names = new Map<number, Set<string>>();
      for (const m of this.meta()) {
        if (m.taxon && m.species_taxon) map.set(Number(m.taxon), Number(m.species_taxon));
        const entries: Array<{ taxon: number; organism?: string }> = [
          ...(m.taxon ? [{ taxon: Number(m.taxon), organism: m.organism as string | undefined }] : []),
          ...((m.summary as { taxa?: Array<{ taxon: number; organism?: string }> } | undefined)?.taxa ?? []),
        ];
        for (const e of entries) if (e.organism) (names.get(e.taxon) ?? names.set(e.taxon, new Set()).get(e.taxon)!).add(e.organism);
      }
      const plain = (n: string) => n.replace(/\s*\(.*\)\s*$/, "").trim();
      const byName = new Map<string, number>();
      for (const [t, ns] of names) for (const n of ns) if (/^[A-Z][a-z]+ [a-z][a-z-]+$/.test(plain(n))) byName.set(plain(n), t);
      for (const [t, ns] of names) {
        if (map.has(t)) continue;
        for (const n of ns) {
          const m = /^([A-Z][a-z]+ [a-z][a-z-]+) (?:subsp\.|var\.|f\.|str\.|strain|substr\.|serovar|biovar|pv\.|cv\.)/.exec(plain(n));
          const species = m ? byName.get(m[1]!) : undefined;
          if (species !== undefined && species !== t) map.set(t, species);
        }
      }
      this.#speciesOf = map;
    }
    return this.#speciesOf.get(taxon) ?? taxon;
  }

  #assembly(ref: string): string | undefined {
    if (this.category(ref) !== "genome") return undefined;
    const meta = this.meta();
    for (const [i, s] of this.stores.entries()) {
      const a = meta[i]!.assembly;
      if (typeof a === "string" && s.sequence(ref)) return a;
    }
    // A sequence the assembly names (its report), though no store of that assembly holds a record of it.
    return this.assemblies().find((x) => x.refs.has(ref))?.name;
  }

  /** Loaded species with their names and genome assemblies (from store metadata and content summaries). */
  species(): Species[] {
    if (this.#species) return this.#species;
    const out = new Map<number, Species>();
    // Names as recorded ("Mus musculus", "Mus musculus (house mouse)"); the shortest is shown.
    const add = (raw: number, organism?: string, assembly?: string) => {
      const taxon = this.speciesTaxon(raw);
      const e = out.get(taxon) ?? { taxon, names: [], assemblies: [] };
      if (raw !== taxon && !(e.taxa ??= []).includes(raw)) e.taxa.push(raw);
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
    // The default genome assembly of a species: an annotated one (where transcripts and proteins are), the newest
    // release first (MpTak_v7.1 over Marchanta_polymorpha_v1), then the most annotated (GRCh38 over GRCh37).
    for (const e of out.values()) {
      const [best] = this.assemblies()
        .filter((a) => a.taxon === e.taxon)
        .sort((a, b) => Number(b.annotated > 0) - Number(a.annotated > 0) || (b.released ?? "").localeCompare(a.released ?? "") || b.annotated - a.annotated);
      if (best) e.defaultAssembly = best.name;
    }
    this.#species = [...out.values()];
    return this.#species;
  }

  /** Genome assemblies of the stores (from --assembly-report): names, UCSC database name and sequence names. */
  assemblies(): Assembly[] {
    if (this.#assemblies) return this.#assemblies;
    const out = new Map<string, Assembly>();
    for (const s of this.stores) {
      const m = s.meta();
      if (!m.assembly) continue;
      const a = out.get(m.assembly) ?? { name: m.assembly, aliases: {}, refs: new Set<string>(), annotated: 0 };
      if (m.taxon) a.taxon = this.speciesTaxon(Number(m.taxon));
      if (m.released) a.released = m.released;
      if (m.ucsc) a.ucsc = m.ucsc;
      if (m.accession) a.accession = m.accession;
      if (m.aliases) {
        // Stores built before INSDC-only assemblies recorded bare RefSeq accessions.
        for (const [name, v] of Object.entries(JSON.parse(m.aliases) as Record<string, string>)) {
          const ref = v.includes(":") ? v : `refseq:${v}`;
          a.aliases[name] ??= ref;
          a.refs.add(ref);
        }
      }
      const edges = (s.summary() as { edges?: Record<string, number> }).edges ?? {};
      a.annotated += edges.annotation ?? 0;
      out.set(m.assembly, a);
    }
    this.#assemblies = [...out.values()];
    return this.#assemblies;
  }

  /** An assembly by name (`GRCh37.p13`, `GRCh37`, `grch37`), UCSC database name (`hg19`) or accession. */
  assembly(name: string): Assembly | undefined {
    const n = name.toLowerCase();
    return this.assemblies().find(
      (a) => a.name.toLowerCase() === n || a.name.toLowerCase().replace(/\.p\d+$/, "") === n || a.ucsc === n || a.accession?.toLowerCase() === n,
    );
  }

  /**
   * Whether a genome sequence belongs to an assembly: its store's assembly, or named in the assembly's report. A
   * sequence can belong to several (the mitochondrial NC_012920.1 is in GRCh37.p13 and GRCh38).
   */
  inAssembly(ref: string, name: string): boolean {
    if (this.assemblyOf(ref) === name) return true;
    return this.assemblies().find((a) => a.name === name)?.refs.has(ref) ?? false;
  }

  /**
   * What the liftOver chains connect: species and assemblies (from a sample of each chain store's edges). Used to
   * tell which species can be reached beyond identical sequences.
   */
  crossings(): Array<{ fromTaxon?: number; toTaxon?: number; fromAssembly?: string; toAssembly?: string }> {
    this.#crossings ??= this.stores.flatMap((s, i) => {
      if (!this.#holds(i, "liftover")) return [];
      const seen = new Set<string>();
      return s.edgeEnds("liftover").flatMap(({ from, to }) => {
        const c = {
          ...(this.taxonOf(from) !== undefined && { fromTaxon: this.taxonOf(from) }),
          ...(this.taxonOf(to) !== undefined && { toTaxon: this.taxonOf(to) }),
          ...(this.assemblyOf(from) && { fromAssembly: this.assemblyOf(from) }),
          ...(this.assemblyOf(to) && { toAssembly: this.assemblyOf(to) }),
        };
        const k = JSON.stringify(c);
        if (seen.has(k)) return [];
        seen.add(k);
        return [c];
      });
    });
    return this.#crossings;
  }

  /**
   * Whether a store holds edges of a kind, from the summary recorded at build time. Asking the store instead means a
   * query that scans its whole edge table when the kind is absent (`edge` has no index on `kind`): on the demo set
   * that cost ten seconds per `/v1/meta`, most of it in stores with no alignment at all.
   */
  #holds(store: number, kind: string): boolean {
    const edges = (this.meta()[store]?.summary as { edges?: Record<string, number> } | undefined)?.edges;
    return edges ? (edges[kind] ?? 0) > 0 : true;
  }

  /** Tags (e.g. MANE Select) and the species whose sequences carry them. */
  tagSpecies(): Array<{ tag: string; taxa: number[] }> {
    if (this.#tagSpecies) return this.#tagSpecies;
    const out = new Map<string, Set<number>>();
    for (const s of this.stores) {
      for (const [tag, refs] of s.tagged()) {
        const taxa = out.get(tag) ?? new Set<number>();
        for (const r of refs) {
          const t = this.taxonOf(r);
          if (t !== undefined) taxa.add(t);
        }
        out.set(tag, taxa);
      }
    }
    this.#tagSpecies = [...out].map(([tag, taxa]) => ({ tag, taxa: [...taxa] }));
    return this.#tagSpecies;
  }

  defaultAssembly(taxon: number | undefined): string | undefined {
    return taxon === undefined ? undefined : this.species().find((s) => s.taxon === taxon)?.defaultAssembly;
  }

  organismName(taxon: number): string | undefined {
    return this.species().find((s) => s.taxon === taxon)?.organism;
  }

  /** Mismatching aligned bases of a genome-alignment edge (`<store>:<edge id>`) within [start, end) of its source. */
  mismatches(key: string, ref: string, start: number, end: number): Array<{ pos: number; a: string; b: string }> {
    const [store, id] = key.split(":");
    return this.stores[Number(store)]?.mismatches(Number(id), ref, start, end) ?? [];
  }

  edge(key: string): StoredEdge | undefined {
    const [store, id] = key.split(":");
    return this.stores[Number(store)]?.edge(Number(id));
  }

  close(): void {
    for (const s of this.stores) s.close();
  }
}
