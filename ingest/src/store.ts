// SQLite store (node:sqlite) with R*Tree interval indexes (scaling §4).
import { existsSync, rmSync } from "node:fs";
import { DatabaseSync, type StatementSync } from "node:sqlite";
import {
  createContext,
  Mapping,
  mapLocation,
  NamespaceRegistry,
  type Block,
  type CoordContext,
  type Location,
  type MapResult,
  type Unit,
} from "@togocoord/core";
import { ownString } from "./common.ts";
import { Lru } from "./lru.ts";
import type { Annotation, Edge, Provenance, SequenceRecord, Sink, Validation } from "./model.ts";

export const STORE_SCHEMA_VERSION = "6";
/**
 * Versions the reader accepts (5 added chunked storage for directional edges, 6 the mismatching bases of genome
 * alignments; older stores simply have none).
 */
export const READABLE_SCHEMA_VERSIONS: readonly string[] = ["4", "5", "6"];

/** Blocks per chunk of a directional edge (spec-ingest §14). */
export const CHUNK_SIZE = 256;

const SCHEMA = `
CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE sequence(
  id INTEGER PRIMARY KEY, ref TEXT UNIQUE NOT NULL, moltype TEXT, unit TEXT, length INTEGER,
  topology TEXT, taxon INTEGER, organism TEXT, digest TEXT, md5 TEXT, tags TEXT, gene TEXT, provenance TEXT);
CREATE TABLE edge(
  id INTEGER PRIMARY KEY, kind TEXT NOT NULL, from_seq INTEGER NOT NULL, to_seq INTEGER NOT NULL,
  location TEXT, attributes TEXT, provenance TEXT, status TEXT, detail TEXT, basis TEXT);
CREATE TABLE block(
  id INTEGER PRIMARY KEY, edge INTEGER NOT NULL, src_seq INTEGER NOT NULL, src INTEGER NOT NULL,
  tgt_seq INTEGER NOT NULL, tgt INTEGER NOT NULL, len INTEGER NOT NULL, rev INTEGER NOT NULL);
CREATE TABLE annotation(
  id INTEGER PRIMARY KEY, seq INTEGER NOT NULL, start INTEGER NOT NULL, "end" INTEGER NOT NULL,
  type TEXT, location TEXT, attributes TEXT, provenance TEXT);
CREATE TABLE warning(id INTEGER PRIMARY KEY, message TEXT);
CREATE TABLE mismatch(edge INTEGER NOT NULL, seq INTEGER NOT NULL, pos INTEGER NOT NULL, a TEXT NOT NULL, b TEXT NOT NULL);
CREATE TABLE chunk(
  id INTEGER PRIMARY KEY, edge INTEGER NOT NULL, src_seq INTEGER NOT NULL, tgt_seq INTEGER NOT NULL,
  lo INTEGER NOT NULL, hi INTEGER NOT NULL, src0 INTEGER NOT NULL, tgt0 INTEGER NOT NULL, rev INTEGER NOT NULL,
  n INTEGER NOT NULL, data BLOB NOT NULL);
`;

/** Built after bulk loading, from rows sorted by position. Intervals are stored half-open as [lo, hi). */
const INDEXES = `
CREATE VIRTUAL TABLE block_src USING rtree_i32(id, seq_lo, seq_hi, lo, hi);
INSERT INTO block_src SELECT id, src_seq, src_seq, src, src + len FROM block ORDER BY src_seq, src;
CREATE VIRTUAL TABLE block_tgt USING rtree_i32(id, seq_lo, seq_hi, lo, hi);
INSERT INTO block_tgt SELECT id, tgt_seq, tgt_seq, tgt, tgt + len FROM block ORDER BY tgt_seq, tgt;
CREATE VIRTUAL TABLE chunk_src USING rtree_i32(id, seq_lo, seq_hi, lo, hi);
INSERT INTO chunk_src SELECT id, src_seq, src_seq, lo, hi FROM chunk ORDER BY src_seq, lo;
CREATE VIRTUAL TABLE annotation_idx USING rtree_i32(id, seq_lo, seq_hi, lo, hi);
INSERT INTO annotation_idx SELECT id, seq, seq, start, max("end", start + 1) FROM annotation ORDER BY seq, start;
CREATE INDEX mismatch_pos ON mismatch(seq, pos);
CREATE INDEX sequence_digest ON sequence(digest) WHERE digest IS NOT NULL;
CREATE INDEX sequence_md5 ON sequence(md5) WHERE md5 IS NOT NULL;
CREATE INDEX edge_from ON edge(from_seq);
CREATE INDEX edge_to ON edge(to_seq);
CREATE INDEX block_edge ON block(edge);
ANALYZE;
`;

export interface SqliteSinkOptions {
  /** Replace an existing file (default: refuse). */
  overwrite?: boolean;
  /** Statements per transaction. */
  batch?: number;
}

/** Bulk loader. Call `close()` to commit and build the indexes. */
export class SqliteSink implements Sink {
  readonly #db: DatabaseSync;
  readonly #seqIds = new Map<string, number>();
  readonly #batch: number;
  #pending = 0;
  readonly counts = { sequence: 0, edge: 0, block: 0, chunk: 0, annotation: 0, warning: 0 };
  readonly validation: Record<string, number> = {};
  readonly mismatches: string[] = [];
  readonly #s: Record<"newSeq" | "fillSeq" | "edge" | "block" | "chunk" | "mismatch" | "annotation" | "warning", StatementSync>;

  constructor(path: string, options: SqliteSinkOptions = {}) {
    if (existsSync(path)) {
      if (!options.overwrite) throw new Error(`${path} exists; pass overwrite to replace it`);
      rmSync(path);
    }
    this.#batch = options.batch ?? 50_000;
    this.#db = new DatabaseSync(path);
    this.#db.exec("PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA temp_store=FILE; PRAGMA cache_size=-65536;");
    this.#db.exec(SCHEMA);
    this.#s = {
      newSeq: this.#db.prepare("INSERT INTO sequence(ref) VALUES (?) RETURNING id"),
      fillSeq: this.#db.prepare(
        // Later records fill only missing fields (e.g. a FASTA record adds the checksums of a protein first seen in a GFF3).
        "UPDATE sequence SET moltype=COALESCE(moltype,?), unit=COALESCE(unit,?), length=COALESCE(length,?), topology=COALESCE(topology,?), " +
          "taxon=COALESCE(taxon,?), organism=COALESCE(organism,?), digest=COALESCE(digest,?), md5=COALESCE(md5,?), " +
          "tags=COALESCE(tags,?), gene=COALESCE(gene,?), provenance=COALESCE(provenance,?) WHERE id=?",
      ),
      edge: this.#db.prepare(
        "INSERT INTO edge(kind, from_seq, to_seq, location, attributes, provenance, status, detail, basis) VALUES (?,?,?,?,?,?,?,?,?) RETURNING id",
      ),
      block: this.#db.prepare("INSERT INTO block(edge, src_seq, src, tgt_seq, tgt, len, rev) VALUES (?,?,?,?,?,?,?)"),
      mismatch: this.#db.prepare("INSERT INTO mismatch(edge, seq, pos, a, b) VALUES (?,?,?,?,?)"),
      chunk: this.#db.prepare("INSERT INTO chunk(edge, src_seq, tgt_seq, lo, hi, src0, tgt0, rev, n, data) VALUES (?,?,?,?,?,?,?,?,?,?)"),
      annotation: this.#db.prepare('INSERT INTO annotation(seq, start, "end", type, location, attributes, provenance) VALUES (?,?,?,?,?,?,?)'),
      warning: this.#db.prepare("INSERT INTO warning(message) VALUES (?)"),
    };
    this.#db.exec("BEGIN");
  }

  #seq(ref: string): number {
    let id = this.#seqIds.get(ref);
    if (id === undefined) {
      id = Number((this.#s.newSeq.get(ref) as { id: number }).id);
      this.#seqIds.set(ownString(ref), id);
    }
    return id;
  }

  #tick(n = 1): void {
    this.#pending += n;
    if (this.#pending >= this.#batch) {
      this.#db.exec("COMMIT; BEGIN");
      this.#pending = 0;
    }
  }

  sequence(r: SequenceRecord): void {
    this.counts.sequence++;
    this.#s.fillSeq.run(
      r.moltype, r.unit, r.length, r.topology ?? null, r.taxon ?? null, r.organism ?? null, r.digest ?? null, r.md5 ?? null,
      r.tags?.length ? JSON.stringify(r.tags) : null, r.gene ?? null,
      JSON.stringify(r.provenance), this.#seq(r.ref),
    );
    this.#tick();
  }

  edge(e: Edge): void {
    this.counts.edge++;
    this.validation[e.validation.status] = (this.validation[e.validation.status] ?? 0) + 1;
    if (e.validation.status === "mismatch") this.mismatches.push(`${e.from} -> ${e.to}: ${e.validation.detail}`);
    const row = this.#s.edge.get(
      e.kind, this.#seq(e.from), this.#seq(e.to), e.location ?? null, JSON.stringify(e.attributes), JSON.stringify(e.provenance),
      e.validation.status, e.validation.detail ?? null, e.validation.basis ?? null,
    ) as { id: number };
    if (e.mismatches?.length) {
      const seq = this.#seq(e.from);
      for (const [pos, a, b] of e.mismatches) this.#s.mismatch.run(row.id, seq, pos, a, b);
      this.#tick(e.mismatches.length);
    }
    if (e.directional) {
      // Millions of blocks (whole-genome chains): compact chunks indexed by their source extent.
      const sorted = [...e.blocks].sort((a, b) => a.src - b.src);
      for (let i = 0; i < sorted.length; i += CHUNK_SIZE) {
        const part = sorted.slice(i, i + CHUNK_SIZE);
        const first = part[0]!;
        if (part.some((b) => b.srcRef !== first.srcRef || b.tgtRef !== first.tgtRef || b.rev !== first.rev)) {
          throw new Error(`directional edge ${e.from} -> ${e.to}: blocks must share sequences and orientation`);
        }
        const hi = Math.max(...part.map((b) => b.src + b.len));
        this.#s.chunk.run(row.id, this.#seq(first.srcRef), this.#seq(first.tgtRef), first.src, hi, first.src, first.tgt, first.rev ? 1 : 0, part.length, encodeChunk(part));
        this.counts.chunk++;
      }
      this.counts.block += e.blocks.length;
      this.#tick(1 + Math.ceil(e.blocks.length / CHUNK_SIZE));
      return;
    }
    for (const b of e.blocks) {
      this.#s.block.run(row.id, this.#seq(b.srcRef), b.src, this.#seq(b.tgtRef), b.tgt, b.len, b.rev ? 1 : 0);
    }
    this.counts.block += e.blocks.length;
    this.#tick(1 + e.blocks.length);
  }

  annotation(a: Annotation): void {
    this.counts.annotation++;
    this.#s.annotation.run(
      this.#seq(a.extent.ref), a.extent.start, a.extent.end, a.type, a.location, JSON.stringify(a.attributes), JSON.stringify(a.provenance),
    );
    this.#tick();
  }

  warning(message: string): void {
    this.counts.warning++;
    this.#s.warning.run(message);
    this.#tick();
  }

  /** Commit, build R*Tree and B-tree indexes, record metadata and a content summary. */
  close(meta: Record<string, string> = {}): void {
    this.#db.exec("COMMIT");
    this.#db.exec("BEGIN");
    this.#db.exec(INDEXES);
    // Annotations addressable by their ID in a namespace of their own (e.g. fanta:FCHS_1, spec-ingest §17).
    if (meta.id_namespace) {
      this.#db.exec(
        "CREATE TABLE annotation_id(name TEXT PRIMARY KEY, annotation INTEGER NOT NULL) WITHOUT ROWID;" +
          "INSERT OR IGNORE INTO annotation_id SELECT json_extract(attributes, '$.ID[0]'), id FROM annotation WHERE json_extract(attributes, '$.ID[0]') IS NOT NULL;",
      );
    }
    const put = this.#db.prepare("INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)");
    const summary = JSON.stringify(summarize(this.#db));
    for (const [k, v] of Object.entries({ schema: STORE_SCHEMA_VERSION, created: new Date().toISOString(), summary, ...meta })) put.run(k, v);
    this.#db.exec("COMMIT");
    this.#db.close();
  }
}

// ---- chunk encoding: per block zigzag LEB128 varints of (src - previous src, tgt - previous tgt, len) ----------------
function encodeChunk(blocks: ReadonlyArray<{ src: number; tgt: number; len: number }>): Uint8Array {
  const out: number[] = [];
  const put = (v: number) => {
    let z = v >= 0 ? v * 2 : -v * 2 - 1;
    while (z >= 0x80) {
      out.push((z % 0x80) | 0x80);
      z = Math.floor(z / 0x80);
    }
    out.push(z);
  };
  let src = blocks[0]!.src;
  let tgt = blocks[0]!.tgt;
  for (const b of blocks) {
    put(b.src - src);
    put(b.tgt - tgt);
    put(b.len);
    src = b.src;
    tgt = b.tgt;
  }
  return Uint8Array.from(out);
}

export function decodeChunk(data: Uint8Array, n: number, src0: number, tgt0: number): Array<{ src: number; tgt: number; len: number }> {
  let i = 0;
  const get = () => {
    let z = 0;
    let scale = 1;
    for (;;) {
      const byte = data[i++]!;
      z += (byte & 0x7f) * scale;
      if (byte < 0x80) break;
      scale *= 0x80;
    }
    return z % 2 === 0 ? z / 2 : -(z + 1) / 2;
  };
  const out: Array<{ src: number; tgt: number; len: number }> = [];
  let src = src0;
  let tgt = tgt0;
  for (let k = 0; k < n; k++) {
    src += get();
    tgt += get();
    out.push({ src, tgt, len: get() });
  }
  return out;
}

/** What a store holds: counts, organisms and example locations to try (shown by the service's /v1/meta). */
export interface StoreSummary {
  sequences: Record<string, number>;
  edges: Record<string, number>;
  blocks: number;
  annotations: number;
  taxa: Array<{ taxon: number; organism?: string; sequences: number }>;
  examples: string[];
}

export function summarize(db: DatabaseSync): StoreSummary {
  const all = <T>(sql: string) => db.prepare(sql).all() as T[];
  const one = (sql: string) => Number((db.prepare(sql).get() as { n: number }).n);
  const sequences = Object.fromEntries(
    all<{ moltype: string | null; n: number }>("SELECT moltype, count(*) AS n FROM sequence GROUP BY moltype").map((r) => [r.moltype ?? "unknown", Number(r.n)]),
  );
  const edges = Object.fromEntries(all<{ kind: string; n: number }>("SELECT kind, count(*) AS n FROM edge GROUP BY kind").map((r) => [r.kind, Number(r.n)]));
  const taxa = all<{ taxon: number; organism: string | null; n: number }>(
    "SELECT taxon, max(organism) AS organism, count(*) AS n FROM sequence WHERE taxon IS NOT NULL GROUP BY taxon ORDER BY n DESC",
  ).map((r) => ({ taxon: Number(r.taxon), ...(r.organism && { organism: r.organism }), sequences: Number(r.n) }));
  // Examples: a protein residue on a verified CDS edge, a residue on an alignment (e.g. structure), a whole protein.
  const examples: string[] = [];
  const cds = db
    .prepare(
      "SELECT f.ref FROM edge e JOIN sequence f ON f.id = e.from_seq WHERE e.kind = 'annotation' AND e.status = 'ok' AND f.moltype = 'protein' AND f.length > 60 ORDER BY e.id LIMIT 1 OFFSET 100",
    )
    .get() as { ref: string } | undefined;
  if (cds) examples.push(`${cds.ref}:50`, cds.ref);
  const aln = db
    .prepare("SELECT f.ref, b.src FROM edge e JOIN sequence f ON f.id = e.from_seq JOIN block b ON b.edge = e.id WHERE e.kind = 'alignment' AND e.status = 'ok' ORDER BY e.id LIMIT 1 OFFSET 100")
    .get() as { ref: string; src: number } | undefined;
  if (aln) {
    const aa = (db.prepare("SELECT unit FROM sequence WHERE ref = ?").get(aln.ref) as { unit: string | null } | undefined)?.unit === "aa";
    examples.push(`${aln.ref}:${aa ? Math.floor(aln.src / 3) + 1 : aln.src + 1}`);
  }
  const hasChunks = !!db.prepare("SELECT 1 FROM sqlite_master WHERE name = 'chunk'").get();
  if (hasChunks) {
    const c = db.prepare("SELECT s.ref, c.lo FROM chunk c JOIN sequence s ON s.id = c.src_seq ORDER BY c.id LIMIT 1 OFFSET 1000").get() as { ref: string; lo: number } | undefined;
    if (c) examples.push(`${c.ref}:${c.lo + 1}..${c.lo + 30}`);
  }
  if (!examples.length) {
    const seq = db.prepare("SELECT ref FROM sequence WHERE length > 0 ORDER BY id LIMIT 1").get() as { ref: string } | undefined;
    if (seq) examples.push(seq.ref);
  }
  const chunked = hasChunks ? one("SELECT coalesce(sum(n), 0) AS n FROM chunk") : 0;
  return { sequences, edges, blocks: one("SELECT count(*) AS n FROM block") + chunked, annotations: one("SELECT count(*) AS n FROM annotation"), taxa, examples };
}

export interface StoredEdge {
  id: number;
  kind: Edge["kind"];
  from: string;
  to: string;
  location?: string;
  attributes: Record<string, string>;
  provenance: Provenance;
  validation: Validation;
}

export interface StoredBlock extends Block {
  edge: number;
}

interface BlockRow {
  edge: number;
  src_seq: number;
  src: number;
  tgt_seq: number;
  tgt: number;
  len: number;
  rev: number;
}

export interface StoreOptions {
  /** Entries per lookup cache (ref <-> id, units, edges); least recently used entries are dropped. */
  cacheSize?: number;
  /** SQLite page cache in KiB (reads beyond it are served by the OS file cache). */
  pageCacheKiB?: number;
}

/** Read-only access for queries. Caches are bounded, so memory does not grow with the number of queries. */
export class TogoCoordStore {
  readonly path: string;
  readonly #db: DatabaseSync;
  readonly #refs: Lru<number, string>;
  readonly #ids: Lru<string, number | null>;
  readonly #units: Lru<string, Unit | null>;
  readonly #edges: Lru<number, StoredEdge>;
  readonly #s: Record<string, StatementSync>;

  constructor(path: string, options: StoreOptions = {}) {
    if (!existsSync(path)) throw new Error(`${path} does not exist`);
    this.path = path;
    this.#db = new DatabaseSync(path, { readOnly: true });
    const schema = (this.#db.prepare("SELECT value FROM meta WHERE key = 'schema'").get() as { value: string } | undefined)?.value;
    if (!schema || !READABLE_SCHEMA_VERSIONS.includes(schema)) {
      this.#db.close();
      throw new Error(`${path}: store schema ${schema ?? "?"} is not supported (expected ${STORE_SCHEMA_VERSION}); rebuild it with togocoord-ingest`);
    }
    this.#db.exec(`PRAGMA cache_size=-${options.pageCacheKiB ?? 8192};`);
    const n = options.cacheSize ?? 50_000;
    this.#refs = new Lru(n);
    this.#ids = new Lru(n);
    this.#units = new Lru(n);
    this.#edges = new Lru(n);
    this.#s = {
      edgeById: this.#db.prepare("SELECT * FROM edge WHERE id = ?"),
      id: this.#db.prepare("SELECT id FROM sequence WHERE ref = ?"),
      ref: this.#db.prepare("SELECT ref FROM sequence WHERE id = ?"),
      seq: this.#db.prepare("SELECT * FROM sequence WHERE ref = ?"),
      edgesOf: this.#db.prepare("SELECT * FROM edge WHERE from_seq = ? UNION ALL SELECT * FROM edge WHERE to_seq = ? AND from_seq <> to_seq"),
      bySrc: this.#db.prepare(
        "SELECT b.edge, b.src_seq, b.src, b.tgt_seq, b.tgt, b.len, b.rev FROM block_src r JOIN block b ON b.id = r.id WHERE r.seq_lo = ? AND r.seq_hi = ? AND r.lo < ? AND r.hi > ?",
      ),
      byTgt: this.#db.prepare(
        "SELECT b.edge, b.src_seq, b.src, b.tgt_seq, b.tgt, b.len, b.rev FROM block_tgt r JOIN block b ON b.id = r.id WHERE r.seq_lo = ? AND r.seq_hi = ? AND r.lo < ? AND r.hi > ?",
      ),
      annotations: this.#db.prepare(
        'SELECT a.* FROM annotation_idx r JOIN annotation a ON a.id = r.id WHERE r.seq_lo = ? AND r.seq_hi = ? AND r.lo < ? AND r.hi > ? ORDER BY a.start',
      ),
      meta: this.#db.prepare("SELECT key, value FROM meta"),
      ...(this.#db.prepare("SELECT 1 FROM sqlite_master WHERE name = 'chunk_src'").get() && {
        chunks: this.#db.prepare(
          "SELECT c.edge, c.src_seq, c.tgt_seq, c.src0, c.tgt0, c.rev, c.n, c.data FROM chunk_src r JOIN chunk c ON c.id = r.id WHERE r.seq_lo = ? AND r.seq_hi = ? AND r.lo < ? AND r.hi > ?",
        ),
      }),
      byDigest: this.#db.prepare("SELECT ref FROM sequence WHERE digest = ?"),
      byMd5: this.#db.prepare("SELECT ref FROM sequence WHERE md5 = ?"),
    };
  }

  #id(ref: string): number | undefined {
    let id = this.#ids.get(ref);
    if (id === undefined) {
      id = (this.#s.id!.get(ref) as { id: number } | undefined)?.id ?? null;
      this.#ids.set(ref, id);
    }
    return id ?? undefined;
  }

  #ref(id: number): string {
    let ref = this.#refs.get(id);
    if (ref === undefined) {
      ref = (this.#s.ref!.get(id) as { ref: string }).ref;
      this.#refs.set(id, ref);
    }
    return ref;
  }

  /** Content summary recorded at build time (computed now for stores built before summaries existed). */
  summary(): StoreSummary {
    const recorded = this.meta().summary;
    return recorded ? (JSON.parse(recorded) as StoreSummary) : summarize(this.#db);
  }

  meta(): Record<string, string> {
    return Object.fromEntries((this.#s.meta!.all() as Array<{ key: string; value: string }>).map((r) => [r.key, r.value]));
  }

  sequence(ref: string): (Partial<SequenceRecord> & { ref: string }) | undefined {
    const row = this.#s.seq!.get(ref) as Record<string, unknown> | undefined;
    if (!row) return undefined;
    const out: Partial<SequenceRecord> & { ref: string } = { ref };
    for (const k of ["moltype", "unit", "length", "topology", "taxon", "organism", "digest", "md5", "gene"] as const) {
      if (row[k] !== null) (out as Record<string, unknown>)[k] = row[k];
    }
    if (typeof row.tags === "string") out.tags = JSON.parse(row.tags);
    if (typeof row.provenance === "string") out.provenance = JSON.parse(row.provenance);
    return out;
  }

  unitOf(ref: string): Unit | undefined {
    let unit = this.#units.get(ref);
    if (unit === undefined) {
      unit = (this.sequence(ref)?.unit as Unit | undefined) ?? null;
      this.#units.set(ref, unit);
    }
    return unit ?? undefined;
  }

  /** Units from the store, then namespace defaults. */
  context(registry = new NamespaceRegistry()): CoordContext {
    return createContext({ registry, units: (ref) => this.unitOf(ref) });
  }

  /** Sequences with this refget digest (identical residues). */
  refsByDigest(digest: string): string[] {
    return (this.#s.byDigest!.all(digest) as Array<{ ref: string }>).map((r) => r.ref);
  }

  /** Sequences with this MD5 (hex, any case), e.g. to join UniParc. */
  refsByMd5(md5: string): string[] {
    return (this.#s.byMd5!.all(md5.toLowerCase()) as Array<{ ref: string }>).map((r) => r.ref);
  }

  /** Edges from or to `ref` (without blocks). */
  edges(ref: string): StoredEdge[] {
    const id = this.#id(ref);
    if (id === undefined) return [];
    return (this.#s.edgesOf!.all(id, id) as Array<Record<string, unknown>>).map((r) => this.#edgeRow(r));
  }

  /** Sequence pairs joined by edges of a kind (a sample, e.g. which assemblies liftOver chains connect). */
  edgeEnds(kind: Edge["kind"], limit = 20): Array<{ from: string; to: string }> {
    return this.#db
      .prepare(
        "SELECT f.ref AS f, t.ref AS t FROM edge e JOIN sequence f ON f.id = e.from_seq JOIN sequence t ON t.id = e.to_seq " +
          "WHERE e.kind = ? LIMIT ?",
      )
      .all(kind, limit)
      .map((r) => ({ from: String((r as { f: string }).f), to: String((r as { t: string }).t) }));
  }

  /** Tags in the store, each with a few sequences carrying it. */
  tagged(perTag = 5): Map<string, string[]> {
    const out = new Map<string, string[]>();
    for (const r of this.#db.prepare("SELECT ref, tags FROM sequence WHERE tags IS NOT NULL").iterate() as Iterable<{ ref: string; tags: string }>) {
      for (const tag of JSON.parse(r.tags) as string[]) {
        const list = out.get(tag) ?? [];
        if (list.length < perTag) list.push(r.ref);
        out.set(tag, list);
      }
    }
    return out;
  }

  /** One edge by id (cached). */
  edge(id: number): StoredEdge | undefined {
    const cached = this.#edges.get(id);
    if (cached) return cached;
    const row = this.#s.edgeById!.get(id) as Record<string, unknown> | undefined;
    return row && this.#edgeRow(row);
  }

  #edgeRow(r: Record<string, unknown>): StoredEdge {
    const e: StoredEdge = {
      id: Number(r.id),
      kind: r.kind as Edge["kind"],
      from: this.#ref(Number(r.from_seq)),
      to: this.#ref(Number(r.to_seq)),
      attributes: JSON.parse(String(r.attributes)),
      provenance: JSON.parse(String(r.provenance)),
      validation: {
        status: r.status as Validation["status"],
        ...(r.detail !== null && { detail: String(r.detail) }),
        ...(r.basis !== null && r.basis !== undefined && { basis: r.basis as NonNullable<Validation["basis"]> }),
      },
    };
    if (r.location !== null) e.location = String(r.location);
    this.#edges.set(e.id, e);
    return e;
  }

  /**
   * Blocks touching `[start, end)` of `ref`: `forward` blocks start on `ref`; `inverse` blocks end on `ref` and are
   * returned inverted, so both kinds map away from `ref`.
   */
  blocksAt(ref: string, start: number, end: number, direction: "forward" | "inverse" | "both" = "both"): StoredBlock[] {
    const id = this.#id(ref);
    if (id === undefined) return [];
    const out: StoredBlock[] = [];
    if (direction !== "inverse") {
      for (const r of this.#s.bySrc!.all(id, id, end, start) as unknown as BlockRow[]) {
        out.push({ edge: r.edge, srcRef: this.#ref(r.src_seq), src: r.src, tgtRef: this.#ref(r.tgt_seq), tgt: r.tgt, len: r.len, rev: r.rev === 1 });
      }
    }
    // Directional edges (chunked) are only used forward.
    if (direction !== "inverse" && this.#s.chunks) {
      type ChunkRow = { edge: number; src_seq: number; tgt_seq: number; src0: number; tgt0: number; rev: number; n: number; data: Uint8Array };
      for (const c of this.#s.chunks.all(id, id, end, start) as unknown as ChunkRow[]) {
        const srcRef = this.#ref(c.src_seq);
        const tgtRef = this.#ref(c.tgt_seq);
        for (const b of decodeChunk(c.data, c.n, c.src0, c.tgt0)) {
          if (b.src < end && b.src + b.len > start) out.push({ edge: c.edge, srcRef, src: b.src, tgtRef, tgt: b.tgt, len: b.len, rev: c.rev === 1 });
        }
      }
    }
    if (direction !== "forward") {
      for (const r of this.#s.byTgt!.all(id, id, end, start) as unknown as BlockRow[]) {
        out.push({ edge: r.edge, srcRef: this.#ref(r.tgt_seq), src: r.tgt, tgtRef: this.#ref(r.src_seq), tgt: r.src, len: r.len, rev: r.rev === 1 });
      }
    }
    return out;
  }

  /** Mapping restricted to the blocks relevant for `loc` (one hop from its sequences). */
  mappingFor(loc: Location, options: { direction?: "forward" | "inverse" | "both"; edges?: ReadonlySet<number> } = {}): Mapping {
    const seen = new Set<string>();
    const blocks: StoredBlock[] = [];
    for (const seg of loc.segments) {
      // A between-position needs both flanking units.
      const [a, b] = seg.start === seg.end ? [seg.start - 1, seg.start + 1] : [seg.start, seg.end];
      for (const blk of this.blocksAt(seg.ref, a, b, options.direction)) {
        if (options.edges && !options.edges.has(blk.edge)) continue;
        const key = `${blk.edge}|${blk.srcRef}|${blk.src}|${blk.tgtRef}|${blk.tgt}|${blk.len}|${blk.rev}`;
        if (seen.has(key)) continue;
        seen.add(key);
        blocks.push(blk);
      }
    }
    return new Mapping(blocks);
  }

  /** Every sequence directly connected to `loc`, with the converted locations. */
  neighbors(loc: Location, ctx: CoordContext = this.context()): MapResult {
    return mapLocation(loc, this.mappingFor(loc), ctx);
  }

  /** Mismatching aligned bases of an edge within [start, end) of its source sequence (schema 6). */
  mismatches(edge: number, ref: string, start: number, end: number): Array<{ pos: number; a: string; b: string }> {
    if (this.#hasMismatches === undefined) {
      this.#hasMismatches = this.#db.prepare("SELECT 1 FROM sqlite_master WHERE name = 'mismatch'").get() !== undefined;
    }
    const id = this.#hasMismatches ? this.#id(ref) : undefined;
    if (id === undefined) return [];
    return this.#db
      .prepare("SELECT pos, a, b FROM mismatch WHERE seq = ? AND pos >= ? AND pos < ? AND edge = ? ORDER BY pos")
      .all(id, start, end, edge)
      .map((r) => ({ pos: Number((r as { pos: number }).pos), a: String((r as { a: string }).a), b: String((r as { b: string }).b) }));
  }
  #hasMismatches: boolean | undefined;

  /** An annotation by its ID, in a store built with --id-namespace (e.g. the CRE FCHS_1). */
  annotationById(name: string): Omit<Annotation, "extent"> | undefined {
    if (this.#hasAnnotationIds === undefined) {
      this.#hasAnnotationIds = this.#db.prepare("SELECT 1 FROM sqlite_master WHERE name = 'annotation_id'").get() !== undefined;
    }
    if (!this.#hasAnnotationIds) return undefined;
    const r = this.#db.prepare("SELECT a.* FROM annotation_id i JOIN annotation a ON a.id = i.annotation WHERE i.name = ?").get(name) as
      | Record<string, unknown>
      | undefined;
    return r && {
      location: String(r.location),
      type: String(r.type),
      attributes: JSON.parse(String(r.attributes)),
      provenance: JSON.parse(String(r.provenance)),
    };
  }
  #hasAnnotationIds: boolean | undefined;

  annotations(ref: string, start: number, end: number): Array<Omit<Annotation, "extent">> {
    const id = this.#id(ref);
    if (id === undefined) return [];
    return (this.#s.annotations!.all(id, id, Math.max(end, start + 1), start) as Array<Record<string, unknown>>).map((r) => ({
      location: String(r.location),
      type: String(r.type),
      attributes: JSON.parse(String(r.attributes)),
      provenance: JSON.parse(String(r.provenance)),
    }));
  }

  counts(): Record<string, number> {
    const out: Record<string, number> = {};
    for (const t of ["sequence", "edge", "block", "annotation", "warning"]) {
      out[t] = Number((this.#db.prepare(`SELECT count(*) AS n FROM ${t}`).get() as { n: number }).n);
    }
    return out;
  }

  close(): void {
    this.#db.close();
  }
}
