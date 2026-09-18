import assert from "node:assert/strict";
import { mkdtempSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { describe, it } from "node:test";
import { parseLocationId, residueBlock } from "@togocoord/core";
import {
  extract,
  ingestGenBank,
  ingestGff3,
  Lru,
  MemorySequenceSource,
  parseFasta,
  parseGenBank,
  qualifier,
  SqliteSink,
  TogoCoordStore,
  translate,
  type Edge,
  type IngestResult,
} from "@togocoord/ingest";
import { categoryOf, convert, StoreSet } from "../src/index.ts";

const dir = mkdtempSync(join(tmpdir(), "togocoord-service-"));
const fixture = (name: string) =>
  // Real records are shared with the ingest package.
  new URL(`../../ingest/test/fixtures/${name}`, import.meta.url);
const read = async (name: string) => (await import("node:fs")).readFileSync(fixture(name), "utf8");

function store(name: string, ...results: IngestResult[]): TogoCoordStore {
  const path = join(dir, `${name}.sqlite`);
  const sink = new SqliteSink(path);
  for (const r of results) {
    for (const s of r.sequences) sink.sequence(s);
    for (const e of r.edges) sink.edge(e);
    for (const a of r.annotations) sink.annotation(a);
  }
  sink.close();
  return new TogoCoordStore(path);
}

async function chr3Source(): Promise<MemorySequenceSource> {
  const source = new MemorySequenceSource();
  for (const [id, s] of parseFasta(await read("NM_012234.7.fa"))) source.add(`refseq:${id}`, s);
  for (const [, s] of parseFasta(await read("NC_000003.12_72374597-72446623.fa"))) source.add("refseq:NC_000003.12", s, 72374596);
  return source;
}

describe("protein -> RefSeq transcript -> GRCh38 through two stores (RYBP, real data)", async () => {
  const gb = await read("NM_012234.7.gb");
  const source = await chr3Source();
  const stores = new StoreSet()
    .add(store("rybp_rna", ingestGenBank(gb)))
    .add(store("rybp_genome", ingestGff3(await read("NC_000003.12_cDNA_match.gff3"), { source })));
  const ctx = stores.context();
  const cds = parseGenBank(gb)[0]!.features.find((f) => f.key === "CDS")!;
  const translation = qualifier(cds, "translation")!;

  it("maps every residue whose codon is aligned to the genome, and the codon encodes it", () => {
    assert.equal(translation.length, 228);
    const mapped: number[] = [];
    for (let r = 1; r <= translation.length; r++) {
      const results = convert(stores, parseLocationId(`refseq:NP_036366.3:${r}`, ctx), { to: { category: "genome" } }, ctx);
      if (results.length === 0) continue;
      const [hit] = results;
      assert.equal(results.length, 1);
      assert.deepEqual(hit!.path.map((s) => [s.kind, s.direction]), [["annotation", "forward"], ["alignment", "forward"]]);
      assert.equal(translate(extract(hit!.location, source)!), translation[r - 1], `residue ${r}`);
      mapped.push(r);
    }
    // NM_012234.7 1..213 is not aligned to GRCh38 (cDNA_match starts at 214); the CDS starts at 184.
    assert.deepEqual([mapped[0], mapped.at(-1), mapped.length], [11, 228, 218]);
  });

  it("reports what a step could not map and marks the truncation", () => {
    const [hit] = convert(stores, parseLocationId("refseq:NP_036366.3:9..12", ctx), { to: { category: "genome" } }, ctx);
    assert.equal(hit!.id, "refseq:NC_000003.12:complement(72446618..>72446623)");
    assert.equal(hit!.path[1]!.input, "refseq:NM_012234.7:208..219");
    assert.equal(hit!.path[1]!.unmapped, "refseq:NM_012234.7:208..213");
  });

  it("goes back from the genome to the protein", () => {
    const results = convert(stores, parseLocationId("refseq:NC_000003.12:complement(72446621..72446623)", ctx), { to: { ref: "refseq:NP_036366.3" } }, ctx);
    assert.deepEqual(results.map((r) => r.id), ["refseq:NP_036366.3:11"]);
    assert.deepEqual(results[0]!.path.map((s) => s.direction), ["inverse", "inverse"]);
  });
});

describe("path search within one genome (human mtDNA, real data)", async () => {
  const stores = new StoreSet().add(store("mt", ingestGenBank(await read("NC_012920.1.gb"))));
  const ctx = stores.context();
  const ids = (text: string, options: Parameters<typeof convert>[2]) => convert(stores, parseLocationId(text, ctx), options, ctx).map((r) => r.id);

  it("reaches an overlapping protein in another frame through the genome (ATP8 -> ATP6)", () => {
    const [hit] = convert(stores, parseLocationId("refseq:YP_003024030.1:55", ctx), { to: { ref: "refseq:YP_003024031.1" } }, ctx);
    assert.equal(hit!.id, "refseq:YP_003024031.1:1c2..2c1");
    assert.deepEqual(hit!.path.map((s) => `${s.from}>${s.to}:${s.direction}`), [
      "refseq:YP_003024030.1>refseq:NC_012920.1:forward",
      "refseq:YP_003024031.1>refseq:NC_012920.1:inverse",
    ]);
  });

  it("honours targets and hop limits", () => {
    assert.deepEqual(ids("refseq:YP_003024030.1:55", { to: { category: "genome" } }), ["refseq:NC_012920.1:8528..8530"]);
    assert.deepEqual(ids("refseq:YP_003024030.1:55", {}), ["refseq:NC_012920.1:8528..8530"]);
    assert.deepEqual(ids("refseq:YP_003024030.1:55", { to: { category: "protein" }, maxHops: 1 }), []);
    assert.deepEqual(ids("refseq:YP_003024030.1:55", { to: { category: "protein" } }), ["refseq:YP_003024031.1:1c2..2c1"]);
    assert.deepEqual(ids("refseq:NC_012920.1:8528", { to: { namespace: "refseq" } }).sort(), ["refseq:YP_003024030.1:55c1", "refseq:YP_003024031.1:1c2"]);
  });
});

describe("edge preferences", () => {
  const prov = { adapter: "gff3" as const };
  const ok = { status: "ok" as const };
  const edge = (e: Partial<Edge> & Pick<Edge, "from" | "to" | "blocks">): Edge => ({ kind: "annotation", attributes: {}, provenance: prov, validation: ok, ...e });

  it("uses the transcript alignment instead of the genome model of the same transcript", () => {
    const gff = [
      "##gff-version 3",
      "##sequence-region NC_000099.1 1 1000",
      "NC_000099.1\tt\tmRNA\t101\t200\t.\t+\t.\tID=rna-1;transcript_id=NM_000999.1",
      "NC_000099.1\tt\tcDNA_match\t101\t198\t.\t+\t.\tID=aln-1;Target=NM_000999.1 1 100 +;Gap=M50 I2 M48",
    ].join("\n");
    const stores = new StoreSet().add(store("pref", ingestGff3(gff)));
    const ctx = stores.context();
    const [hit] = convert(stores, parseLocationId("refseq:NM_000999.1:60", ctx), { to: { category: "genome" } }, ctx);
    assert.equal(hit!.id, "refseq:NC_000099.1:158");
    assert.equal(hit!.path[0]!.kind, "alignment");
  });

  it("takes an edge flagged with /exception only when nothing else maps", () => {
    const P = "refseq:NP_000001.1";
    const T = "refseq:NM_000001.1";
    const G = "refseq:NC_000001.1";
    const stores = new StoreSet().add(
      store("exception", {
        sequences: [],
        annotations: [],
        warnings: [],
        edges: [
          edge({ from: P, to: G, blocks: [{ srcRef: P, src: 0, tgtRef: G, tgt: 0, len: 30, rev: false }], attributes: { exception: "annotated by transcript or proteomic data" } }),
          edge({ from: P, to: T, blocks: [residueBlock({ srcRef: P, srcBegin: 1, tgtRef: T, tgtBegin: 1, length: 5 })] }),
          edge({ from: T, to: G, kind: "alignment", blocks: [{ srcRef: T, src: 0, tgtRef: G, tgt: 100, len: 15, rev: false }] }),
        ],
      }),
    );
    const ctx = stores.context();
    const first = convert(stores, parseLocationId(`${P}:1`, ctx), { to: { category: "genome" } }, ctx);
    assert.deepEqual(first.map((r) => [r.id, r.cost, r.path.length]), [[`${G}:101..103`, 2, 2]]);
    const eighth = convert(stores, parseLocationId(`${P}:8`, ctx), { to: { category: "genome" } }, ctx);
    assert.deepEqual(eighth.map((r) => [r.id, r.cost]), [[`${G}:22..24`, 11]]);
    // The same flag costs nothing once the edge was verified against the actual sequences.
    const verified = new StoreSet().add(
      store("exception-verified", {
        sequences: [],
        annotations: [],
        warnings: [],
        edges: [
          edge({
            from: P,
            to: G,
            blocks: [{ srcRef: P, src: 0, tgtRef: G, tgt: 0, len: 30, rev: false }],
            attributes: { exception: "annotated by transcript or proteomic data" },
            validation: { status: "ok", basis: "full", detail: "translation matches" },
          }),
        ],
      }),
    );
    assert.deepEqual(convert(verified, parseLocationId(`${P}:8`, ctx), { to: { category: "genome" } }, ctx).map((r) => r.cost), [1]);
    assert.equal(eighth[0]!.path[0]!.attributes.exception, "annotated by transcript or proteomic data");
    assert.deepEqual([first[0]!.approximate, eighth[0]!.approximate], [false, true]);
  });
});

describe("helpers", () => {
  it("categorises sequence keys", () => {
    assert.equal(categoryOf("refseq:NC_000001.11"), "genome");
    assert.equal(categoryOf("refseq:NM_012234.7"), "transcript");
    assert.equal(categoryOf("refseq:NG_044083.1"), "gene_region");
    assert.equal(categoryOf("refseq:NP_036366.3"), "protein");
    assert.equal(categoryOf("ensembl:ENST00000531224.5"), "transcript");
    assert.equal(categoryOf("insdc:AB000001.1", { moltype: "mRNA" }), "transcript");
    assert.equal(categoryOf("uniprot:P12883"), "protein");
  });
  it("bounds caches", () => {
    const lru = new Lru<string, number>(2).set("a", 1).set("b", 2);
    lru.get("a");
    lru.set("c", 3);
    assert.deepEqual([lru.has("a"), lru.has("b"), lru.has("c"), lru.size], [true, false, true, 2]);
  });
});
