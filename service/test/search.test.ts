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
import { ingestFastaFile, ingestManeSummary, MemorySink } from "@togocoord/ingest";
import { fileURLToPath } from "node:url";
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

describe("identity by refget digest (UniProt P07203 = RefSeq NP_000572.2, real data)", async () => {
  const uniprot = new MemorySink();
  await ingestFastaFile(fileURLToPath(fixture("uniprot_P07203.fa")), uniprot);
  const stores = new StoreSet().add(store("gpx1_rna", ingestGenBank(await read("NM_000581.4.gb")))).add(store("uniprot", uniprot.result));
  const ctx = stores.context();

  it("reaches the transcript through the identical RefSeq protein at no extra cost", () => {
    const [hit] = convert(stores, parseLocationId("uniprot:P07203:49", ctx), { to: { category: "transcript" } }, ctx);
    // Residue 49 is the selenocysteine, encoded by the UGA at NM_000581.4:220..222 (/transl_except).
    assert.equal(hit!.id, "refseq:NM_000581.4:220..222");
    assert.equal(hit!.cost, 1);
    assert.deepEqual(hit!.path.map((s) => [s.kind, s.to]), [["identity", "refseq:NP_000572.2"], ["annotation", "refseq:NM_000581.4"]]);
    assert.equal(hit!.approximate, false);
  });

  it("finds identical sequences in both directions", () => {
    assert.deepEqual(stores.identical("uniprot:P07203"), ["refseq:NP_000572.2"]);
    assert.deepEqual(stores.identical("refseq:NP_000572.2"), ["uniprot:P07203"]);
    const back = convert(stores, parseLocationId("refseq:NM_000581.4:220..222", ctx), { to: { namespace: "uniprot" } }, ctx);
    assert.deepEqual(back.map((r) => r.id), ["uniprot:P07203:49"]);
  });
});

describe("genome -> PDB structure: Ensembl CDS, identical UniProt protein, SIFTS (GPX1, real data)", async () => {
  const { assemblyReportSeqids, defaultFastaRef, ingestSiftsFile, parseFastaHeaders } = await import("@togocoord/ingest");
  const { NamespaceRegistry } = await import("@togocoord/core");
  const registry = new NamespaceRegistry();
  const residues = new MemorySequenceSource();
  for (const name of ["ensembl_GRCh38.116_GPX1.pep.fa", "uniprot_P07203.fa", "pdb_seqres_P07203.txt"]) {
    for (const [header, s] of parseFastaHeaders(await read(name))) residues.add(defaultFastaRef(header, registry)!, s);
  }
  for (const [, s] of parseFasta(await read("NC_000003.12_49357176-49358353.fa"))) residues.add("refseq:NC_000003.12", s, 49357175);
  const seqids = assemblyReportSeqids(await read("GRCh38.p14_assembly_report_chr3.txt"));
  const uniprot = new MemorySink();
  await ingestFastaFile(fileURLToPath(fixture("uniprot_P07203.fa")), uniprot);
  const mane = new MemorySink();
  await ingestManeSummary(fileURLToPath(fixture("MANE.GRCh38.v1.5.summary_GPX1_RYBP.txt")), mane);
  const sifts = new MemorySink();
  await ingestSiftsFile(fileURLToPath(fixture("sifts_P07203.tsv")), sifts, { source: residues });
  const stores = new StoreSet()
    .add(store("ensembl_gpx1", ingestGff3(await read("ensembl_GRCh38.116_GPX1.gff3"), { source: residues, seqidToRef: (id) => seqids.get(id) })))
    .add(store("uniprot_gpx1", uniprot.result))
    .add(store("sifts_gpx1", sifts.result))
    .add(store("mane_gpx1", mane.result));
  const ctx = stores.context();

  it("maps the selenocysteine codon to residue 59 of both 2F8A chains (the U49G mutant)", () => {
    const [codon] = convert(stores, parseLocationId("ensembl:ENSP00000407375.1:49", ctx), { to: { category: "genome" } }, ctx);
    // CDS starts at 49358278 on the minus strand: residue 49 = CDS nt 145..147 = 49358134..49358132.
    assert.equal(codon!.id, "refseq:NC_000003.12:complement(49358132..49358134)");
    assert.equal(translate(extract(codon!.location, residues)!), "*"); // UGA read as selenocysteine
    const hits = convert(stores, codon!.location, { to: { namespace: "pdb" } }, ctx);
    assert.deepEqual(hits.map((h) => h.id).sort(), ["pdb:2F8A.A:59", "pdb:2F8A.B:59"]);
    const path = hits[0]!.path;
    assert.deepEqual(path.map((s) => [s.kind, s.to]), [
      ["annotation", "refseq:NC_000003.12"],
      ["identity", "uniprot:P07203"],
      ["alignment", hits[0]!.location.outer],
    ]);
    // Several Ensembl proteins of GPX1 are identical to P07203; any of them gives the same, exact path.
    assert.match(path[0]!.from, /^ensembl:ENSP\d{11}\.\d+$/);
    assert.ok(stores.identical("uniprot:P07203").includes(path[0]!.from));
    // Preferring MANE picks the MANE Select protein among them, at the same cost.
    const [mane] = convert(stores, codon!.location, { to: { ref: "pdb:2F8A.A" }, prefer: ["MANE Select"] }, ctx);
    assert.equal(mane!.path[0]!.from, "ensembl:ENSP00000407375.1");
    assert.equal(mane!.cost, hits[0]!.cost);
    assert.deepEqual([hits[0]!.cost, hits[0]!.approximate], [2, false]);
  });

  it("does not expand structures when converting a protein to the genome (layer rule 1)", () => {
    const visited: string[] = [];
    const blocksAt = stores.blocksAt.bind(stores);
    stores.blocksAt = (ref, a, b) => (visited.push(ref), blocksAt(ref, a, b));
    const hits = convert(stores, parseLocationId("uniprot:P07203:49", ctx), { to: { category: "genome" } }, ctx);
    stores.blocksAt = blocksAt;
    assert.deepEqual(hits.map((h) => h.id), ["refseq:NC_000003.12:complement(49358132..49358134)"]);
    assert.ok(visited.length > 0 && visited.every((r) => !r.startsWith("pdb:")), visited.join(" "));
    assert.deepEqual(convert(stores, parseLocationId("uniprot:P07203:49", ctx), { to: { category: "structure" } }, ctx).map((h) => h.id).sort(), ["pdb:2F8A.A:59", "pdb:2F8A.B:59"]);
  });

  it("lists the MANE Select protein first among equal-cost results and reports its tags", () => {
    const proteins = convert(stores, parseLocationId("refseq:NC_000003.12:complement(49358132..49358134)", ctx), { to: { category: "protein" }, prefer: ["MANE Select"] }, ctx);
    assert.equal(proteins[0]!.location.outer, "ensembl:ENSP00000407375.1");
    assert.deepEqual(proteins[0]!.tags, ["MANE Select"]);
    assert.ok(proteins.length > 5 && proteins.slice(1).every((p) => p.tags.length === 0 || p.location.outer === "uniprot:P07203"));
  });

  it("maps a single base inside a codon to a codon position on the structure", () => {
    const hits = convert(stores, parseLocationId("refseq:NC_000003.12:49358133", ctx), { to: { ref: "pdb:2F8A.A" } }, ctx);
    assert.deepEqual(hits.map((h) => h.id), ["pdb:2F8A.A:59c2"]);
  });
});

describe("layer rules (spec-service §2.1)", () => {
  it("allows one U-turn (protein -> genome -> protein) but not a zigzag to a second genome", () => {
    const A = "refseq:NP_000011.1";
    const B = "refseq:NP_000012.1";
    const G1 = "refseq:NC_000011.1";
    const G2 = "refseq:NW_000012.1";
    const e = (from: string, to: string, tgt: number): Edge => ({
      kind: "annotation",
      from,
      to,
      blocks: [{ srcRef: from, src: 0, tgtRef: to, tgt, len: 30, rev: false }],
      attributes: {},
      provenance: { adapter: "gff3" },
      validation: { status: "ok" },
    });
    const stores = new StoreSet().add(store("zigzag", { sequences: [], annotations: [], warnings: [], edges: [e(A, G1, 100), e(B, G1, 100), e(B, G2, 500)] }));
    const ctx = stores.context();
    const ids = (to: Parameters<typeof convert>[2]) => convert(stores, parseLocationId(`${A}:2`, ctx), to, ctx).map((r) => r.id);
    assert.deepEqual(ids({ to: { category: "protein" } }), [`${B}:2`]); // A -> G1 -> B: one U-turn
    assert.deepEqual(ids({ to: { category: "genome" } }), [`${G1}:104..106`]); // A -> G1 -> B -> G2 would turn twice
  });
});

describe("MANE transcripts: RefSeq NM and Ensembl ENST are identical sequences (GPX1, real data)", async () => {
  const mane = new MemorySink();
  await ingestManeSummary(fileURLToPath(fixture("MANE.GRCh38.v1.5.summary_GPX1_RYBP.txt")), mane);
  const rna = new MemorySink();
  await ingestFastaFile(fileURLToPath(fixture("MANE_ensembl_rna_ENST00000419783.3.fa")), rna);
  await ingestFastaFile(fileURLToPath(fixture("NM_000581.4.fa")), rna);
  const stores = new StoreSet()
    .add(store("gpx1_mane", mane.result))
    .add(store("gpx1_mane_rna", rna.result))
    .add(store("gpx1_nm", ingestGenBank(await read("NM_000581.4.gb"))));
  const ctx = stores.context();

  it("tags both pairs, keeps the real lengths and converts NM <-> ENST by identity", () => {
    assert.deepEqual(stores.sequence("ensembl:ENST00000419783.3")?.tags, ["MANE Select"]);
    assert.equal(stores.sequence("refseq:NM_000581.4")?.length, 899); // MANE's unknown length (0) does not win
    assert.equal(stores.sequence("refseq:NM_000581.4")?.gene, "GPX1");
    const [hit] = convert(stores, parseLocationId("refseq:NP_000572.2:49", ctx), { to: { ref: "ensembl:ENST00000419783.3" } }, ctx);
    assert.equal(hit!.id, "ensembl:ENST00000419783.3:220..222");
    assert.deepEqual(hit!.path.map((s) => s.kind), ["annotation", "identity"]);
    assert.deepEqual(hit!.tags, ["MANE Select"]);
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
    assert.equal(categoryOf("pdb:2F8A.A", { unit: "aa" }), "structure");
  });
  it("bounds caches", () => {
    const lru = new Lru<string, number>(2).set("a", 1).set("b", 2);
    lru.get("a");
    lru.set("c", 3);
    assert.deepEqual([lru.has("a"), lru.has("b"), lru.has("c"), lru.size], [true, false, true, 2]);
  });
});
