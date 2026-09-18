// Streaming readers, FASTA index and SQLite store (scaling plan): each must reproduce the in-memory results.
import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import { mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { gzipSync } from "node:zlib";
import { describe, it } from "node:test";
import { createContext, formatLocationId, mapLocation, NamespaceRegistry, parseLocationId } from "@togocoord/core";
import { ingestGenBank } from "../src/adapter-gbff.ts";
import { ingestGff3 } from "../src/adapter-gff3.ts";
import { buildFai, FaiSequenceSource } from "../src/fasta-index.ts";
import { FeatureGrouper, parseGffLine } from "../src/gff3.ts";
import { edgeMapping, MemorySink } from "../src/model.ts";
import { SqliteSink, TogoCoordStore } from "../src/store.ts";
import { ingestGenBankFile, ingestGff3File } from "../src/stream.ts";
import { defaultFastaRef, ingestFastaFile } from "../src/adapter-fasta.ts";
import { ingestSiftsFile } from "../src/adapter-sifts.ts";
import { ingestChainFile } from "../src/adapter-chain.ts";
import { assemblyReportSeqids } from "../src/common.ts";
import { parseFastaHeaders } from "../src/fasta.ts";
import { MemorySequenceSource } from "../src/sequence.ts";
import { fixture, genbankSource } from "./helpers.ts";

const dir = mkdtempSync(join(tmpdir(), "togocoord-"));
const fixturePath = (name: string) => fileURLToPath(new URL(`./fixtures/${name}`, import.meta.url));

describe("streaming readers reproduce the in-memory adapters", () => {
  for (const name of ["NC_012920.1.gff3", "NC_001405.1.gff3", "NC_000003.12_cDNA_match.gff3", "NC_000003.12_RYBP_GPX1.gff3"]) {
    it(`GFF3 ${name} (plain and gzip)`, async () => {
      const expected = ingestGff3(fixture(name));
      const gz = join(dir, `${name}.gz`);
      writeFileSync(gz, gzipSync(fixture(name)));
      for (const path of [fixturePath(name), gz]) {
        const sink = new MemorySink();
        await ingestGff3File(path, sink);
        assert.deepEqual(sink.result, expected);
      }
    });
  }

  it("GenBank, several records in one gzipped file", async () => {
    const text = ["NC_045512.2.gb", "NM_000581.4.gb", "NP_055554.1.gp"].map(fixture).join("");
    const gz = join(dir, "multi.gb.gz");
    writeFileSync(gz, gzipSync(text));
    const sink = new MemorySink();
    const stats = await ingestGenBankFile(gz, sink);
    assert.deepEqual(sink.result, ingestGenBank(text));
    assert.ok(stats.features > 100);
  });

  it("exon features are not stored as annotations by default", () => {
    const withExons = ingestGff3(fixture("NC_001405.1.gff3"), { excludeAnnotations: new Set() });
    const without = ingestGff3(fixture("NC_001405.1.gff3"));
    const exons = withExons.annotations.filter((a) => a.type === "exon").length;
    assert.ok(exons > 0);
    assert.equal(without.annotations.length, withExons.annotations.length - exons);
    assert.deepEqual(without.edges, withExons.edges);
  });
});

describe("FeatureGrouper", () => {
  const row = (id: string, n: number) => parseGffLine(`s\tx\tCDS\t${n}\t${n + 1}\t.\t+\t0\tID=${id}`, n)!;
  it("groups adjacent rows and reports rows of an already emitted feature", () => {
    const g = new FeatureGrouper(1);
    const out = [...g.push(row("a", 1)), ...g.push(row("a", 3)), ...g.push(row("b", 5)), ...g.push(row("a", 7)), ...g.end()];
    assert.deepEqual(out.map((f) => [f.id, f.rows.length]), [["a", 2], ["b", 1], ["a", 1]]);
    assert.equal(g.violations.length, 1);
    assert.match(g.violations[0]!, /line 7: rows of s a are not adjacent/);
  });
});

describe("FASTA index", () => {
  const residues = "ACGTACGTTTGCAACGTAGGGCATCGATCGAT"; // 32
  const wrap = (s: string, w: number, eol: string) => s.match(new RegExp(`.{1,${w}}`, "g"))!.join(eol) + eol;
  for (const eol of ["\n", "\r\n"]) {
    it(`random access equals the in-memory sequence (${JSON.stringify(eol)} line ends)`, () => {
      const fa = join(dir, `t${eol.length}.fa`);
      writeFileSync(fa, `>NC_000001.1 first${eol}${wrap(residues, 7, eol)}>NC_000002.1${eol}${wrap("NNACGT", 4, eol)}`);
      const entries = buildFai(fa);
      assert.deepEqual(entries.map((e) => [e.name, e.length, e.lineBases, e.lineBytes]), [
        ["NC_000001.1", 32, 7, 7 + eol.length],
        ["NC_000002.1", 6, 4, 4 + eol.length],
      ]);
      const src = new FaiSequenceSource(fa, (n) => `refseq:${n}`, entries);
      for (let a = 0; a <= 32; a++) for (let b = a; b <= 32; b++) assert.equal(src.get("refseq:NC_000001.1", a, b), residues.slice(a, b));
      assert.equal(src.get("refseq:NC_000002.1", 1, 6), "NACGT");
      assert.equal(src.get("refseq:NC_000001.1", 30, 33), undefined);
      src.close();
    });
  }
  it("rejects irregular line lengths", () => {
    const fa = join(dir, "bad.fa");
    writeFileSync(fa, ">x\nACGT\nAC\nACGT\n");
    assert.throws(() => buildFai(fa), /irregular line length/);
  });
});

describe("SQLite store", () => {
  const path = join(dir, "mt.sqlite");
  const sink = new SqliteSink(path);
  const text = fixture("NC_012920.1.gb");
  const memory = ingestGenBank(text);
  for (const s of memory.sequences) sink.sequence(s);
  for (const e of memory.edges) sink.edge(e);
  for (const a of memory.annotations) sink.annotation(a);
  sink.close({ source: "NC_012920.1.gb" });
  const store = new TogoCoordStore(path);
  const ctx = store.context();
  const ids = (text: string) => store.neighbors(parseLocationId(text, ctx), ctx).targets.map((t) => formatLocationId(t.location, ctx)).sort();

  it("stores everything and records metadata", () => {
    assert.deepEqual(store.counts(), {
      sequence: memory.sequences.length,
      edge: memory.edges.length,
      block: memory.edges.reduce((n, e) => n + e.blocks.length, 0),
      annotation: memory.annotations.length,
      warning: 0,
    });
    assert.equal(store.meta().source, "NC_012920.1.gb");
    assert.equal(store.sequence("refseq:NC_012920.1")?.topology, "circular");
    assert.equal(store.unitOf("refseq:YP_003024037.1"), "aa");
  });

  it("converts in both directions using only the blocks it fetches", () => {
    assert.deepEqual(ids("refseq:YP_003024037.1:174"), ["refseq:NC_012920.1:complement(14152..14154)"]);
    // ATP8 (8366..8572) and ATP6 (8527..9207) overlap in different reading frames:
    // nt 8527 is unit 161 of ATP8 (residue 54, codon position 3) and the first base of ATP6.
    assert.deepEqual(ids("refseq:NC_012920.1:8527..8529"), ["refseq:YP_003024030.1:54c3..55c2", "refseq:YP_003024031.1:1"]);
    assert.deepEqual(ids("refseq:NC_012920.1:8526^8527"), ["refseq:YP_003024030.1:54c2^54c3"]);
  });

  it("agrees with the in-memory edge mapping for every residue of every protein", () => {
    for (const e of memory.edges) {
      const direct = edgeMapping(e);
      for (let r = 1; r <= Number(e.attributes.aaLength); r++) {
        const loc = parseLocationId(`${e.from}:${r}`, ctx);
        const a = mapLocation(loc, direct, ctx).targets.map((t) => formatLocationId(t.location, ctx));
        const b = mapLocation(loc, store.mappingFor(loc, { direction: "forward" }), ctx).targets.map((t) => formatLocationId(t.location, ctx));
        assert.deepEqual(b, a, `${e.from}:${r}`);
      }
    }
  });

  it("finds annotations by interval, including origin-spanning ones", () => {
    const types = store.annotations("refseq:NC_012920.1", 0, 10).map((a) => a.type);
    assert.ok(types.includes("D-loop"));
    assert.ok(store.edges("refseq:YP_003024037.1").some((e) => e.to === "refseq:NC_012920.1" && e.validation.status === "ok"));
  });
});

it("togocoord-ingest --db builds a store from gzipped GFF3 with FASTA validation", () => {
  const cli = fixturePath("../../src/cli.ts");
  const gff = join(dir, "sars2.gff3.gz");
  writeFileSync(gff, gzipSync(fixture("NC_045512.2.gff3")));
  const fa = join(dir, "sars2.fa");
  writeFileSync(fa, `>NC_045512.2\n${genbankSource("NC_045512.2.gb").get("refseq:NC_045512.2", 0, 29903)}\n`);
  const db = join(dir, "sars2.sqlite");
  const run = spawnSync(process.execPath, ["--no-warnings", cli, "--db", db, "--fasta", fa, gff], { encoding: "utf8" });
  assert.equal(run.status, 0, run.stderr);
  assert.match(run.stderr, /edges 12 \{"ok":12\}/);
  const store = new TogoCoordStore(db);
  const ctx = store.context();
  const r = store.neighbors(parseLocationId("refseq:NC_045512.2:13468", ctx), ctx);
  // nt 13468 is in ORF1ab (both sides of the slippage) and in ORF1a (266..13483).
  assert.deepEqual(r.targets.map((t) => formatLocationId(t.location, ctx)).sort(), [
    "refseq:YP_009724389.1:join(4401c3,4402c1)",
    "refseq:YP_009725295.1:4401c3",
  ]);
  store.close();
});

it("FASTA records add checksums to sequences first seen without residues", async () => {
  const path = join(dir, "fasta.sqlite");
  const sink = new SqliteSink(path);
  sink.sequence({ ref: "uniprot:P07203", moltype: "protein", unit: "aa", length: 203, provenance: { adapter: "gff3" } });
  await ingestFastaFile(fixturePath("uniprot_P07203.fa"), sink);
  sink.close();
  const store = new TogoCoordStore(path);
  const seq = store.sequence("uniprot:P07203")!;
  assert.deepEqual([seq.length, seq.digest, seq.md5], [203, "SQ.qsfE5UDDg5US4PKbR3-mgCvYbsfBsxFH", "b3ac19c6abd503ac0cc8c151701eb02a"]);
  // UniParc reports MD5 in upper case (UPI00001B07C3: B3AC19C6ABD503AC0CC8C151701EB02A).
  assert.deepEqual(store.refsByMd5("B3AC19C6ABD503AC0CC8C151701EB02A"), ["uniprot:P07203"]);
  assert.deepEqual(store.refsByDigest("SQ.qsfE5UDDg5US4PKbR3-mgCvYbsfBsxFH"), ["uniprot:P07203"]);
  store.close();
});

describe("SIFTS (UniProt P07203 <-> PDB 2F8A, real data)", () => {
  it("builds residue blocks on SEQRES numbering and validates against pdb_seqres", async () => {
    const source = new MemorySequenceSource();
    for (const [header, s] of parseFastaHeaders(fixture("uniprot_P07203.fa") + fixture("pdb_seqres_P07203.txt"))) {
      const ref = defaultFastaRef(header, new NamespaceRegistry());
      if (ref) source.add(ref, s);
    }
    const sink = new MemorySink();
    const stats = await ingestSiftsFile(fixturePath("sifts_P07203.tsv"), sink, { source });
    assert.deepEqual(stats, { rows: 2, edges: 2, skipped: 0 });
    const [a] = sink.result.edges;
    assert.deepEqual([a!.from, a!.to, a!.kind], ["uniprot:P07203", "pdb:2F8A.A", "alignment"]);
    assert.deepEqual(a!.blocks, [{ srcRef: "uniprot:P07203", src: 39, tgtRef: "pdb:2F8A.A", tgt: 69, len: 552, rev: false }]);
    // 2F8A is the U49G mutant: 183 of 184 aligned residues identical.
    assert.deepEqual([a!.validation.status, a!.validation.detail], ["ok", "183/184 aligned residues identical"]);
    assert.equal(a!.attributes.authorNumbering, "12-195");
    const ctx = createContext({ units: { "pdb:2F8A.A": "aa" } });
    const r = mapLocation(parseLocationId("uniprot:P07203:49", ctx), edgeMapping(a!), ctx);
    assert.deepEqual(r.targets.map((t) => formatLocationId(t.location, ctx)), ["pdb:2F8A.A:59"]);
    assert.equal(source.get("pdb:2F8A.A", 58, 59), "G");
  });
});

describe("UCSC liftOver chain (hg38 -> mm39, chain 8179: human chr13 + to mouse chr3 -, real data)", () => {
  const human = assemblyReportSeqids(fixture("GRCh38.p14_assembly_report_chr13.txt"));
  const mouse = assemblyReportSeqids(fixture("GRCm39_assembly_report_chr3.txt"));
  const source = new MemorySequenceSource();
  for (const [, s] of parseFastaHeaders(fixture("NC_000013.11_65787513-65788763.fa"))) source.add("refseq:NC_000013.11", s, 65787512);
  for (const [, s] of parseFastaHeaders(fixture("NC_000069.7_54389129-54390378.fa"))) source.add("refseq:NC_000069.7", s, 54389128);

  it("converts '-' query coordinates to forward coordinates and validates by sequence identity", async () => {
    const sink = new MemorySink();
    const stats = await ingestChainFile(fixturePath("hg38ToMm39_chain8179.chain"), sink, {
      source,
      fromRef: (n) => human.get(n),
      toRef: (n) => mouse.get(n),
    });
    assert.deepEqual([stats.chains, stats.blocks, stats.skipped], [1, 2, 0]);
    const [e] = sink.result.edges;
    assert.deepEqual([e!.kind, e!.directional, e!.from, e!.to], ["liftover", true, "refseq:NC_000013.11", "refseq:NC_000069.7"]);
    // qSize 159745316: reverse block [105354938, +16) is forward [54390362, 54390378).
    assert.deepEqual(e!.blocks, [
      { srcRef: "refseq:NC_000013.11", src: 65787512, tgtRef: "refseq:NC_000069.7", tgt: 54390362, len: 16, rev: true },
      { srcRef: "refseq:NC_000013.11", src: 65787529, tgtRef: "refseq:NC_000069.7", tgt: 54389128, len: 1234, rev: true },
    ]);
    assert.equal(e!.validation.status, "ok");
    const [same, total] = e!.validation.detail!.match(/\d+/g)!.map(Number);
    assert.ok(same! / total! > 0.6, e!.validation.detail); // orthologous DNA; misplaced coordinates would give ~0.25

    // Stored in chunks, used forward only.
    const path = join(dir, "chain.sqlite");
    const sq = new SqliteSink(path);
    sq.edge(e!);
    sq.close();
    const store = new TogoCoordStore(path);
    const fwd = store.blocksAt("refseq:NC_000013.11", 65787520, 65787540, "forward");
    assert.deepEqual(fwd.map(({ edge: _e, ...b }) => b), e!.blocks);
    assert.deepEqual(store.blocksAt("refseq:NC_000069.7", 54389128, 54390378), []); // not used in reverse
    assert.equal(store.summary().blocks, 2);
    const ctx = createContext({ units: () => "nt" });
    const r = mapLocation(parseLocationId("refseq:NC_000013.11:65787530..65787532", ctx), store.mappingFor(parseLocationId("refseq:NC_000013.11:65787530..65787532", ctx)), ctx);
    assert.deepEqual(r.targets.map((t) => formatLocationId(t.location, ctx)), ["refseq:NC_000069.7:complement(54390360..54390362)"]);
    store.close();
  });

  it("round-trips chunk encoding for many blocks with negative target steps", async () => {
    const blocks = Array.from({ length: 700 }, (_, i) => ({ srcRef: "refseq:NC_000001.1", src: 1000 + 10 * i, tgtRef: "refseq:NC_000002.1", tgt: 900000 - 12 * i, len: 5 + (i % 4), rev: true }));
    const path = join(dir, "chunks.sqlite");
    const sq = new SqliteSink(path);
    sq.edge({ kind: "liftover", directional: true, from: "refseq:NC_000001.1", to: "refseq:NC_000002.1", blocks, attributes: {}, provenance: { adapter: "chain" }, validation: { status: "skipped" } });
    sq.close();
    const store = new TogoCoordStore(path);
    assert.deepEqual(store.blocksAt("refseq:NC_000001.1", 0, 10_000_000, "forward").map(({ edge: _e, ...b }) => b), blocks);
    assert.deepEqual(store.blocksAt("refseq:NC_000001.1", 3500, 3501, "forward").map((b) => b.src), [3500]);
    store.close();
  });
});
