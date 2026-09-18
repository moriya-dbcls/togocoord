import assert from "node:assert/strict";
import { mkdtempSync, readFileSync } from "node:fs";
import type { AddressInfo } from "node:net";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { after, before, describe, it } from "node:test";
import { ingestFastaFile, ingestGenBank, MemorySink, SqliteSink, TogoCoordStore, type IngestResult } from "@togocoord/ingest";
import { createApi, StoreSet } from "../src/index.ts";

const dir = mkdtempSync(join(tmpdir(), "togocoord-api-"));
const fixture = (name: string) => fileURLToPath(new URL(`../../ingest/test/fixtures/${name}`, import.meta.url));

function store(name: string, r: IngestResult): TogoCoordStore {
  const path = join(dir, `${name}.sqlite`);
  const sink = new SqliteSink(path);
  for (const s of r.sequences) sink.sequence(s);
  for (const e of r.edges) sink.edge(e);
  for (const a of r.annotations) sink.annotation(a);
  sink.close();
  return new TogoCoordStore(path);
}

describe("REST API (GPX1 mRNA + UniProt P07203 + human mtDNA, real data)", () => {
  let base = "";
  let close = () => {};
  before(async () => {
    const uniprot = new MemorySink();
    await ingestFastaFile(fixture("uniprot_P07203.fa"), uniprot);
    const stores = new StoreSet()
      .add(store("gpx1", ingestGenBank(readFileSync(fixture("NM_000581.4.gb"), "utf8"))))
      .add(store("uniprot", uniprot.result))
      .add(store("mt", ingestGenBank(readFileSync(fixture("NC_012920.1.gb"), "utf8"))));
    const server = createApi(stores, { base: "https://t/" });
    await new Promise<void>((r) => server.listen(0, "127.0.0.1", r));
    base = `http://127.0.0.1:${(server.address() as AddressInfo).port}`;
    close = () => {
      server.close();
      stores.close();
    };
  });
  after(() => close());
  const get = async (path: string, accept = "application/json") => {
    const res = await fetch(base + path, { headers: { accept }, redirect: "manual" });
    const text = await res.text();
    return { status: res.status, type: res.headers.get("content-type"), location: res.headers.get("location"), body: text && res.headers.get("content-type")?.includes("json") ? JSON.parse(text) : text };
  };

  it("GET /v1/convert follows identity and CDS edges", async () => {
    const { status, body } = await get(`/v1/convert?loc=${encodeURIComponent("uniprot:P07203:49")}&to=transcript`);
    assert.equal(status, 200);
    assert.equal(body.input, "uniprot:P07203:49");
    assert.deepEqual(
      body.results.map((r: { location: string; cost: number; approximate: boolean }) => [r.location, r.cost, r.approximate]),
      [["refseq:NM_000581.4:220..222", 1, false]],
    );
    assert.deepEqual(body.results[0].path.map((s: { kind: string }) => s.kind), ["identity", "annotation"]);
    assert.equal(body.results[0].iri, "https://t/refseq:NM_000581.4:220..222");
  });

  it("GET /v1/convert without a target lists directly connected sequences; codon=never strips codon positions", async () => {
    const { body } = await get(`/v1/convert?loc=${encodeURIComponent("refseq:NC_012920.1:8528")}`);
    assert.deepEqual(body.results.map((r: { location: string }) => r.location).sort(), ["refseq:YP_003024030.1:55c1", "refseq:YP_003024031.1:1c2"]);
    const never = await get(`/v1/convert?loc=${encodeURIComponent("refseq:NC_012920.1:8528")}&codon=never`);
    assert.deepEqual(never.body.results.map((r: { location: string }) => r.location).sort(), ["refseq:YP_003024030.1:55", "refseq:YP_003024031.1:1"]);
  });

  it("POST /v1/convert handles batches and reports per-location errors", async () => {
    const res = await fetch(`${base}/v1/convert`, {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ locations: ["refseq:YP_003024037.1:174", "refseq:NC_012920.1:12..5"], to: "genome" }),
    });
    const body = await res.json();
    assert.equal(res.status, 200);
    assert.equal(body.results[0].results[0].location, "refseq:NC_012920.1:complement(14152..14154)");
    assert.match(body.results[1].error, /range begin is greater than end/);
    assert.equal(body.results[1].position, 0);
  });

  it("GET /v1/location and /v1/location/faldo describe a location", async () => {
    const { body } = await get(`/v1/location?loc=${encodeURIComponent("uniprot:P07203:48c3..49")}`);
    assert.equal(body.id, "uniprot:P07203:48c3..49");
    assert.deepEqual(body.segments, [
      { sequence: "uniprot:P07203", strand: "+", begin: { residue: 48, codonPosition: 3 }, end: { residue: 49, codonPosition: 3 } },
    ]);
    const ld = await get(`/v1/location/faldo?loc=${encodeURIComponent("refseq:NC_012920.1:complement(14152..14154)")}`);
    assert.equal(ld.type, "application/ld+json");
    assert.equal(ld.body["@type"], "faldo:Region");
    assert.equal(ld.body["faldo:begin"]["faldo:position"], 14154);
  });

  it("resolves location IRIs by content negotiation", async () => {
    const path = "/refseq:NC_012920.1:join(16560..16569,1..%3E5)";
    const ld = await get(path, "application/ld+json");
    assert.equal(ld.body["@id"], "https://t/refseq:NC_012920.1:join(16560..16569,1..%3E5)");
    const html = await get(path, "text/html");
    assert.match(html.body, /<h1>refseq:NC_012920.1:join\(16560..16569,1..&gt;5\)<\/h1>/);
    const json = await get(path, "application/json");
    assert.equal(json.status, 303);
    assert.equal(json.location, `/v1/location?loc=${encodeURIComponent("refseq:NC_012920.1:join(16560..16569,1..>5)")}`);
  });

  it("GET /v1/sequences/{ref} merges records across stores; /edges lists edges", async () => {
    const { body } = await get("/v1/sequences/uniprot%3AP07203");
    assert.deepEqual([body.length, body.category, body.identical], [203, "protein", ["refseq:NP_000572.2"]]);
    const edges = await get("/v1/sequences/refseq%3ANP_000572.2/edges");
    assert.deepEqual(edges.body.edges.map((e: { to: string; validation: { status: string } }) => [e.to, e.validation.status]), [["refseq:NM_000581.4", "ok"]]);
  });

  it("GET /v1/annotations returns features overlapping a location", async () => {
    const { body } = await get(`/v1/annotations?loc=${encodeURIComponent("refseq:NC_012920.1:5")}`);
    assert.ok(body.annotations.some((a: { type: string }) => a.type === "D-loop"));
  });

  it("returns JSON errors", async () => {
    assert.equal((await get("/v1/convert?loc=refseq:NC_012920.1:0")).status, 400);
    assert.equal((await get("/v1/convert?loc=refseq:NC_012920.1:5&to=nowhere")).status, 400);
    assert.equal((await get("/v1/nothing")).status, 404);
    assert.equal((await get("/v1/sequences/refseq%3ANM_999999.1")).status, 404);
    const post = await fetch(`${base}/v1/meta`, { method: "POST" });
    assert.equal(post.status, 405);
    assert.equal(post.headers.get("access-control-allow-origin"), "*");
  });
});
