import assert from "node:assert/strict";
import { mkdtempSync, readFileSync } from "node:fs";
import type { AddressInfo } from "node:net";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { after, before, describe, it } from "node:test";
import { ingestFastaFile, ingestGenBank, ingestManeSummary, MemorySink, SqliteSink, TogoCoordStore, type IngestResult } from "@togocoord/ingest";
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
    const mane = new MemorySink();
    await ingestManeSummary(fixture("MANE.GRCh38.v1.5.summary_GPX1_RYBP.txt"), mane);
    const stores = new StoreSet()
      .add(store("mane", mane.result))
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
    assert.equal(html.status, 303);
    assert.equal(html.location, `/?loc=${encodeURIComponent("refseq:NC_012920.1:join(16560..16569,1..>5)")}`);
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

  it("GET /v1/annotations returns features overlapping a location, not those spanning it with a gap", async () => {
    const { body } = await get(`/v1/annotations?loc=${encodeURIComponent("refseq:NC_012920.1:5")}`);
    assert.ok(body.annotations.some((a: { type: string }) => a.type === "D-loop"));
    // The D-loop is join(16024..16569,1..576): position 10000 lies inside its bounding interval only.
    const inner = await get(`/v1/annotations?loc=${encodeURIComponent("refseq:NC_012920.1:10000")}`);
    assert.ok(!inner.body.annotations.some((a: { type: string }) => a.type === "D-loop"));
    assert.ok(inner.body.annotations.some((a: { type: string; location: string }) => a.type === "tRNA" && a.location === "refseq:NC_012920.1:9991..10058"));
  });

  it("expands namespace:accession to the whole sequence and reports the explicit range", async () => {
    const { body } = await get(`/v1/convert?loc=${encodeURIComponent("uniprot:P07203")}&to=transcript`);
    assert.equal(body.input, "uniprot:P07203:1..203");
    assert.deepEqual(body.results.map((r: { location: string }) => r.location), ["refseq:NM_000581.4:76..684"]);
    const loc = await get(`/v1/location?loc=${encodeURIComponent("refseq:NC_012920.1")}`);
    assert.equal(loc.body.id, "refseq:NC_012920.1:1..16569");
    assert.equal((await get("/v1/convert?loc=refseq:NC_012920.1:16560..16570")).status, 400); // beyond the recorded length
  });

  it("resolves a sequence IRI (no location) to the sequence, not to a range", async () => {
    const json = await get("/refseq:NC_012920.1");
    assert.deepEqual([json.status, json.location], [303, "/v1/sequences/refseq%3ANC_012920.1"]);
    const html = await get("/refseq:NC_012920.1", "text/html");
    assert.equal(html.location, "/?loc=refseq%3ANC_012920.1");
  });

  it("limits the input length and the number of results", async () => {
    const stores = new StoreSet().add(new TogoCoordStore(join(dir, "mt.sqlite")));
    const server = createApi(stores, { maxInputLength: 1000, maxResults: 1 });
    await new Promise<void>((r) => server.listen(0, "127.0.0.1", r));
    const url = `http://127.0.0.1:${(server.address() as AddressInfo).port}`;
    try {
      const whole = await fetch(`${url}/v1/convert?loc=refseq:NC_012920.1`);
      assert.equal(whole.status, 413);
      assert.match((await whole.json()).error, /16569 bases\/residues; at most 1000/);
      const two = await (await fetch(`${url}/v1/convert?loc=${encodeURIComponent("refseq:NC_012920.1:8528")}`)).json();
      assert.deepEqual([two.results.length, two.truncated], [1, true]);
    } finally {
      server.close();
      stores.close();
    }
  });

  it("reports tags and filters by tag", async () => {
    const all = await get(`/v1/convert?loc=${encodeURIComponent("uniprot:P07203:49")}&to=protein`);
    const np = all.body.results.find((r: { sequence: string }) => r.sequence === "refseq:NP_000572.2");
    assert.deepEqual(np.tags, ["MANE Select"]);
    const tagged = await get(`/v1/convert?loc=${encodeURIComponent("refseq:NM_000581.4:220..222")}&to=protein&tag=${encodeURIComponent("MANE Select")}`);
    assert.deepEqual(tagged.body.results.map((r: { location: string }) => r.location), ["refseq:NP_000572.2:49"]);
    const seq = await get("/v1/sequences/refseq%3ANM_000581.4");
    assert.deepEqual([seq.body.length, seq.body.gene, seq.body.tags], [899, "GPX1", ["MANE Select"]]);
  });

  it("serves the web UI", async () => {
    const index = await get("/", "text/html");
    assert.equal(index.status, 200);
    assert.match(index.type!, /text\/html/);
    assert.match(index.body, /<script type="module" src="\/ui\/app.js">/);
    assert.match((await get("/ui/app.js")).type!, /javascript/);
    assert.match((await get("/ui/style.css")).type!, /text\/css/);
    assert.equal((await get("/ui/../src/api.ts")).status, 404);
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
