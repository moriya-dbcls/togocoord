import assert from "node:assert/strict";
import { mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { describe, it } from "node:test";
import { ingestBedFile } from "../src/adapter-bed.ts";
import { assemblyReportSeqids, lookupSeqid } from "../src/common.ts";
import { MemorySink } from "../src/model.ts";
import { fixture } from "./helpers.ts";

const names = assemblyReportSeqids(fixture("GRCh38.p14_assembly_report_chr13.txt"));
const refOf = (c: string) => names.get(c);

describe("BED adapter", () => {
  it("reads fanta.bio CREs (BED9+2, real data) as annotations of the genome", async () => {
    const path = new URL("./fixtures/fanta_human_CREv1.2.1_chr13_head.bed", import.meta.url).pathname;
    const sink = new MemorySink();
    const stats = await ingestBedFile(path, sink, { refOf, type: "CRE", extraColumns: ["Name", "attributes"] });
    assert.deepEqual([stats.records, stats.skipped], [3, 0]);
    const [first] = sink.result.annotations;
    assert.equal(first!.type, "CRE");
    assert.match(first!.location, /^refseq:NC_000013\.11:\d+\.\.\d+$/);
    assert.match(first!.attributes.ID![0]!, /^FCHS_\d+$/);
    assert.ok(first!.attributes.class && first!.attributes.directionality && first!.attributes.Name);
  });

  it("finds patches by the accession inside their UCSC name", () => {
    const m = new Map([["GL456060.1", "refseq:NW_1"], ["KI270706.1", "refseq:NT_2"]]);
    assert.deepEqual(
      ["chr11_GL456060_alt", "chr1_KI270706v1_random", "chr2_XX000001_fix"].map((n) => lookupSeqid(m, n)),
      ["refseq:NW_1", "refseq:NT_2", undefined],
    );
  });

  it("writes BED12 blocks on the minus strand as a complement join, and skips unknown chromosomes", async () => {
    const dir = mkdtempSync(join(tmpdir(), "togocoord-bed-"));
    const path = join(dir, "t.bed");
    writeFileSync(path, ["track name=x", "chr13\t100\t400\ttx1\t0\t-\t100\t400\t0\t2\t50,100,\t0,200,", "chrZ\t1\t2\tnope"].join("\n") + "\n");
    const sink = new MemorySink();
    const stats = await ingestBedFile(path, sink, { refOf });
    assert.deepEqual([stats.records, stats.skipped], [1, 1]);
    assert.equal(sink.result.annotations[0]!.location, "refseq:NC_000013.11:complement(join(101..150,301..400))");
    assert.equal(sink.result.warnings.length, 1);
  });
});
