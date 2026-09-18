import assert from "node:assert/strict";
import { mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { describe, it } from "node:test";
import { ingestPafFile, pafBlocks, parsePafLine } from "../src/adapter-paf.ts";
import { MemorySink } from "../src/model.ts";
import { MemorySequenceSource, reverseComplement } from "../src/sequence.ts";

// A target T and a query Q made from it with a 2-base insertion and a 3-base deletion; R is Q reverse-complemented.
let seed = 7;
const T = Array.from({ length: 200 }, () => "ACGT"[(seed = (seed * 1103515245 + 12345) >>> 0) % 4]).join("");
const Q = T.slice(0, 50) + "GG" + T.slice(50, 100) + T.slice(103, 150);
const R = reverseComplement(Q);
const T2 = T.slice(120, 170); // a second target, for the overlap test
const row = (q: string, qlen: number, qs: number, qe: number, strand: string, t: string, tlen: number, ts: number, te: number, score: number, cigar: string) =>
  [q, qlen, qs, qe, strand, t, tlen, ts, te, score, qe - qs, 60, "tp:A:P", `AS:i:${score}`, `cg:Z:${cigar}`].join("\t");
const ref = (n: string) => `insdc:${n}`;
const source = new MemorySequenceSource();
for (const [n, s] of [["Q", Q], ["R", R], ["T", T], ["T2", T2]] as const) source.add(ref(n), s);
const aligned = (b: { src: number; tgt: number; len: number; rev: boolean; srcRef: string; tgtRef: string }) => {
  const s = source.get(b.srcRef, b.src, b.src + b.len)!;
  return (b.rev ? reverseComplement(s) : s) === source.get(b.tgtRef, b.tgt, b.tgt + b.len);
};

describe("PAF adapter", () => {
  it("turns a CIGAR into ungapped blocks on the forward strand", () => {
    const r = parsePafLine(row("Q", Q.length, 0, 149, "+", "T", 200, 0, 150, 147, "50=2I50=3D47="))!;
    const blocks = pafBlocks(r, ref("Q"), ref("T"));
    assert.deepEqual(blocks.map((b) => [b.src, b.tgt, b.len, b.rev]), [[0, 0, 50, false], [52, 50, 50, false], [102, 103, 47, false]]);
    assert.ok(blocks.every(aligned));
  });

  it("walks the query backwards on the reverse strand", () => {
    const r = parsePafLine(row("R", R.length, 0, 149, "-", "T", 200, 0, 150, 147, "25=25=2I50=3D47="))!; // adjacent matches merge
    const blocks = pafBlocks(r, ref("R"), ref("T"));
    assert.deepEqual(blocks.map((b) => [b.src, b.tgt, b.len, b.rev]), [[99, 0, 50, true], [47, 50, 50, true], [0, 103, 47, true]]);
    assert.ok(blocks.every(aligned));
  });

  it("keeps one alignment per source base, best first, and validates against the sequences", async () => {
    const dir = mkdtempSync(join(tmpdir(), "togocoord-paf-"));
    const path = join(dir, "q_to_t.paf");
    writeFileSync(
      path,
      [
        row("Q", Q.length, 100, 149, "+", "T2", 50, 0, 49, 40, "49="), // weaker, overlaps the first one on Q[100, 120)
        row("Q", Q.length, 0, 120, "+", "T", 200, 0, 121, 110, "50=2I50=3D18="),
        row("Q", Q.length, 0, 50, "+", "T", 200, 0, 50, 45, "50=").replace("tp:A:P", "tp:A:S"), // secondary: ignored
      ].join("\n") + "\n",
    );
    const sink = new MemorySink();
    const stats = await ingestPafFile(path, sink, { fromRef: ref, toRef: ref, source, minLength: 10 });
    assert.deepEqual([stats.alignments, stats.skipped, stats.overlapBases], [2, 1, 20]);
    const [best, trimmed] = sink.result.edges;
    assert.deepEqual([best!.to, best!.validation.status, best!.directional], [ref("T"), "ok", true]);
    assert.deepEqual(trimmed!.blocks.map((b) => [b.src, b.tgt, b.len]), [[120, 20, 29]]); // Q[120, 149) only
    assert.ok([...best!.blocks, ...trimmed!.blocks].every(aligned));
    assert.equal(stats.identicalBases, stats.sampledBases);
  });

  it("requires the CIGAR", () => {
    const r = parsePafLine(row("Q", Q.length, 0, 149, "+", "T", 200, 0, 150, 147, "x").replace("\tcg:Z:x", ""))!;
    assert.throws(() => pafBlocks(r, ref("Q"), ref("T")), /cg:Z/);
  });
});
