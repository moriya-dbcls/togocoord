import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { alignProteins } from "../src/adapter-protein-align.ts";

let seed = 3;
const protein = (n: number) => Array.from({ length: n }, () => "ACDEFGHIKLMNPQRSTVWY"[((seed = (seed * 1103515245 + 12345) >>> 0) >>> 16) % 20]).join("");
const B = protein(189);

describe("protein alignment (T2)", () => {
  it("aligns substitutions residue by residue (like Swiss-Prot P08556 and GRCm39's NP_035067)", () => {
    const a = B.slice(0, 167) + (B[167] === "L" ? "M" : "L") + B.slice(168, 183) + (B[183] === "S" ? "T" : "S") + B.slice(184);
    const al = alignProteins(a, B)!;
    assert.deepEqual(al.blocks, [[0, 0, 189]]);
    assert.deepEqual(al.substitutions.map(([i, j]) => [i, j]), [[168, 168], [184, 184]]);
    assert.equal(al.identical, 187);
  });

  it("splits at an insertion and keeps the numbering on both sides", () => {
    const a = B.slice(0, 100) + "WWWWW" + B.slice(100);
    const al = alignProteins(a, B)!;
    assert.deepEqual(al.blocks, [[0, 0, 100], [105, 100, 89]]);
    assert.equal(al.substitutions.length, 0);
  });

  it("leaves an extra N-terminus unaligned (a different first exon)", () => {
    const a = protein(40) + B.slice(20);
    const al = alignProteins(a, B)!;
    assert.deepEqual(al.blocks.at(-1), [40, 20, 169]);
    assert.ok(al.aligned < 169 + 20); // the unrelated 40 residues do not pair up with B's first 20
  });

  it("finds nothing between unrelated proteins", () => {
    assert.equal(alignProteins(protein(150), protein(150)), undefined);
  });
});
