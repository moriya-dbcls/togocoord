import assert from "node:assert/strict";
import { describe, it } from "node:test";
import {
  canonicalize,
  cdsMapping,
  compose,
  createContext,
  formatLocationId,
  LocationSemanticError,
  LocationSyntaxError,
  Mapping,
  mapLocation,
  MappingError,
  mappingFromLocation,
  parseLocationId,
} from "../src/index.ts";
import { testContext } from "./helpers.ts";

const ctx = testContext({ units: { "test:P": "aa" } });
const canon = (text: string) => canonicalize(text, ctx);

describe("INSDC location examples", () => {
  const same = [
    "test:A:467",
    "test:A:340..565",
    "test:A:<345..500",
    "test:A:<1..888",
    "test:A:1..>888",
    "test:A:102.110",
    "test:A:123^124",
    "test:A:join(12..78,134..202)",
    "test:A:complement(34..126)",
    "test:A:complement(join(2691..4571,4918..5163))",
    "test:A:join(1..100,B:100..202)",
    "test:A:order(1..10,20..30)",
  ];
  for (const text of same) it(`${text} is canonical`, () => assert.equal(canon(text), text));

  const rewrites: Array<[string, string]> = [
    ["test:A: join( 12..78 , 134..202 )", "test:A:join(12..78,134..202)"],
    ["test:A:5..5", "test:A:5"],
    ["test:A:join(complement(4918..5163),complement(2691..4571))", "test:A:complement(join(2691..4571,4918..5163))"],
    ["test:A:complement(join(complement(1..5),6..9))", "test:A:join(complement(6..9),1..5)"],
    ["test:A:join(join(1..2,3..4),5..6)", "test:A:join(1..2,3..4,5..6)"],
    ["test:A:A:1..5", "test:A:1..5"],
    ["test:A:join(1..5)", "test:A:1..5"],
    ["test:A:complement(<1..5)", "test:A:complement(<1..5)"],
    ["test:A:complement(B:1..5)", "test:A:complement(B:1..5)"],
    ["test:A:join(B:complement(1..5),B:complement(8..9))", "test:A:complement(join(B:8..9,B:1..5))"],
  ];
  for (const [input, expected] of rewrites) it(`${input} -> ${expected}`, () => assert.equal(canon(input), expected));
});

describe("protein locations and the codon extension", () => {
  const cases: Array<[string, string]> = [
    ["test:P:60", "test:P:60"],
    ["test:P:60c1..60c3", "test:P:60"],
    ["test:P:60c2", "test:P:60c2"],
    ["test:P:60c2..63c1", "test:P:60c2..63c1"],
    ["test:P:60c2..63c3", "test:P:60c2..63"],
    ["test:P:60c1..63", "test:P:60..63"],
    ["test:P:60^61", "test:P:60^61"],
    ["test:P:60c3^61c1", "test:P:60^61"],
    ["test:P:60c1^60c2", "test:P:60c1^60c2"],
    ["test:P:<1..5c2", "test:P:<1..5c2"],
  ];
  for (const [input, expected] of cases) it(`${input} -> ${expected}`, () => assert.equal(canon(input), expected));

  it("codon=never strips to covering residues", () => {
    const never = (t: string) => formatLocationId(parseLocationId(t, ctx), ctx, "never");
    assert.equal(never("test:P:60c2..63c1"), "test:P:60..63");
    assert.equal(never("test:P:60c2"), "test:P:60");
    assert.equal(never("test:P:60c1^60c2"), "test:P:60");
    assert.equal(never("test:P:60c3^61c1"), "test:P:60^61");
  });

  it("uses codon units internally", () => {
    assert.deepEqual(parseLocationId("test:P:2c2..3", ctx).segments, [{ ref: "test:P", start: 4, end: 9, strand: 1 }]);
  });
});

describe("errors", () => {
  const syntax = ["test:A:0", "test:A:-5", "test:A:12..5", "test:A:1..<5", "test:A:>1..5", "test:A:<5", "test:P:5c4", "test:A:join(1..2", "test:A:1..2)", "test:P:5c1.7"];
  for (const text of syntax) it(`syntax error: ${text}`, () => assert.throws(() => parseLocationId(text, ctx), LocationSyntaxError));

  const semantic = [
    "test:A:5^7",
    "test:A:5c1",
    "test:P:complement(1..5)",
    "test:A:join(order(1..2),3..4)",
    "nope:X:1",
    "refseq:not an accession:1",
    "insdc:NC_000001.11:join(1..5,test:B:1..2)",
  ];
  for (const text of semantic) it(`semantic error: ${text}`, () => assert.throws(() => parseLocationId(text, ctx), LocationSemanticError));

  it("range check when lengths are known", () => {
    const sized = testContext({ units: { "test:P": "aa" }, lengths: { "test:A": 100, "test:P": 10 } });
    assert.throws(() => parseLocationId("test:A:90..101", sized), LocationSemanticError);
    assert.throws(() => parseLocationId("test:P:11", sized), LocationSemanticError);
    assert.doesNotThrow(() => parseLocationId("test:P:10", sized));
  });

  it("refget sequences need an explicit unit", () => {
    const plain = createContext();
    assert.throws(() => parseLocationId("refget:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2:1..5", plain), LocationSemanticError);
    const withUnit = createContext({ units: { "refget:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2": "nt" } });
    assert.equal(
      canonicalize("refget:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2:1..5", withUnit),
      "refget:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2:1..5",
    );
  });
});

describe("namespaces", () => {
  const plain = createContext();
  it("infers units from accessions", () => {
    assert.equal(plain.unitOf("refseq:NP_055554.1"), "aa");
    assert.equal(plain.unitOf("refseq:NM_014739.3"), "nt");
    assert.equal(plain.unitOf("insdc:AAF99721.1"), "aa");
    assert.equal(plain.unitOf("insdc:J00194.1"), "nt");
    assert.equal(plain.unitOf("ensembl:ENSMUSP00000073626"), "aa");
    assert.equal(plain.unitOf("ensembl:ENST00000531224.5"), "nt");
  });
  it("accepts UniProt accessions and isoforms", () => {
    assert.equal(canonicalize("uniprot:Q9BYF1-1:60", plain), "uniprot:Q9BYF1-1:60");
    assert.equal(canonicalize("uniprot:A0A023GPI8:1..5", plain), "uniprot:A0A023GPI8:1..5");
  });
  it("treats '.' as the chain separator for PDB and normalises the ID", () => {
    assert.equal(canonicalize("pdb:4hhb.A:1..141", plain), "pdb:4HHB.A:1..141");
    assert.equal(canonicalize("pdb:pdb_00004hhb.H:5", plain), "pdb:4HHB.H:5");
    assert.equal(canonicalize("pdb:4HHB.h:5", plain), "pdb:4HHB.h:5");
    assert.throws(() => parseLocationId("pdb:4HHB:5", plain), LocationSemanticError);
  });
  it("lower-cases the prefix", () => {
    assert.equal(canonicalize("RefSeq:NM_014739.3:1..10", plain), "refseq:NM_014739.3:1..10");
  });
});

describe("composition through an mRNA", () => {
  it("protein -> mRNA -> genome equals the direct CDS mapping", () => {
    const c = testContext({ units: { "test:P": "aa" } });
    const mrnaToGenome = mappingFromLocation("test:M", parseLocationId("test:G:complement(join(1..20,31..40,51..70))", c));
    const proteinToMrna = cdsMapping({ protein: "test:P", cds: parseLocationId("test:M:11..37", c), aaLength: 8 });
    const direct = cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:complement(join(14..20,31..40,51..60))", c), aaLength: 8 });
    const viaMrna = compose(proteinToMrna, mrnaToGenome);
    for (let aa = 1; aa <= 8; aa++) {
      const loc = parseLocationId(`test:P:${aa}`, c);
      const a = mapLocation(loc, viaMrna, c).targets.map((t) => formatLocationId(t.location, c));
      const b = mapLocation(loc, direct, c).targets.map((t) => formatLocationId(t.location, c));
      assert.deepEqual(a, b, `aa ${aa}`);
    }
  });

  it("rejects a protein longer than its CDS", () => {
    const c = testContext({ units: { "test:P": "aa" } });
    assert.throws(() => cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:1..9", c), aaLength: 4 }), MappingError);
    assert.throws(() => cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:1..9", c), codonStart: 2, aaLength: 3 }), MappingError);
    assert.throws(() => new Mapping([{ srcRef: "a", src: 0, tgtRef: "b", tgt: 0, len: 0, rev: false }]), MappingError);
  });
});

describe("between-positions over overlapping features", () => {
  it("are mapped separately for each target sequence", () => {
    const c = testContext({ units: { "test:P": "aa", "test:Q": "aa" } });
    // Two overlapping CDSs in different frames on test:G.
    const m = cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:1..30", c), aaLength: 9 }).concat(
      cdsMapping({ protein: "test:Q", cds: parseLocationId("test:G:5..34", c), aaLength: 9 }),
    );
    const r = mapLocation(parseLocationId("test:G:6^7", c), m.inverse(), c);
    assert.deepEqual(r.targets.map((t) => formatLocationId(t.location, c)), ["test:P:2^3", "test:Q:1c2^1c3"]);
    assert.equal(r.unmapped, null);
    // Ambiguous within one sequence (ribosomal slippage re-reads a base): not mapped there.
    const slip = cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:join(1..6,6..17)", c), aaLength: 4 });
    const s = mapLocation(parseLocationId("test:G:5^6", c), slip.inverse(), c);
    assert.deepEqual(s.targets, []);
  });
});

describe("merging across join segments", () => {
  const c = testContext({ units: { "test:P": "aa" } });
  const cds = cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:join(11..20,31..40,51..57)", c), aaLength: 8 });
  const ids = (text: string, m = cds.inverse()) => mapLocation(parseLocationId(text, c), m, c).targets.map((t) => formatLocationId(t.location, c));
  it("maps a CDS location written per exon to one protein interval", () => {
    assert.deepEqual(ids("test:G:join(11..20,31..40,51..54)"), ["test:P:1..8"]);
  });
  it("round-trips a codon split across an exon junction to the same residue", () => {
    const genome = mapLocation(parseLocationId("test:P:4", c), cds, c).targets[0]!.location;
    assert.equal(formatLocationId(genome, c), "test:G:join(20,31..32)");
    assert.deepEqual(mapLocation(genome, cds.inverse(), c).targets.map((t) => formatLocationId(t.location, c)), ["test:P:4"]);
  });
  it("keeps order() segments and one-of segments apart", () => {
    assert.deepEqual(ids("test:G:order(11..20,31..40)"), ["test:P:order(1..4c1,4c2..7c2)"]);
  });
});

describe("cdsMapping with a leading partial codon (Ensembl X)", () => {
  it("maps residue 1 to the bases before the first complete codon", () => {
    const c = testContext({ units: { "test:P": "aa" } });
    // CDS of 11 nt with 2 leading bases (phase 2): X + 3 complete codons.
    const m = cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:101..111", c), codonStart: 3, aaLength: 4, leadingPartialCodon: true });
    const ids = (t: string) => mapLocation(parseLocationId(t, c), m, c).targets.map((x) => formatLocationId(x.location, c));
    // Only the last two codon positions exist in the CDS, so the begin is truncated.
    assert.deepEqual(ids("test:P:1"), ["test:G:<101..102"]);
    assert.deepEqual(ids("test:P:2"), ["test:G:103..105"]);
    assert.deepEqual(ids("test:P:4"), ["test:G:109..111"]);
    assert.deepEqual(mapLocation(parseLocationId("test:G:101", c), m.inverse(), c).targets.map((x) => formatLocationId(x.location, c)), ["test:P:1c2"]);
    assert.throws(() => cdsMapping({ protein: "test:P", cds: parseLocationId("test:G:101..111", c), aaLength: 3, leadingPartialCodon: true }), MappingError);
  });
});

describe("whole-sequence shorthand (namespace:accession)", () => {
  const sized = testContext({ units: { "test:P": "aa" }, lengths: { "test:A": 120, "test:P": 40 } });
  it("expands to 1..length when allowed and a length is known; formats explicitly", () => {
    assert.equal(formatLocationId(parseLocationId("test:A", sized, { wholeSequence: true }), sized), "test:A:1..120");
    assert.equal(formatLocationId(parseLocationId("test:P", sized, { wholeSequence: true }), sized), "test:P:1..40");
  });
  it("is rejected without the option or without a length", () => {
    assert.throws(() => parseLocationId("test:A", sized), LocationSyntaxError);
    assert.throws(() => parseLocationId("test:B", sized, { wholeSequence: true }), /length of 'test:B' is unknown/);
    assert.equal(formatLocationId(parseLocationId("test:A:5", sized, { wholeSequence: true }), sized), "test:A:5");
  });
});

