import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { createContext, encodeLocationId, parseLocationId, toFaldo } from "../src/index.ts";
import { testContext } from "./helpers.ts";

const ctx = testContext({ units: { "test:P": "aa" } });
const faldo = (text: string, c = ctx) => {
  const { "@context": _, ...node } = toFaldo(parseLocationId(text, c), c, { base: "https://t/" }) as Record<string, unknown>;
  return node;
};
const ref = (r: string) => ({ "@id": `https://identifiers.org/${r}` });
const pos = (types: string[], position: number, r = "test:A", extra: Record<string, unknown> = {}) => ({
  "@type": types,
  "faldo:position": position,
  "faldo:reference": ref(r),
  ...extra,
});
const F = "faldo:ForwardStrandPosition";
const R = "faldo:ReverseStrandPosition";
const E = "faldo:ExactPosition";

describe("FALDO JSON-LD", () => {
  it("reproduces the reverse-strand example of the FALDO README (cheY, complement(1965072..1965461))", () => {
    const c = createContext();
    assert.deepEqual(faldo("refseq:NC_000913.2:complement(1965072..1965461)", c), {
      "@id": "https://t/refseq:NC_000913.2:complement(1965072..1965461)",
      "@type": "faldo:Region",
      "faldo:begin": pos([E, R], 1965461, "refseq:NC_000913.2"),
      "faldo:end": pos([E, R], 1965072, "refseq:NC_000913.2"),
    });
  });

  it("writes single positions, fuzzy ends and uncertain positions", () => {
    assert.deepEqual(faldo("test:A:467"), { "@id": "https://t/test:A:467", ...pos([E, F], 467) });
    assert.deepEqual(faldo("test:A:<1..888")["faldo:begin"], pos(["faldo:FuzzyPosition", F], 1));
    assert.deepEqual(faldo("test:A:102.110"), {
      "@id": "https://t/test:A:102.110",
      "@type": "faldo:InRangePosition",
      "faldo:begin": pos([E, F], 102),
      "faldo:end": pos([E, F], 110),
    });
  });

  it("writes between-positions with after/before in the direction of the strand", () => {
    assert.deepEqual(faldo("test:A:123^124"), {
      "@id": "https://t/test:A:123%5E124",
      "@type": "faldo:InBetweenPosition",
      "faldo:after": pos([E, F], 123),
      "faldo:before": pos([E, F], 124),
    });
    const rev = faldo("test:A:complement(123^124)");
    assert.deepEqual([rev["faldo:after"], rev["faldo:before"]], [pos([E, R], 124), pos([E, R], 123)]);
  });

  it("writes join() as an ordered ListOfRegions in biological order and order() as a BagOfRegions", () => {
    const list = faldo("test:A:complement(join(2691..4571,4918..5163))");
    assert.equal(list["@type"], "faldo:ListOfRegions");
    assert.deepEqual(list["rdf:_1"], { "@type": "faldo:Region", "faldo:begin": pos([E, R], 5163), "faldo:end": pos([E, R], 4918) });
    assert.deepEqual(list["rdf:_2"], { "@type": "faldo:Region", "faldo:begin": pos([E, R], 4571), "faldo:end": pos([E, R], 2691) });
    assert.equal(faldo("test:A:order(1..10,20..30)")["@type"], "faldo:BagOfRegions");
  });

  it("writes protein positions without strand and with codon positions", () => {
    assert.deepEqual(faldo("test:P:60"), { "@id": "https://t/test:P:60", ...pos([E], 60, "test:P") });
    assert.deepEqual(faldo("test:P:60c2"), { "@id": "https://t/test:P:60c2", ...pos([E], 60, "test:P", { "tgc:codonPosition": 2 }) });
    const r = faldo("test:P:60c2..63c1");
    assert.deepEqual([r["faldo:begin"], r["faldo:end"]], [pos([E], 60, "test:P", { "tgc:codonPosition": 2 }), pos([E], 63, "test:P", { "tgc:codonPosition": 1 })]);
    const b = faldo("test:P:60c1^60c2");
    assert.deepEqual([b["faldo:after"], b["faldo:before"]], [pos([E], 60, "test:P", { "tgc:codonPosition": 1 }), pos([E], 60, "test:P", { "tgc:codonPosition": 2 })]);
    const whole = faldo("test:P:60^61");
    assert.deepEqual([whole["faldo:after"], whole["faldo:before"]], [pos([E], 60, "test:P"), pos([E], 61, "test:P")]);
  });

  it("declares its prefixes and percent-encodes only < > ^ in IRIs", () => {
    const node = toFaldo(parseLocationId("test:A:<1..>9", ctx), ctx);
    assert.deepEqual(Object.keys(node["@context"] as object), ["faldo", "rdf", "tgc"]);
    assert.equal(encodeLocationId("refseq:NM_1.1:join(<1..5,6^7,8..>9)"), "refseq:NM_1.1:join(%3C1..5,6%5E7,8..%3E9)");
  });
});
