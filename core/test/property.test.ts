// Property-based checks against the naive oracle (design §10).
import assert from "node:assert/strict";
import { describe, it } from "node:test";
import fc from "fast-check";
import {
  compose,
  formatLocationId,
  invert,
  Mapping,
  mapLocation,
  parseLocationId,
  type Block,
  type Location,
  type Segment,
} from "../src/index.ts";
import { testContext } from "./helpers.ts";
import { actualPairs, expectedPairs, mapUnit, traversal } from "./oracle.ts";

const RUNS = { numRuns: 500 };
const ctx = testContext({ units: { "test:P": "aa", "test:Q": "aa" } });

const arbBlock = (srcRefs: string[], tgtRefs: string[]): fc.Arbitrary<Block> =>
  fc.record({
    srcRef: fc.constantFrom(...srcRefs),
    src: fc.nat(40),
    tgtRef: fc.constantFrom(...tgtRefs),
    tgt: fc.nat(40),
    len: fc.integer({ min: 1, max: 12 }),
    rev: fc.boolean(),
  });

const arbMapping = (srcRefs: string[], tgtRefs: string[]) =>
  fc.array(arbBlock(srcRefs, tgtRefs), { minLength: 0, maxLength: 8 }).map((bs) => new Mapping(bs));

const arbSegment = (ref: string, stranded = true): fc.Arbitrary<Segment> =>
  fc
    .record({
      start: fc.nat(50),
      len: fc.integer({ min: 1, max: 15 }),
      strand: stranded ? fc.constantFrom<1 | -1>(1, -1) : fc.constant<1 | -1>(1),
      fuzzyLow: fc.boolean(),
      fuzzyHigh: fc.boolean(),
    })
    .map(({ start, len, strand, fuzzyLow, fuzzyHigh }) => {
      const seg: Segment = { ref, start, end: start + len, strand };
      if (fuzzyLow) seg.fuzzyLow = true;
      if (fuzzyHigh) seg.fuzzyHigh = true;
      return seg;
    });

const arbLocation = (ref: string, stranded = true): fc.Arbitrary<Location> =>
  fc
    .array(arbSegment(ref, stranded), { minLength: 1, maxLength: 4 })
    .map((segments) => ({ outer: ref, kind: "join" as const, segments }));

const sortImages = (xs: Array<{ ref: string; pos: number; rev: boolean }>) =>
  xs.map((x) => `${x.ref}|${x.pos}|${x.rev}`).sort();

describe("mapLocation agrees with the per-unit oracle", () => {
  it("nt -> nt, arbitrary overlapping blocks", () => {
    fc.assert(
      fc.property(arbMapping(["test:A"], ["test:X", "test:Y"]), arbLocation("test:A"), (m, loc) => {
        const result = mapLocation(loc, m, ctx);
        assert.deepEqual(actualPairs(result), expectedPairs(loc, m));
      }),
      RUNS,
    );
  });

  it("targets describe exactly the image units", () => {
    fc.assert(
      fc.property(arbMapping(["test:A"], ["test:X"]), arbLocation("test:A"), (m, loc) => {
        const result = mapLocation(loc, m, ctx);
        const fromTargets = result.targets.flatMap((t) => t.location.segments.flatMap((s) => traversal(s).map((u) => `${s.ref}|${u}`)));
        const fromPieces = result.pieces.flatMap((p) => traversal(p.target).map((u) => `${p.target.ref}|${u}`));
        assert.deepEqual(new Set(fromTargets), new Set(fromPieces));
      }),
      RUNS,
    );
  });

  it("a segment truncated at an end yields a fuzzy boundary", () => {
    // Internal gaps (e.g. a base inserted in a transcript) do not truncate: the target may stay one contiguous
    // interval (spec-core §5.4) and the gap is reported as unmapped. Only unmapped ends are truncations.
    fc.assert(
      fc.property(arbMapping(["test:A"], ["test:X"]), arbSegment("test:A"), (m, seg) => {
        const plain: Segment = { ref: seg.ref, start: seg.start, end: seg.end, strand: seg.strand };
        const result = mapLocation({ outer: "test:A", kind: "join", segments: [plain] }, m, ctx);
        const segs = result.targets.flatMap((t) => t.location.segments);
        const anyFuzzy = segs.some((s) => s.fuzzyLow || s.fuzzyHigh);
        const unmappedUnits = new Set((result.unmapped?.segments ?? []).flatMap(traversal));
        const endCut = unmappedUnits.has(plain.start) || unmappedUnits.has(plain.end - 1);
        if (segs.length > 0 && endCut) assert.ok(anyFuzzy, "expected a fuzzy boundary");
        if (!result.unmapped) assert.ok(!anyFuzzy, "no truncation, no fuzzy boundary");
      }),
      RUNS,
    );
  });

  it("an internal gap keeps a contiguous target without truncation marks (regression: seed -346057382)", () => {
    const m = new Mapping([
      { srcRef: "test:A", src: 0, tgtRef: "test:X", tgt: 23, len: 6, rev: false },
      { srcRef: "test:A", src: 7, tgtRef: "test:X", tgt: 29, len: 1, rev: false },
    ]);
    const r = mapLocation({ outer: "test:A", kind: "join", segments: [{ ref: "test:A", start: 0, end: 8, strand: 1 }] }, m, ctx);
    assert.deepEqual(r.targets.map((t) => formatLocationId(t.location, ctx)), ["test:X:24..30"]);
    assert.equal(r.unmapped && formatLocationId(r.unmapped, ctx), "test:A:7");
  });
});

describe("merging is per target sequence", () => {
  it("adding blocks to another target does not change a target's result", () => {
    fc.assert(
      fc.property(arbMapping(["test:A"], ["test:X"]), arbMapping(["test:A"], ["test:Y"]), arbLocation("test:A"), (mx, my, loc) => {
        const alone = mapLocation(loc, mx, ctx).targets.map((t) => formatLocationId(t.location, ctx));
        const both = mapLocation(loc, mx.concat(my), ctx)
          .targets.filter((t) => t.location.outer === "test:X")
          .map((t) => formatLocationId(t.location, ctx));
        assert.deepEqual(both, alone);
      }),
      RUNS,
    );
  });
});

describe("mapping algebra", () => {
  it("invert is an involution", () => {
    fc.assert(
      fc.property(arbMapping(["test:A", "test:B"], ["test:X"]), (m) => {
        assert.deepEqual(invert(invert(m)).blocks, m.blocks);
      }),
      RUNS,
    );
  });

  it("invert swaps the unit relation", () => {
    fc.assert(
      fc.property(arbMapping(["test:A"], ["test:X"]), fc.nat(60), (m, u) => {
        const inv = invert(m);
        for (const img of mapUnit(m, "test:A", u)) {
          assert.ok(mapUnit(inv, img.ref, img.pos).some((back) => back.ref === "test:A" && back.pos === u));
        }
      }),
      RUNS,
    );
  });

  it("compose equals unit-wise composition", () => {
    fc.assert(
      fc.property(arbMapping(["test:A"], ["test:B"]), arbMapping(["test:B"], ["test:C"]), fc.nat(60), (ab, bc, u) => {
        const expected = mapUnit(ab, "test:A", u).flatMap((y) =>
          mapUnit(bc, y.ref, y.pos).map((z) => ({ ref: z.ref, pos: z.pos, rev: y.rev !== z.rev })),
        );
        assert.deepEqual(sortImages(mapUnit(compose(ab, bc), "test:A", u)), sortImages(expected));
      }),
      RUNS,
    );
  });

  it("compose is associative (as a unit relation)", () => {
    fc.assert(
      fc.property(
        arbMapping(["test:A"], ["test:B"]),
        arbMapping(["test:B"], ["test:C"]),
        arbMapping(["test:C"], ["test:D"]),
        fc.nat(60),
        (ab, bc, cd, u) => {
          const left = compose(compose(ab, bc), cd);
          const right = compose(ab, compose(bc, cd));
          assert.deepEqual(sortImages(mapUnit(left, "test:A", u)), sortImages(mapUnit(right, "test:A", u)));
        },
      ),
      RUNS,
    );
  });

  it("round trip through a mapping contains the mapped part of the input", () => {
    fc.assert(
      fc.property(arbMapping(["test:A"], ["test:X"]), arbSegment("test:A"), (m, seg) => {
        const loc: Location = { outer: "test:A", kind: "join", segments: [seg] };
        const forward = mapLocation(loc, m, ctx);
        const mappedUnits = new Set(forward.pieces.flatMap((p) => traversal(p.source)));
        const inv = invert(m);
        const back = new Set<number>();
        for (const t of forward.targets) {
          for (const p of mapLocation(t.location, inv, ctx).pieces) {
            if (p.target.ref === "test:A") traversal(p.target).forEach((u) => back.add(u));
          }
        }
        for (const u of mappedUnits) assert.ok(back.has(u), `unit ${u} lost in round trip`);
      }),
      RUNS,
    );
  });
});

describe("canonical text", () => {
  const roundTrip = (loc: Location) => {
    const text = formatLocationId(loc, ctx);
    const reparsed = parseLocationId(text, ctx);
    return { text, reparsed, again: formatLocationId(reparsed, ctx) };
  };

  it("nt locations: format -> parse is the identity and formatting is idempotent", () => {
    fc.assert(
      fc.property(arbLocation("test:A"), (loc) => {
        const { text, reparsed, again } = roundTrip(loc);
        assert.deepEqual(reparsed.segments, loc.segments, text);
        assert.equal(again, text);
      }),
      RUNS,
    );
  });

  it("protein locations with codon-level boundaries round-trip", () => {
    fc.assert(
      fc.property(arbLocation("test:P", false), (loc) => {
        const { text, reparsed, again } = roundTrip(loc);
        assert.deepEqual(reparsed.segments, loc.segments, text);
        assert.equal(again, text);
      }),
      RUNS,
    );
  });

  it("codon=never covers whole residues", () => {
    fc.assert(
      fc.property(arbLocation("test:P", false), (loc) => {
        const covered = parseLocationId(formatLocationId(loc, ctx, "never"), ctx);
        loc.segments.forEach((s, i) => {
          const c = covered.segments[i]!;
          assert.equal(c.start, s.start - (s.start % 3));
          assert.equal(c.end, s.end + ((3 - (s.end % 3)) % 3));
        });
      }),
      RUNS,
    );
  });
});
