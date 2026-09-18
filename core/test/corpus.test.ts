// Runs the Phase 0 corpus (spec-core §7).
import assert from "node:assert/strict";
import { readdirSync, readFileSync } from "node:fs";
import { describe, it } from "node:test";
import {
  cdsMapping,
  formatLocationId,
  Mapping,
  mapLocation,
  mappingFromLocation,
  parseLocationId,
  residueBlock,
  type CodonMode,
  type ResidueBlockSpec,
  type Unit,
} from "../src/index.ts";
import { testContext } from "./helpers.ts";

type MappingSpec =
  | { kind: "cds"; protein: string; location: string; codonStart?: number; aaLength: number }
  | { kind: "location"; from: string; location: string }
  | { kind: "residueBlocks"; blocks: ResidueBlockSpec[] };

interface CorpusFile {
  id: string;
  title: string;
  units?: Record<string, Unit>;
  mappings: MappingSpec[];
  cases: Array<{
    note?: string;
    input: string;
    direction?: "forward" | "inverse";
    codon?: CodonMode;
    expect: { targets: string[]; unmapped: string | null; orientation?: string[] };
  }>;
}

const dir = new URL("./corpus/", import.meta.url);

for (const file of readdirSync(dir).filter((f) => f.endsWith(".json")).sort()) {
  const corpus = JSON.parse(readFileSync(new URL(file, dir), "utf8")) as CorpusFile;
  const ctx = testContext({ units: corpus.units ?? {} });
  const mapping = new Mapping(
    corpus.mappings.flatMap((spec) => {
      switch (spec.kind) {
        case "cds":
          return cdsMapping({
            protein: spec.protein,
            cds: parseLocationId(spec.location, ctx),
            codonStart: spec.codonStart ?? 1,
            aaLength: spec.aaLength,
          }).blocks;
        case "location":
          return mappingFromLocation(spec.from, parseLocationId(spec.location, ctx)).blocks;
        case "residueBlocks":
          return spec.blocks.map(residueBlock);
      }
    }),
  );
  const inverse = mapping.inverse();

  describe(`corpus ${corpus.id}: ${corpus.title}`, () => {
    for (const c of corpus.cases) {
      const dirLabel = c.direction === "inverse" ? "inverse" : "forward";
      it(`${dirLabel} ${c.input}${c.codon ? ` [codon=${c.codon}]` : ""}${c.note ? ` — ${c.note}` : ""}`, () => {
        const result = mapLocation(parseLocationId(c.input, ctx), c.direction === "inverse" ? inverse : mapping, ctx);
        assert.deepEqual(
          result.targets.map((t) => formatLocationId(t.location, ctx, c.codon)),
          c.expect.targets,
        );
        assert.equal(result.unmapped && formatLocationId(result.unmapped, ctx, c.codon), c.expect.unmapped);
        if (c.expect.orientation) {
          assert.deepEqual(result.targets.map((t) => t.orientation), c.expect.orientation);
        }
      });
    }
  });
}
