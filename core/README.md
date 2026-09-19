# @togocoord/core

English | [日本語](README.ja.md)

The TogoCoord core library (v0.1). It provides parsing and normalization of Location IDs, and coordinate mapping with block lists. It has no I/O.

- Design: [../docs/design.md](../docs/design.md)
- Rules: [../docs/spec-core.md](../docs/spec-core.md)

## Usage

```ts
import { cdsMapping, createContext, formatLocationId, mapLocation, parseLocationId } from "@togocoord/core";

const ctx = createContext();
// SARS-CoV-2 ORF1ab (NC_045512.2): CDS with a -1 ribosomal slippage
const toGenome = cdsMapping({
  protein: "refseq:YP_009724389.1",
  cds: parseLocationId("refseq:NC_045512.2:join(266..13468,13468..21555)", ctx),
  aaLength: 7096,
});

mapLocation(parseLocationId("refseq:YP_009724389.1:4401..4402", ctx), toGenome, ctx)
  .targets.map((t) => formatLocationId(t.location, ctx));
// => ["refseq:NC_045512.2:join(13466..13468,13468..13470)"]

mapLocation(parseLocationId("refseq:NC_045512.2:13468", ctx), toGenome.inverse(), ctx)
  .targets.map((t) => formatLocationId(t.location, ctx));
// => ["refseq:YP_009724389.1:join(4401c3,4402c1)"]   (join(4401,4402) in "never" mode)
```

## Layout

| File | Description |
|---|---|
| `src/parser.ts` | Location grammar (spec-core §3.2) |
| `src/namespace.ts` | Namespace registration, accession checks, default units (§3.4) |
| `src/location.ts` | Conversion to semantic values (segment lists), output in normal form (§3.3, §3.5) |
| `src/mapping.ts` | Blocks, `invert`, `compose`, `mappingFromLocation`, `cdsMapping` (§4) |
| `src/convert.ts` | `mapLocation` (§5) |
| `test/corpus/` | Phase 0 test data (4 real, 3 synthetic) |
| `test/oracle.ts` | Naive per-residue implementation (the expected answers for tests) |

## Development

Node.js 23.6 or later (TypeScript is run directly with type stripping only).

```sh
npm install
npm test          # node --test (test data + property tests + unit tests)
npm run typecheck # tsc --noEmit
```
