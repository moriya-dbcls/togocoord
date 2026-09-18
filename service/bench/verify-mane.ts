// MANE preference check on real data.
//   node service/bench/verify-mane.ts MANE_SUMMARY SAMPLES SEED STORE.sqlite...
// For random MANE genes: residue r of the MANE Select RefSeq protein -> genome; then genome -> protein with MANE
// preference. The first result should be a MANE Select protein (RefSeq NP or Ensembl ENSP of the same gene) at the
// same residue; the RefSeq and Ensembl MANE proteins must give the same genomic codon.
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { formatLocationId, parseLocationId } from "@togocoord/core";
import { convert, DEFAULT_PREFER, StoreSet } from "../src/index.ts";

const [summary, samplesArg, seedArg, ...paths] = process.argv.slice(2);
if (!summary || !samplesArg || !seedArg || paths.length === 0) {
  process.stderr.write("usage: verify-mane.ts MANE_SUMMARY SAMPLES SEED STORE.sqlite...\n");
  process.exit(1);
}
const bytes = readFileSync(summary);
const rows = (summary.endsWith(".gz") ? gunzipSync(bytes) : bytes).toString("utf8").split("\n").filter((l) => l && !l.startsWith("#")).map((l) => l.split("\t"));
const select = rows.filter((c) => c[9] === "MANE Select").map((c) => ({ gene: c[3]!, np: `refseq:${c[6]}`, ensp: `ensembl:${c[8]}` }));

const stores = new StoreSet(paths);
const ctx = stores.context();
let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);
const t = { samples: 0, noLength: 0, npNoGenome: 0, sameCodonNpEnsp: 0, differentCodon: 0, firstIsMane: 0, firstNotMane: 0, sameResidueBack: 0 };
const examples: string[] = [];
for (let i = 0; i < Number(samplesArg); i++) {
  const g = select[Math.floor(random() * select.length)]!;
  const len = stores.sequence(g.np)?.length;
  if (!len) {
    t.noLength++;
    continue;
  }
  t.samples++;
  const r = 1 + Math.floor(random() * len);
  const toGenome = (ref: string) =>
    convert(stores, parseLocationId(`${ref}:${r}`, ctx), { to: { category: "genome" }, prefer: DEFAULT_PREFER }, ctx).find((h) => /:NC_/.test(h.location.outer));
  const a = toGenome(g.np);
  if (!a) {
    t.npNoGenome++;
    continue;
  }
  const b = toGenome(g.ensp);
  if (b && b.id === a.id) t.sameCodonNpEnsp++;
  else {
    t.differentCodon++;
    if (examples.length < 10) examples.push(`${g.gene} ${g.np}:${r} -> ${a.id} but ${g.ensp}:${r} -> ${b?.id ?? "(none)"}`);
  }
  const back = convert(stores, a.location, { to: { category: "protein" }, prefer: DEFAULT_PREFER }, ctx);
  const first = back[0];
  if (first && (first.location.outer === g.np || first.location.outer === g.ensp)) {
    t.firstIsMane++;
    if (formatLocationId(first.location, ctx, "never").endsWith(`:${r}`)) t.sameResidueBack++;
  } else {
    t.firstNotMane++;
    if (examples.length < 20) examples.push(`${g.gene}: ${a.id} -> first ${first?.id ?? "(none)"} [${first?.tags.join(",") ?? ""}], expected ${g.np} or ${g.ensp}`);
  }
}
console.log(JSON.stringify(t));
for (const e of examples) console.log(`  ${e}`);
