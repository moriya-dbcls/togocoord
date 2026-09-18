// GRCh38 <-> GRCh37 through UCSC liftOver chains (spec-service §2.2).
//   node service/bench/verify-assembly.ts MANE_SUMMARY GRCH38.fna GRCH37.fna SAMPLES SEED STORE.sqlite...
// For random human MANE Select proteins and residues: protein -> genome on GRCh38, and on GRCh37 (through GRCh38 and a
// chain). The two codons must encode the same amino acid (the assemblies rarely differ in coding sequence), and the
// GRCh37 location must convert back to the same residue (GRCh37 -> chain -> GRCh38 -> CDS).
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { parseLocationId } from "@togocoord/core";
import { accessionRef, extract, FaiSequenceSource, translate } from "@togocoord/ingest";
import { convert, DEFAULT_PREFER, StoreSet } from "../src/index.ts";

const [summary, grch38Fna, grch37Fna, samplesArg, seedArg, ...paths] = process.argv.slice(2);
if (!summary || !grch38Fna || !grch37Fna || !samplesArg || !seedArg || paths.length === 0) {
  process.stderr.write("usage: verify-assembly.ts MANE_SUMMARY GRCH38.fna GRCH37.fna SAMPLES SEED STORE.sqlite...\n");
  process.exit(1);
}
const bytes = readFileSync(summary);
const genes = (summary.endsWith(".gz") ? gunzipSync(bytes) : bytes)
  .toString("utf8")
  .split("\n")
  .filter((l) => l && !l.startsWith("#"))
  .map((l) => l.split("\t"))
  .filter((c) => c[9] === "MANE Select")
  .map((c) => ({ symbol: c[3]!, np: `refseq:${c[6]}` }));

const stores = new StoreSet(paths);
const ctx = stores.context();
const toRef = (id: string) => accessionRef(id, stores.registry);
const grch38 = new FaiSequenceSource(grch38Fna, toRef);
const grch37 = new FaiSequenceSource(grch37Fna, toRef);
const GRCH37 = stores.assembly("GRCh37")?.name;
if (!GRCH37) throw new Error("no GRCh37 assembly among the stores");

let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);
const t = { samples: 0, grch37: 0, sameAminoAcid: 0, roundTrip: 0, notLifted: 0 };
const ms: number[] = [];
const examples: string[] = [];

for (let i = 0; i < Number(samplesArg); i++) {
  const g = genes[Math.floor(random() * genes.length)]!;
  const len = stores.sequence(g.np)?.length;
  if (!len) continue;
  const r = 1 + Math.floor(random() * len);
  t.samples++;
  const residue = parseLocationId(`${g.np}:${r}`, ctx);
  const on38 = convert(stores, residue, { to: { category: "genome" }, prefer: DEFAULT_PREFER }, ctx)[0];
  const start = performance.now();
  const on37 = convert(stores, residue, { to: { category: "genome" }, assembly: GRCH37, prefer: DEFAULT_PREFER }, ctx)[0];
  ms.push(performance.now() - start);
  if (!on37 || !on38) {
    t.notLifted++;
    if (examples.length < 10) examples.push(`not lifted: ${g.symbol} ${g.np}:${r}`);
    continue;
  }
  t.grch37++;
  const aa38 = translate(extract(on38.location, grch38) ?? "");
  const aa37 = translate(extract(on37.location, grch37) ?? "");
  if (aa38 && aa38 === aa37) t.sameAminoAcid++;
  else if (examples.length < 10) examples.push(`amino acid: ${g.symbol} ${g.np}:${r} ${on38.id} ${aa38} / ${on37.id} ${aa37}`);
  const back = convert(stores, on37.location, { to: { ref: g.np }, prefer: DEFAULT_PREFER }, ctx)[0];
  if (back?.id === `${g.np}:${r}`) t.roundTrip++;
  else if (examples.length < 10) examples.push(`round trip: ${g.np}:${r} -> ${on37.id} -> ${back?.id ?? "none"}`);
}

const pct = (p: number) => [...ms].sort((a, b) => a - b)[Math.min(ms.length - 1, Math.floor((p / 100) * ms.length))]?.toFixed(1);
const rate = (n: number, d: number) => `${n}/${d} (${d ? ((100 * n) / d).toFixed(1) : "-"}%)`;
console.log(JSON.stringify(t));
console.log(`protein residue -> GRCh37 codon: ${rate(t.grch37, t.samples)}`);
console.log(`  same amino acid as the GRCh38 codon: ${rate(t.sameAminoAcid, t.grch37)}`);
console.log(`  GRCh37 codon -> the same residue: ${rate(t.roundTrip, t.grch37)}`);
console.log(`latency (protein -> GRCh37): p50 ${pct(50)} ms, p95 ${pct(95)} ms, p99 ${pct(99)} ms`);
for (const e of examples) console.log(`  ${e}`);
