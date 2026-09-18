// Two assemblies of one species joined by a genome alignment (PAF / chain), e.g. Marchantia MpTak v3.1 -> v7.1.
//   node service/bench/verify-genome-pair.ts FROM.fna TO.fna FROM_ASSEMBLY TO_ASSEMBLY FROM_STORE.sqlite SAMPLES SEED STORE.sqlite...
// 1. Random 21-base windows of the FROM genome -> TO genome (alignment only): lifted, and identical residues.
// 2. Random residues of FROM proteins with an identical protein on TO: the codon lifted by the alignment alone must be
//    the codon reached through the identical protein (genome -> CDS -> identical protein -> CDS -> genome).
import { DatabaseSync } from "node:sqlite";
import { parseLocationId, formatLocationId } from "@togocoord/core";
import { accessionRef, extract, FaiSequenceSource } from "@togocoord/ingest";
import { convert, StoreSet } from "../src/index.ts";

const [fromFna, toFna, fromAssembly, toAssembly, fromStore, samplesArg, seedArg, ...paths] = process.argv.slice(2);
if (!fromFna || !toFna || !fromAssembly || !toAssembly || !fromStore || !samplesArg || !seedArg || paths.length === 0) {
  process.stderr.write("usage: verify-genome-pair.ts FROM.fna TO.fna FROM_ASSEMBLY TO_ASSEMBLY FROM_STORE.sqlite SAMPLES SEED STORE.sqlite...\n");
  process.exit(1);
}
const stores = new StoreSet(paths);
const ctx = stores.context();
const toRef = (id: string) => accessionRef(id, stores.registry);
const fromSeq = new FaiSequenceSource(fromFna, toRef);
const toSeq = new FaiSequenceSource(toFna, toRef);
const samples = Number(samplesArg);
let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);
const pct = (xs: number[], p: number) => [...xs].sort((a, b) => a - b)[Math.min(xs.length - 1, Math.floor((p / 100) * xs.length))]?.toFixed(1);
const rate = (n: number, d: number) => `${n}/${d} (${d ? ((100 * n) / d).toFixed(1) : "-"}%)`;

// 1. Genome windows, weighted by sequence length.
const genomes = stores.assemblies().find((a) => a.name === fromAssembly)!;
const refs = [...genomes.refs].map((r) => ({ ref: r, len: fromSeq.length(r) ?? 0 })).filter((x) => x.len > 1000);
const total = refs.reduce((n, x) => n + x.len, 0);
const w = { samples: 0, lifted: 0, identical: 0, split: 0 };
const ms: number[] = [];
for (let i = 0; i < samples; i++) {
  let at = Math.floor(random() * total);
  const g = refs.find((x) => (at -= x.len) < 0)!;
  const start = Math.floor(random() * (g.len - 21));
  const loc = parseLocationId(`${g.ref}:${start + 1}..${start + 21}`, ctx);
  w.samples++;
  const t0 = performance.now();
  const hit = convert(stores, loc, { to: { category: "genome" }, assembly: toAssembly, maxHops: 1 }, ctx)[0];
  ms.push(performance.now() - t0);
  if (!hit) continue;
  w.lifted++;
  if (hit.location.segments.length !== 1) {
    w.split++;
    continue;
  }
  if (extract(hit.location, toSeq) === extract(loc, fromSeq)) w.identical++;
}

// 2. Codons of identical proteins: alignment vs identical protein.
const db = new DatabaseSync(fromStore, { readOnly: true });
const proteins = (db.prepare("SELECT ref, length FROM sequence WHERE moltype = 'protein' AND length > 0").all() as Array<{ ref: string; length: number }>).filter(
  (p) => stores.identical(p.ref).some((r) => stores.category(r) === "protein" && r !== p.ref && !r.startsWith("uniprot:")),
);
const c = { samples: 0, bothReached: 0, same: 0 };
const examples: string[] = [];
for (let i = 0; i < samples && proteins.length; i++) {
  const p = proteins[Math.floor(random() * proteins.length)]!;
  const residue = parseLocationId(`${p.ref}:${1 + Math.floor(random() * p.length)}`, ctx);
  const onFrom = convert(stores, residue, { to: { category: "genome" }, assembly: fromAssembly }, ctx)[0];
  const onTo = convert(stores, residue, { to: { category: "genome" }, assembly: toAssembly }, ctx).find((h) => h.path.every((s) => s.kind !== "liftover"));
  if (!onFrom || !onTo) continue;
  c.samples++;
  const lifted = convert(stores, onFrom.location, { to: { category: "genome" }, assembly: toAssembly, maxHops: 1 }, ctx)[0];
  if (!lifted) continue;
  c.bothReached++;
  if (lifted.id === onTo.id) c.same++;
  else if (examples.length < 8) examples.push(`${formatLocationId(residue, ctx)}: ${onFrom.id} -> alignment ${lifted.id} / protein ${onTo.id}`);
}

console.log(JSON.stringify({ windows: w, codons: c }));
console.log(`21-base windows of ${fromAssembly} lifted to ${toAssembly}: ${rate(w.lifted, w.samples)}`);
console.log(`  identical residues: ${rate(w.identical, w.lifted - w.split)} (split by an indel: ${w.split})`);
console.log(`codons of identical proteins lifted by the alignment: ${rate(c.bothReached, c.samples)}`);
console.log(`  same codon as through the identical protein: ${rate(c.same, c.bothReached)}`);
console.log(`latency (window lift): p50 ${pct(ms, 50)} ms, p99 ${pct(ms, 99)} ms`);
for (const e of examples) console.log(`  ${e}`);
