// Cross-species check of liftOver chains with orthologous coding sequences.
//   node service/bench/verify-liftover.ts MANE_SUMMARY HUMAN.fna MOUSE.fna SAMPLES SEED STORE.sqlite...
// For random human MANE Select proteins and residues: human protein -> mouse genome (through the human genome and a
// chain) and -> mouse protein. Orthologous codons translate to the same amino acid most of the time (~85-90%), while
// misplaced coordinates give ~5%; the mouse protein reached should be encoded by the gene of the same symbol.
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { parseLocationId } from "@togocoord/core";
import { accessionRef, extract, FaiSequenceSource, translate } from "@togocoord/ingest";
import { convert, DEFAULT_PREFER, StoreSet } from "../src/index.ts";

const [summary, humanFna, mouseFna, samplesArg, seedArg, ...paths] = process.argv.slice(2);
if (!summary || !humanFna || !mouseFna || !samplesArg || !seedArg || paths.length === 0) {
  process.stderr.write("usage: verify-liftover.ts MANE_SUMMARY HUMAN.fna MOUSE.fna SAMPLES SEED STORE.sqlite...\n");
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
const human = new FaiSequenceSource(humanFna, toRef);
const mouse = new FaiSequenceSource(mouseFna, toRef);
const isMouse = (ref: string) => mouse.length(ref) !== undefined;
const MOUSE = 10090;

let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);
const t = { samples: 0, mouseGenome: 0, sameAminoAcid: 0, mouseProtein: 0, sameGeneSymbol: 0, otherGene: 0, noHumanGenome: 0 };
const ms: number[] = [];
const examples: string[] = [];

for (let i = 0; i < Number(samplesArg); i++) {
  const g = genes[Math.floor(random() * genes.length)]!;
  const len = stores.sequence(g.np)?.length;
  if (!len) continue;
  const r = 1 + Math.floor(random() * len);
  t.samples++;
  const humanCodon = convert(stores, parseLocationId(`${g.np}:${r}`, ctx), { to: { category: "genome" }, prefer: DEFAULT_PREFER }, ctx).find(
    (h) => /:NC_/.test(h.location.outer) && !isMouse(h.location.outer),
  );
  if (!humanCodon) {
    t.noHumanGenome++;
    continue;
  }
  const aaHuman = translate(extract(humanCodon.location, human) ?? "");
  const start = performance.now();
  const lifted = convert(stores, humanCodon.location, { to: { category: "genome" }, taxon: MOUSE, prefer: DEFAULT_PREFER }, ctx).find((h) => isMouse(h.location.outer));
  const proteins = convert(stores, parseLocationId(`${g.np}:${r}`, ctx), { to: { category: "protein" }, taxon: MOUSE, prefer: DEFAULT_PREFER }, ctx);
  ms.push(performance.now() - start);
  if (lifted) {
    t.mouseGenome++;
    const codon = extract(lifted.location, mouse);
    if (codon && codon.length === 3 && translate(codon) === aaHuman) t.sameAminoAcid++;
  }
  if (proteins.length) {
    t.mouseProtein++;
    // The gene of the reached mouse protein, from its CDS edge (GFF3 `gene` attribute).
    const symbols = new Set(proteins.flatMap((p) => stores.edges(p.location.outer).map((e) => e.attributes.gene?.toUpperCase()).filter(Boolean)));
    if (symbols.has(g.symbol.toUpperCase())) t.sameGeneSymbol++;
    else {
      t.otherGene++;
      if (examples.length < 12) examples.push(`${g.symbol} ${g.np}:${r} -> ${proteins[0]!.id} (${[...symbols].join(",")}) via ${proteins[0]!.path.map((s) => s.kind).join(">")}`);
    }
  }
}

const pct = (p: number) => [...ms].sort((a, b) => a - b)[Math.min(ms.length - 1, Math.floor((p / 100) * ms.length))]?.toFixed(1);
const rate = (n: number, d: number) => `${n}/${d} (${d ? ((100 * n) / d).toFixed(1) : "-"}%)`;
console.log(JSON.stringify(t));
console.log(`human codon lifted to mouse genome (taxon=10090): ${rate(t.mouseGenome, t.samples - t.noHumanGenome)}`);
console.log(`  same amino acid at the lifted codon: ${rate(t.sameAminoAcid, t.mouseGenome)}`);
console.log(`human residue reaches a mouse protein: ${rate(t.mouseProtein, t.samples - t.noHumanGenome)}`);
console.log(`  of the same gene symbol: ${rate(t.sameGeneSymbol, t.mouseProtein)}`);
console.log(`latency (lift + protein): p50 ${pct(50)} ms, p95 ${pct(95)} ms, p99 ${pct(99)} ms`);
for (const e of examples) console.log(`  ${e}`);
