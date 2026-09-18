// Revised genes: proteins of a newer assembly with no identical protein anywhere (the gene model or its residues
// changed), so they reach an older protein and its UniProt entry only through the genome alignment:
//   newer residue -> codon on the newer genome -> (alignment) -> older genome -> older CDS -> identical UniProt entry.
//   node service/bench/verify-revised-genes.ts NEW.fna OLD.fna NEW_STORE.sqlite OLD_ASSEMBLY SAMPLES SEED STORE.sqlite...
// The UniProt residue reached must encode, on the older genome, the amino acid of the newer codon (unless the residue
// itself was revised).
import { DatabaseSync } from "node:sqlite";
import { parseLocationId } from "@togocoord/core";
import { accessionRef, extract, FaiSequenceSource, translate } from "@togocoord/ingest";
import { convert, StoreSet } from "../src/index.ts";

const [newFna, oldFna, newStore, oldAssembly, samplesArg, seedArg, ...paths] = process.argv.slice(2);
if (!newFna || !oldFna || !newStore || !oldAssembly || !samplesArg || !seedArg || paths.length === 0) {
  process.stderr.write("usage: verify-revised-genes.ts NEW.fna OLD.fna NEW_STORE.sqlite OLD_ASSEMBLY SAMPLES SEED STORE.sqlite...\n");
  process.exit(1);
}
const stores = new StoreSet(paths);
const ctx = stores.context();
const toRef = (id: string) => accessionRef(id, stores.registry);
const newer = new FaiSequenceSource(newFna, toRef);
const older = new FaiSequenceSource(oldFna, toRef);
const db = new DatabaseSync(newStore, { readOnly: true });
const revised = (db.prepare("SELECT ref, length FROM sequence WHERE moltype = 'protein' AND length > 50").all() as Array<{ ref: string; length: number }>).filter(
  (p) => stores.identical(p.ref).length === 0,
);

let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);
const t = { revisedProteins: revised.length, samples: 0, onNewGenome: 0, uniprotViaAlignment: 0, sameAminoAcid: 0, uniprotOtherGeneOnly: 0 };
const good: string[] = [];
const differs: string[] = [];

for (let i = 0; i < Number(samplesArg); i++) {
  const p = revised[Math.floor(random() * revised.length)]!;
  const residue = parseLocationId(`${p.ref}:${1 + Math.floor(random() * p.length)}`, ctx);
  t.samples++;
  const codon = convert(stores, residue, { to: { category: "genome" } }, ctx)[0];
  if (!codon) continue;
  t.onNewGenome++;
  const hits = convert(stores, codon.location, { to: { category: "protein" }, prefer: [] }, ctx).filter((h) => h.location.outer.startsWith("uniprot:"));
  // Through the alignment: the same gene in the older assembly. Without it: another gene overlapping the codon on the
  // newer genome (another frame or strand) whose protein happens to be unchanged.
  const hit = hits.find((h) => h.path.some((s) => s.kind === "liftover"));
  if (!hit) {
    if (hits.length) t.uniprotOtherGeneOnly++;
    continue;
  }
  t.uniprotViaAlignment++;
  const aaNew = translate(extract(codon.location, newer) ?? "");
  const old = convert(stores, hit.location, { to: { category: "genome" }, assembly: oldAssembly }, ctx)[0];
  const aaOld = old ? translate(extract(old.location, older) ?? "") : "";
  const line = `${p.ref}:${residue.segments[0]!.start / 3 + 1} (${aaNew}) -> ${codon.id} -> ${hit.id} (${aaOld || "?"})`;
  if (aaOld && aaOld === aaNew) {
    t.sameAminoAcid++;
    if (good.length < 5) good.push(line);
  } else if (differs.length < 8) differs.push(line);
}

const rate = (n: number, d: number) => `${n}/${d} (${d ? ((100 * n) / d).toFixed(1) : "-"}%)`;
console.log(JSON.stringify(t));
console.log(`revised proteins (no identical protein anywhere): ${t.revisedProteins}`);
console.log(`residues reaching a UniProt entry through the alignment: ${rate(t.uniprotViaAlignment, t.onNewGenome)}`);
console.log(`  same amino acid: ${rate(t.sameAminoAcid, t.uniprotViaAlignment)}`);
console.log(`residues reaching only another gene's UniProt entry: ${t.uniprotOtherGeneOnly}`);
for (const g of good) console.log(`  ok: ${g}`);
for (const d of differs) console.log(`  differs: ${d}`);
