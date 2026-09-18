// UniProt identity coverage and genome <-> PDB structure verification.
//   node service/bench/verify-structure.ts GENOME.fna PDB_SEQRES UNIPROT_FASTA[,UNIPROT_FASTA] UNIPROT.sqlite SIFTS.sqlite SAMPLES SEED STORE.sqlite...
// 1. For every UniProt sequence: is there an identical (same refget digest) Ensembl / RefSeq protein?
// 2. For random SIFTS residues: UniProt residue -> genome through the cheapest path; the codon must encode the UniProt
//    residue, and the genome location must convert back to the same PDB residue.
import { readFileSync } from "node:fs";
import { DatabaseSync } from "node:sqlite";
import { gunzipSync } from "node:zlib";
import { formatLocationId, parseLocationId } from "@togocoord/core";
import { accessionRef, defaultFastaRef, extract, FaiSequenceSource, parseFastaHeaders, translate } from "@togocoord/ingest";
import { convert, StoreSet } from "../src/index.ts";

const [fna, seqres, uniprotFasta, uniprotDb, siftsDb, samplesArg, seedArg, ...others] = process.argv.slice(2);
if (!fna || !seqres || !uniprotFasta || !uniprotDb || !siftsDb || !samplesArg || !seedArg) {
  process.stderr.write("usage: verify-structure.ts GENOME.fna PDB_SEQRES UNIPROT_FASTA[,UNIPROT_FASTA] UNIPROT.sqlite SIFTS.sqlite SAMPLES SEED STORE.sqlite...\n");
  process.exit(1);
}

// ---- 1. identity coverage -------------------------------------------------------------------------------------
const digests = (path: string, where: string) => {
  const db = new DatabaseSync(path, { readOnly: true });
  const rows = db.prepare(`SELECT ref, digest, provenance FROM sequence WHERE digest IS NOT NULL AND ${where}`).all() as Array<{ ref: string; digest: string; provenance: string }>;
  db.close();
  return rows;
};
const byNamespace = new Map<string, Set<string>>();
for (const path of others) {
  for (const r of digests(path, "moltype = 'protein'")) {
    const ns = r.ref.slice(0, r.ref.indexOf(":"));
    if (!byNamespace.has(ns)) byNamespace.set(ns, new Set());
    byNamespace.get(ns)!.add(r.digest);
  }
}
const ens = byNamespace.get("ensembl") ?? new Set<string>();
const rs = byNamespace.get("refseq") ?? new Set<string>();
const siftsAccessions = new Set(
  (new DatabaseSync(siftsDb, { readOnly: true }).prepare("SELECT DISTINCT s.ref FROM edge e JOIN sequence s ON s.id = e.from_seq").all() as Array<{ ref: string }>).map((r) => r.ref),
);
type Cov = { total: number; ensembl: number; refseq: number; either: number; none: number };
const cov: Record<string, Cov> = {};
for (const r of digests(uniprotDb, "ref LIKE 'uniprot:%'")) {
  const file = String(JSON.parse(r.provenance).file ?? "");
  const set = `${/additional/.test(file) ? "additional (isoforms etc.)" : "canonical proteome"}${siftsAccessions.has(r.ref) ? " / with structure" : ""}`;
  const c = (cov[set] ??= { total: 0, ensembl: 0, refseq: 0, either: 0, none: 0 });
  const e = ens.has(r.digest);
  const n = rs.has(r.digest);
  c.total++;
  if (e) c.ensembl++;
  if (n) c.refseq++;
  if (e || n) c.either++;
  else c.none++;
}
console.log("identity coverage (UniProt sequences with an identical Ensembl / RefSeq protein):");
for (const [k, c] of Object.entries(cov).sort()) {
  const p = (x: number) => `${x} (${((100 * x) / c.total).toFixed(1)}%)`;
  console.log(`  ${k.padEnd(42)} total ${c.total}  ensembl ${p(c.ensembl)}  refseq ${p(c.refseq)}  either ${p(c.either)}  none ${p(c.none)}`);
}

// ---- 2. structure round trip ----------------------------------------------------------------------------------
const stores = new StoreSet([uniprotDb, siftsDb, ...others]);
const ctx = stores.context();
const toRef = (id: string) => accessionRef(id, stores.registry);
const genome = new FaiSequenceSource(fna, toRef);
const pdbSeq = new FaiSequenceSource(seqres, (name) => {
  const m = /^([0-9][A-Za-z0-9]{3})_(\S+)$/.exec(name);
  try {
    return m ? stores.registry.refKey("pdb", `${m[1]}.${m[2]}`) : undefined;
  } catch {
    return undefined;
  }
});
const uniprotResidues = new Map<string, string>();
for (const file of uniprotFasta.split(",")) {
  const bytes = readFileSync(file);
  for (const [header, seq] of parseFastaHeaders((file.endsWith(".gz") ? gunzipSync(bytes) : bytes).toString("utf8"))) {
    const ref = defaultFastaRef(header, stores.registry);
    if (ref) uniprotResidues.set(ref, seq);
  }
}
const sdb = new DatabaseSync(siftsDb, { readOnly: true });
const blocks = sdb.prepare(
  "SELECT u.ref AS uniprot, p.ref AS pdb, b.src, b.tgt, b.len FROM block b JOIN sequence u ON u.id = b.src_seq JOIN sequence p ON p.id = b.tgt_seq",
).all() as Array<{ uniprot: string; pdb: string; src: number; tgt: number; len: number }>;
sdb.close();

let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);
const t = {
  samples: 0,
  noGenome: 0,
  exact: 0,
  approximate: 0,
  codonMatchesUniprot: 0,
  codonMismatch: 0,
  pdbResidueMatchesUniprot: 0,
  roundTripSamePdbResidue: 0,
  roundTripLost: 0,
};
const ms: number[] = [];
const examples: string[] = [];
for (let i = 0; i < Number(samplesArg); i++) {
  const b = blocks[Math.floor(random() * blocks.length)]!;
  const k = Math.floor(random() * (b.len / 3));
  const residue = b.src / 3 + k + 1;
  const pdbResidue = b.tgt / 3 + k + 1;
  t.samples++;
  const start = performance.now();
  const toGenome = convert(stores, parseLocationId(`${b.uniprot}:${residue}`, ctx), { to: { category: "genome" } }, ctx).filter(
    (h) => genome.length(h.location.outer) !== undefined && /:NC_/.test(h.location.outer),
  );
  ms.push(performance.now() - start);
  const hit = toGenome.find((h) => !h.approximate) ?? toGenome[0];
  if (!hit) {
    t.noGenome++;
    continue;
  }
  if (hit.approximate) t.approximate++;
  else t.exact++;
  const codon = extract(hit.location, genome);
  const aa = codon && codon.length === 3 ? translate(codon) : undefined;
  const pdbAa = pdbSeq.get(b.pdb, pdbResidue - 1, pdbResidue);
  const expected = uniprotResidues.get(b.uniprot)?.[residue - 1];
  if (aa === expected || (aa === "*" && (expected === "U" || expected === "O"))) t.codonMatchesUniprot++;
  else {
    t.codonMismatch++;
    if (examples.length < 10) examples.push(`${b.uniprot}:${residue} -> ${hit.id} ${codon} (${aa}) expected ${expected} via ${hit.path.map((s) => s.kind).join(">")}${hit.approximate ? " [approximate]" : ""}`);
  }
  if (pdbAa === expected) t.pdbResidueMatchesUniprot++;
  const back = convert(stores, hit.location, { to: { ref: b.pdb } }, ctx);
  // A chain can contain the same UniProt residue twice (tandem constructs), so the residue must be among the targets.
  const wanted = pdbResidue;
  const backResidues = back.flatMap((x) => x.location.segments.map((sg) => [Math.floor(sg.start / 3) + 1, Math.floor((sg.end - 1) / 3) + 1]));
  if (backResidues.some(([lo, hi]) => lo! <= wanted && wanted <= hi!)) t.roundTripSamePdbResidue++;
  else {
    t.roundTripLost++;
    if (examples.length < 20) examples.push(`round trip ${b.uniprot}:${residue} -> ${hit.id} -> ${back.map((x) => x.id).join(" ") || "(none)"}; expected ${b.pdb}:${pdbResidue}`);
  }
}

const pct = (p: number) => [...ms].sort((a, b) => a - b)[Math.min(ms.length - 1, Math.floor((p / 100) * ms.length))]!.toFixed(2);
console.log("structure round trip:", JSON.stringify(t));
console.log(`UniProt residue -> genome: p50 ${pct(50)} ms, p95 ${pct(95)} ms, p99 ${pct(99)} ms`);
for (const e of examples) console.log(`  ${e}`);
