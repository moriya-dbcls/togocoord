// Verify path search against published proteins.
//   node service/bench/verify-search.ts GENOME.fna PROTEINS.faa[.gz] SAMPLES SEED [--exception-only] STORE.sqlite...
// For random (protein, residue) pairs: convert to the genome through the cheapest path, translate the genomic codon
// and compare with the published residue. Results are grouped by the kind of path taken.
import { readFileSync } from "node:fs";
import { DatabaseSync } from "node:sqlite";
import { gunzipSync } from "node:zlib";
import { formatLocationId, parseLocationId } from "@togocoord/core";
import { accessionRef, extract, FaiSequenceSource, parseFasta, translate } from "@togocoord/ingest";
import { convert, StoreSet } from "../src/index.ts";

const args = process.argv.slice(2);
const exceptionOnly = args.includes("--exception-only");
const [fna, faa, samplesArg, seedArg, ...paths] = args.filter((a) => a !== "--exception-only");
if (!fna || !faa || !samplesArg || !seedArg || paths.length === 0) {
  process.stderr.write("usage: verify-search.ts GENOME.fna PROTEINS.faa[.gz] SAMPLES SEED [--exception-only] STORE.sqlite...\n");
  process.exit(1);
}

const stores = new StoreSet(paths);
const ctx = stores.context();
const toRef = (id: string) => accessionRef(id, stores.registry);
const proteins = new Map<string, string>();
const faaBytes = readFileSync(faa);
for (const [id, seq] of parseFasta((faa.endsWith(".gz") ? gunzipSync(faaBytes) : faaBytes).toString("utf8"))) {
  const ref = toRef(id);
  if (ref) proteins.set(ref, seq);
}
const genome = new FaiSequenceSource(fna, toRef);

// Proteins with a genome CDS edge (optionally only those NCBI flags with /exception).
const candidates = new Set<string>();
for (const path of paths) {
  const db = new DatabaseSync(path, { readOnly: true });
  const sql = `SELECT DISTINCT s.ref FROM edge e JOIN sequence s ON s.id = e.from_seq
    WHERE s.moltype = 'protein' ${exceptionOnly ? `AND e.attributes LIKE '%"exception"%'` : ""}`;
  for (const { ref } of db.prepare(sql).all() as Array<{ ref: string }>) if (proteins.has(ref)) candidates.add(ref);
  db.close();
}
const pool = [...candidates].sort();

let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);
type Kind = "direct" | "via transcript" | "exception edge";
/** Path kind; an /exception flag counts only when the edge was not verified against the actual sequence. */
const kindOf = (path: Array<{ attributes: Record<string, string>; validation: { status: string; basis?: string } }>): Kind =>
  path.some((s) => s.attributes.exception !== undefined && !(s.validation.status === "ok" && s.validation.basis === "full"))
    ? "exception edge"
    : path.length > 1
      ? "via transcript"
      : "direct";
// agree/disagree: the codon at the mapped position encodes the published residue. Substitutions between genome and
// RefSeq sequences also cause disagreement, so coordinates are judged separately on a 15-residue window: a colinear
// mapping translates to >= 80% identity with the published window, a shifted one to near-random identity.
type Counts = { agree: number; disagree: number; windowColinear: number; windowShifted: number };
const zero = (): Counts => ({ agree: 0, disagree: 0, windowColinear: 0, windowShifted: 0 });
// Primary assembly (NC_) and alternate loci / patches (NT_, NW_) are counted separately: alternate haplotypes such
// as the MHC region differ from the reference proteins in many residues without any coordinate error.
const tally: Record<string, Counts> = {};
const bucket = (kind: Kind, genomeRef: string) => (tally[`${kind} / ${/:NC_/.test(genomeRef) ? "primary" : "alternate"}`] ??= zero());
let residues = 0;
let noResult = 0;
const ms: number[] = [];
const examples: string[] = [];

for (let i = 0; i < Number(samplesArg); i++) {
  const protein = pool[Math.floor(random() * pool.length)]!;
  const seq = proteins.get(protein)!;
  const r = 1 + Math.floor(random() * seq.length);
  residues++;
  const t = performance.now();
  const results = convert(stores, parseLocationId(`${protein}:${r}`, ctx), { to: { category: "genome" } }, ctx);
  ms.push(performance.now() - t);
  const hits = results.filter((h) => genome.length(h.location.outer) !== undefined);
  if (hits.length === 0) noResult++;
  for (const hit of hits) {
    const kind = kindOf(hit.path);
    const codon = extract(hit.location, genome);
    const table = Number(hit.path[0]!.attributes.translTable ?? 1);
    const aa = codon && codon.length === 3 ? translate(codon, table) : undefined;
    const expected = seq[r - 1]!;
    if (aa === expected || (aa === "*" && (expected === "U" || expected === "O")) || expected === "X") bucket(kind, hit.location.outer).agree++;
    else bucket(kind, hit.location.outer).disagree++;
  }

  // Window check (same path search, residues r-7..r+7, genome targets reached with the same path kind).
  if (r > 7 && r + 7 <= seq.length) {
    const window = convert(stores, parseLocationId(`${protein}:${r - 7}..${r + 7}`, ctx), { to: { category: "genome" } }, ctx);
    for (const hit of window.filter((h) => genome.length(h.location.outer) !== undefined)) {
      const kind = kindOf(hit.path);
      const nt = extract(hit.location, genome);
      if (!nt || nt.length !== 45) continue;
      const aa = translate(nt, Number(hit.path[0]!.attributes.translTable ?? 1));
      const want = seq.slice(r - 8, r + 7);
      let same = 0;
      for (let k = 0; k < 15; k++) if (aa[k] === want[k] || want[k] === "X") same++;
      if (same >= 12) bucket(kind, hit.location.outer).windowColinear++;
      else {
        bucket(kind, hit.location.outer).windowShifted++;
        if (kind !== "exception edge" && /:NC_/.test(hit.location.outer) && examples.length < 15) {
          examples.push(`${protein}:${r - 7}..${r + 7} -> ${formatLocationId(hit.location, ctx)} ${aa} vs ${want} via ${hit.path.map((s) => `${s.kind}:${s.to}[${s.validation.status}/${s.validation.basis ?? "-"}]`).join(" > ")}`);
        }
      }
    }
  }
}

const pct = (p: number) => [...ms].sort((a, b) => a - b)[Math.min(ms.length - 1, Math.floor((p / 100) * ms.length))]!.toFixed(2);
console.log(JSON.stringify({ pool: pool.length, residues, noResult }));
for (const [k, v] of Object.entries(tally).sort()) console.log(`${k.padEnd(28)} ${JSON.stringify(v)}`);
console.log(`convert to genome: p50 ${pct(50)} ms, p95 ${pct(95)} ms, p99 ${pct(99)} ms`);
for (const e of examples) console.log(`disagree: ${e}`);
