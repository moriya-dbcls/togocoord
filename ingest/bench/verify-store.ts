// Verify and benchmark a store against independent protein sequences.
//   node ingest/bench/verify-store.ts STORE.sqlite GENOME.fna PROTEINS.faa[.gz] [SAMPLES=2000] [SEED=1]
// For random (protein, residue) pairs: protein -> genome through the store, translate the genomic codon, compare with
// the published protein residue; then genome -> protein back (round trip). Reports agreement and query latency.
import { readFileSync } from "node:fs";
import { DatabaseSync } from "node:sqlite";
import { gunzipSync } from "node:zlib";
import { formatLocationId, NamespaceRegistry, parseLocationId } from "@togocoord/core";
import { accessionRef, extract, FaiSequenceSource, parseFasta, TogoCoordStore, translate } from "../src/index.ts";

const [dbPath, fnaPath, faaPath, samplesArg = "2000", seedArg = "1"] = process.argv.slice(2);
if (!dbPath || !fnaPath || !faaPath) {
  process.stderr.write("usage: verify-store.ts STORE.sqlite GENOME.fna PROTEINS.faa[.gz] [SAMPLES] [SEED]\n");
  process.exit(1);
}
const registry = new NamespaceRegistry();
const toRef = (id: string) => accessionRef(id, registry);

const faa = readFileSync(faaPath);
const proteins = new Map<string, string>();
for (const [id, seq] of parseFasta((faaPath.endsWith(".gz") ? gunzipSync(faa) : faa).toString("utf8"))) {
  const ref = toRef(id);
  if (ref) proteins.set(ref, seq);
}
const genome = new FaiSequenceSource(fnaPath, toRef);
const store = new TogoCoordStore(dbPath);
const ctx = store.context(registry);

// Proteins that have a CDS edge in the store.
const db = new DatabaseSync(dbPath, { readOnly: true });
const candidates = (db.prepare("SELECT ref FROM sequence WHERE moltype = 'protein' ORDER BY id").all() as Array<{ ref: string }>)
  .map((r) => r.ref)
  .filter((ref) => proteins.has(ref));
db.close();

let seed = Number(seedArg) >>> 0 || 1;
const random = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 2 ** 32);

const samples = Number(samplesArg);
const forwardMs: number[] = [];
const reverseMs: number[] = [];
// Edges carrying an INSDC /exception (e.g. "annotated by transcript or proteomic data") are genome models that NCBI
// itself marks as not reproducing the protein; they are counted separately.
const tally = {
  residues: 0,
  genomicHits: 0,
  agree: 0,
  disagree: 0,
  disagreeOnExceptionEdges: 0,
  noHit: 0,
  noHitOnExceptionEdges: 0,
  roundTrip: 0,
  roundTripMissing: 0,
  roundTripMissingOnExceptionEdges: 0,
};
const disagreements: string[] = [];

for (let i = 0; i < samples; i++) {
  const protein = candidates[Math.floor(random() * candidates.length)]!;
  const seq = proteins.get(protein)!;
  const r = 1 + Math.floor(random() * seq.length);
  const cdsEdges = store.edges(protein).filter((e) => e.from === protein && e.attributes.codonStart !== undefined);
  const edge = cdsEdges[0];
  if (!edge) continue;
  tally.residues++;

  let t = performance.now();
  const result = store.neighbors(parseLocationId(`${protein}:${r}`, ctx), ctx);
  forwardMs.push(performance.now() - t);
  const hits = result.targets.filter((x) => genome.length(x.location.outer) !== undefined);
  if (hits.length === 0) {
    // A genome model shorter than the published protein (flagged by NCBI) leaves residues without a genomic hit.
    if (cdsEdges.some((e) => e.attributes.exception !== undefined)) tally.noHitOnExceptionEdges++;
    else {
      tally.noHit++;
      if (disagreements.length < 20) disagreements.push(`no genomic hit for ${protein}:${r}`);
    }
  }
  for (const hit of hits) {
    tally.genomicHits++;
    const hitEdge = cdsEdges.find((e) => e.to === hit.location.outer) ?? edge;
    const flagged = hitEdge.attributes.exception !== undefined;
    const codon = extract(hit.location, genome);
    const aa = codon === undefined ? undefined : translate(codon, Number(edge.attributes.translTable ?? 1));
    const expected = seq[r - 1];
    // Selenocysteine / pyrrolysine are encoded by stop codons (transl_except).
    if (aa === expected || (aa === "*" && (expected === "U" || expected === "O")) || expected === "X") tally.agree++;
    else if (flagged) tally.disagreeOnExceptionEdges++;
    else {
      tally.disagree++;
      if (disagreements.length < 20) disagreements.push(`${protein}:${r} ${formatLocationId(hit.location, ctx)} ${codon} -> ${aa}, expected ${expected}; ${hitEdge.validation.detail ?? ""}`);
    }

    t = performance.now();
    const back = store.neighbors(hit.location, ctx);
    reverseMs.push(performance.now() - t);
    const ok = back.targets.some((x) => x.location.outer === protein && formatLocationId(x.location, ctx) === `${protein}:${r}`);
    if (ok) tally.roundTrip++;
    else if (flagged) tally.roundTripMissingOnExceptionEdges++;
    else {
      tally.roundTripMissing++;
      if (disagreements.length < 20) disagreements.push(`round trip ${protein}:${r} -> ${formatLocationId(hit.location, ctx)} -> ${back.targets.map((x) => formatLocationId(x.location, ctx)).join(" ")}`);
    }
  }
}

const pct = (xs: number[], p: number) => {
  const s = [...xs].sort((a, b) => a - b);
  return s[Math.min(s.length - 1, Math.floor((p / 100) * s.length))]!.toFixed(2);
};
console.log(JSON.stringify(tally));
console.log(`protein -> genome: p50 ${pct(forwardMs, 50)} ms, p95 ${pct(forwardMs, 95)} ms, p99 ${pct(forwardMs, 99)} ms (n=${forwardMs.length})`);
console.log(`genome -> *:       p50 ${pct(reverseMs, 50)} ms, p95 ${pct(reverseMs, 95)} ms, p99 ${pct(reverseMs, 99)} ms (n=${reverseMs.length})`);
for (const d of disagreements) console.log(`disagree: ${d}`);
