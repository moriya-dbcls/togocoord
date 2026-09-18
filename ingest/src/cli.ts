#!/usr/bin/env node
// togocoord-ingest [--db OUT.sqlite [--overwrite]] [--fasta FILE]... [--all-annotations] FILE...
// FILE may also be FASTA (.fa .faa .fasta, e.g. UniProt or Ensembl proteins): sequence records with checksums only.
// Streams GenBank/GenPept flat files (.gb .gbk .gbff .gp .gpff) and GFF3 (.gff .gff3), optionally gzipped.
// Without --db, writes JSON Lines to stdout ({"record":"sequence"|"edge"|"annotation"|"warning", ...}).
// Large or indexed FASTA (--fasta) is read on demand through a .fai index (built next to the file when missing).
import { readFileSync, statSync, existsSync } from "node:fs";
import { basename } from "node:path";
import { gunzipSync } from "node:zlib";
import { NamespaceRegistry } from "@togocoord/core";
import { accessionRef, assemblyReportSeqids, ChainedSource, DEFAULT_EXCLUDED_ANNOTATIONS, VersionResolver } from "./common.ts";
import { parseFastaHeaders } from "./fasta.ts";
import { FaiSequenceSource } from "./fasta-index.ts";
import { MemorySequenceSource, type SequenceSource } from "./sequence.ts";
import { defaultFastaRef, ingestFastaFile } from "./adapter-fasta.ts";
import { ingestSiftsFile } from "./adapter-sifts.ts";
import { SqliteSink } from "./store.ts";
import { ingestGenBankFile, ingestGff3File, JsonlSink } from "./stream.ts";

const USAGE =
  "usage: togocoord-ingest [--db OUT.sqlite [--overwrite]] [--fasta FILE]... [--seqid-map ASSEMBLY_REPORT] [--sifts-known-only] [--all-annotations] FILE...\n";
const args = process.argv.slice(2);
if (args.length === 0 || args.includes("-h") || args.includes("--help")) {
  process.stderr.write(USAGE);
  process.exit(args.length === 0 ? 1 : 0);
}

const registry = new NamespaceRegistry();
/** FASTA name (first header word) -> key; pdb_seqres names like `101m_A` become pdb chains. */
const nameToRef = (name: string) => {
  const pdb = /^([0-9][A-Za-z0-9]{3})_(\S+)$/.exec(name);
  if (pdb) {
    try {
      return registry.refKey("pdb", `${pdb[1]}.${pdb[2]}`);
    } catch {
      return undefined;
    }
  }
  return defaultFastaRef(name, registry);
};
const sources: SequenceSource[] = [];
const inputs: string[] = [];
let db: string | undefined;
let overwrite = false;
let allAnnotations = false;
let siftsKnownOnly = false;
let seqids: Map<string, string> | undefined;
for (let i = 0; i < args.length; i++) {
  const a = args[i]!;
  if (a === "--db") db = args[++i];
  else if (a === "--overwrite") overwrite = true;
  else if (a === "--all-annotations") allAnnotations = true;
  else if (a === "--sifts-known-only") siftsKnownOnly = true;
  else if (a === "--seqid-map") seqids = assemblyReportSeqids(readFileSync(args[++i]!, "utf8"));
  else if (a === "--fasta") sources.push(openFasta(args[++i]));
  else if (a.startsWith("--")) {
    process.stderr.write(`unknown option ${a}\n${USAGE}`);
    process.exit(1);
  } else inputs.push(a);
}

function openFasta(file: string | undefined): SequenceSource {
  if (!file) throw new Error("--fasta needs a file");
  if (!file.endsWith(".gz") && (existsSync(`${file}.fai`) || statSync(file).size > 64 * 1024 * 1024)) {
    return new FaiSequenceSource(file, nameToRef);
  }
  const text = file.endsWith(".gz") ? gunzipSync(readFileSync(file)).toString("utf8") : readFileSync(file, "utf8");
  const source = new MemorySequenceSource();
  let ignored = 0;
  for (const [id, residues] of parseFastaHeaders(text)) {
    const ref = defaultFastaRef(id, registry);
    if (ref) source.add(ref, residues);
    else ignored++;
  }
  if (ignored) process.stderr.write(`warning: ${file}: ${ignored} FASTA records without a recognised accession ignored\n`);
  return source;
}

// Version-less keys in GFF3 attributes (Ensembl protein_id) are completed from the FASTA sources.
const versions = new VersionResolver();
for (const s of sources) for (const r of s.refs?.() ?? []) versions.add(r);
const source = new ChainedSource(...sources);

const sink = db ? new SqliteSink(db, { overwrite }) : new JsonlSink(1);
const started = performance.now();
for (const file of inputs) {
  const options = {
    file: basename(file),
    registry,
    source,
    excludeAnnotations: allAnnotations ? new Set<string>() : DEFAULT_EXCLUDED_ANNOTATIONS,
    resolveRef: (r: string) => versions.resolve(r),
    ...(seqids && { seqidToRef: (seqid: string) => seqids!.get(seqid) ?? accessionRef(seqid, registry) }),
  };
  const name = file.replace(/\.gz$/, "");
  if (/(?:sifts|uniprot_segments).*\.tsv$/i.test(name)) {
    const accept = siftsKnownOnly ? (acc: string) => source.length(`uniprot:${acc}`) !== undefined : undefined;
    const s = await ingestSiftsFile(file, sink, { file: basename(file), registry, source, ...(accept && { accept }) });
    process.stderr.write(`${file}: ${s.rows} SIFTS rows, ${s.edges} edges, ${s.skipped} skipped\n`);
    continue;
  }
  if (/\.(fa|fasta|faa|fna|pep)$/i.test(name)) {
    const s = await ingestFastaFile(file, sink, { file: basename(file), registry });
    process.stderr.write(`${file}: ${s.records} sequences (checksums only), ${s.skipped} skipped\n`);
    continue;
  }
  const stats = /\.(gff3?|gff3\.txt)$/i.test(name)
    ? await ingestGff3File(file, sink, options)
    : /\.(gb|gbk|gbff|gp|gpff|genbank)$/i.test(name)
      ? await ingestGenBankFile(file, sink, options)
      : undefined;
  if (!stats) throw new Error(`${file}: unknown format (expected .gb/.gbff/.gp/.gff3/.fa, optionally .gz)`);
  process.stderr.write(`${file}: ${stats.lines} lines, ${stats.features} features\n`);
}
if (sink instanceof SqliteSink) sink.close({ inputs: inputs.map((f) => basename(f)).join(",") });
else sink.flush();

const seconds = ((performance.now() - started) / 1000).toFixed(1);
process.stderr.write(
  `sequences ${sink.counts.sequence}, edges ${sink.counts.edge} ${JSON.stringify(sink.validation)}, ` +
    `annotations ${sink.counts.annotation}, warnings ${sink.counts.warning} (${seconds} s, max RSS ${(process.resourceUsage().maxRSS / 1024).toFixed(0)} MB)\n`,
);
const expected = sink.mismatches.filter((m) => m.includes("expected: /exception="));
const unexpected = sink.mismatches.filter((m) => !m.includes("expected: /exception="));
if (sink.mismatches.length) {
  process.stderr.write(`mismatches: ${unexpected.length} unexplained, ${expected.length} with an INSDC /exception\n`);
}
for (const m of unexpected.slice(0, 20)) process.stderr.write(`mismatch: ${m}\n`);
if (unexpected.length > 20) process.stderr.write(`... ${unexpected.length - 20} more unexplained mismatches\n`);
