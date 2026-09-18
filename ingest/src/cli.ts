#!/usr/bin/env node
// togocoord-ingest [--db OUT.sqlite [--overwrite]] [--fasta FILE]... [--all-annotations] FILE...
// FILE may also be FASTA (.fa .faa .fasta, e.g. UniProt or Ensembl proteins): sequence records with checksums only.
// Streams GenBank/GenPept flat files (.gb .gbk .gbff .gp .gpff) and GFF3 (.gff .gff3), optionally gzipped.
// Without --db, writes JSON Lines to stdout ({"record":"sequence"|"edge"|"annotation"|"warning", ...}).
// Large or indexed FASTA (--fasta) is read on demand through a .fai index (built next to the file when missing).
import { execFileSync } from "node:child_process";
import { createReadStream, readFileSync, statSync, existsSync } from "node:fs";
import { createHash } from "node:crypto";
import { basename } from "node:path";
import { gunzipSync } from "node:zlib";
import { NamespaceRegistry } from "@togocoord/core";
import {
  accessionRef,
  assemblyBaseName,
  assemblyReportAliases,
  assemblyReportInfo,
  assemblyReportSeqids,
  assemblyReportMolecules,
  assemblyReportSequences,
  lookupSeqid,
  UCSC_DATABASES,
  ChainedSource, DEFAULT_EXCLUDED_ANNOTATIONS, VersionResolver } from "./common.ts";
import { parseFastaHeaders } from "./fasta.ts";
import { FaiSequenceSource } from "./fasta-index.ts";
import { MemorySequenceSource, type SequenceSource } from "./sequence.ts";
import { defaultFastaRef, ingestFastaFile } from "./adapter-fasta.ts";
import { ingestSiftsFile } from "./adapter-sifts.ts";
import { ingestManeSummary } from "./adapter-mane.ts";
import { ingestChainFile } from "./adapter-chain.ts";
import { ingestPafFile } from "./adapter-paf.ts";
import { ingestBedFile } from "./adapter-bed.ts";
import { SqliteSink } from "./store.ts";
import { ingestGenBankFile, ingestGff3File, JsonlSink } from "./stream.ts";

const USAGE =
  "usage: togocoord-ingest [--db OUT.sqlite [--overwrite]] [--fasta FILE]... [--seqid-map ASSEMBLY_REPORT] [--assembly-report FILE]\n" +
  "                        [--label TEXT] [--species-taxon N] [--taxon ID] [--organism NAME] [--assembly NAME] [--sifts-known-only] [--all-annotations]\n" +
  "                        [--from-report ASSEMBLY_REPORT --to-report ASSEMBLY_REPORT (for .chain / .paf files)] [--method TEXT] FILE...\n" +
  "FILE: .gbff/.gb/.gp, .gff3, .fa/.fna/.faa, SIFTS .tsv, MANE summary, UCSC .chain, PAF with cg:Z CIGAR, NCBI *_assembly_report.txt (optionally .gz)\n" +
  "BED: [--bed-type TYPE] [--bed-columns NAME,...|attributes] [--id-namespace NS (IDs as input, e.g. fanta)] [--link URL_WITH_{id}]\n" +
  "--method: how an input was made (e.g. the aligner, its version and arguments), recorded for reproduction\n";
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
/** Store metadata shown by the service (label, organism, assembly). */
const meta: Record<string, string> = {};
/** Sequence names of the --assembly-report assembly (for BED chromosome names). */
let reportSeqids: Map<string, string> | undefined;
let bedType: string | undefined;
let bedColumns: string[] | undefined;
/** Sequence files given with --fasta (for validation; their MD5 is recorded, e.g. the genomes an alignment used). */
const fastaFiles: string[] = [];
/** Sequence names of the two assemblies of a liftOver chain file. */
let fromNames: Map<string, string> | undefined;
let toNames: Map<string, string> | undefined;
/** Both assemblies' reports: the chain store records the species and length of the sequences it connects. */
const chainReports: string[] = [];
/** Sequences listed in both reports (the same accession in both assemblies): identity, no alignment needed. */
/** Nuclear sequences and organelle genomes do not correspond (organelle DNA inserted in the nucleus is paralogous). */
function sameMolecule(): (from: string, to: string) => boolean {
  const kinds = new Map(chainReports.flatMap((t) => [...assemblyReportMolecules(t)]));
  return (a, b) => !kinds.has(a) || !kinds.has(b) || kinds.get(a) === kinds.get(b);
}
function sharedBetweenReports(): (ref: string) => boolean {
  const [a, b] = chainReports.map((t) => new Set(assemblyReportSequences(t).map((r) => r.ref)));
  return (ref) => !!a && !!b && a.has(ref) && b.has(ref);
}
for (let i = 0; i < args.length; i++) {
  const a = args[i]!;
  if (a === "--db") db = args[++i];
  else if (a === "--overwrite") overwrite = true;
  else if (a === "--all-annotations") allAnnotations = true;
  else if (a === "--sifts-known-only") siftsKnownOnly = true;
  else if (a === "--seqid-map" || a === "--assembly-report") {
    const text = readFileSync(args[++i]!, "utf8");
    if (a === "--seqid-map") seqids = assemblyReportSeqids(text);
    reportSeqids ??= assemblyReportSeqids(text);
    const info = assemblyReportInfo(text);
    for (const [k, v] of Object.entries(info)) meta[k] ??= v;
    // Sequence names of the assembly (chr7, 7, CM000669.2), so input may be written like hg38:chr7:140753336.
    meta.aliases ??= JSON.stringify(assemblyReportAliases(text));
    const ucsc = info.assembly && UCSC_DATABASES[assemblyBaseName(info.assembly)];
    if (ucsc) meta.ucsc ??= ucsc;
  } else if (a === "--from-report" || a === "--to-report") {
    const text = readFileSync(args[++i]!, "utf8");
    const names = assemblyReportSeqids(text);
    chainReports.push(text);
    if (a === "--from-report") fromNames = names;
    else toNames = names;
  } else if (a === "--bed-type") bedType = args[++i];
  else if (a === "--bed-columns") bedColumns = args[++i]!.split(",");
  else if (["--label", "--taxon", "--organism", "--assembly", "--species-taxon", "--method", "--id-namespace", "--link"].includes(a)) meta[a.slice(2).replace("-", "_")] = args[++i]!;
  else if (a === "--fasta") {
    fastaFiles.push(args[i + 1]!);
    sources.push(openFasta(args[++i]));
  }
  else if (a.startsWith("--")) {
    process.stderr.write(`unknown option ${a}\n${USAGE}`);
    process.exit(1);
  } else inputs.push(a);
}

function openFasta(file: string | undefined): SequenceSource {
  if (!file) throw new Error("--fasta needs a file");
  // Large gzipped FASTA cannot be held as one string (V8 limit ~0.5 G characters): decompress it once next to the
  // input and read it through a .fai index like any large FASTA.
  if (file.endsWith(".gz") && statSync(file).size > 32 * 1024 * 1024) {
    const plain = file.slice(0, -3);
    if (!existsSync(plain)) {
      process.stderr.write(`${file}: decompressing to ${plain} for indexed access\n`);
      execFileSync("sh", ["-c", 'gzip -dc "$1" > "$2.tmp" && mv "$2.tmp" "$2"', "sh", file, plain]);
    }
    file = plain;
  }
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
  if (/\.chain$/i.test(name)) {
    if (!fromNames || !toNames) throw new Error(`${file}: chain files need --from-report and --to-report (NCBI assembly reports of both assemblies)`);
    for (const text of chainReports) for (const r of assemblyReportSequences(text, basename(file))) sink.sequence(r);
    const s = await ingestChainFile(file, sink, {
      file: basename(file),
      registry,
      source,
      fromRef: (n) => lookupSeqid(fromNames!, n),
      toRef: (n) => lookupSeqid(toNames!, n),
      shared: sharedBetweenReports(),
      compatible: sameMolecule(),
    });
    const identity = s.sampledBases ? ((100 * s.identicalBases) / s.sampledBases).toFixed(1) : "-";
    process.stderr.write(
      `${file}: ${s.chains} chains, ${s.blocks} blocks, ${s.skipped} skipped (${s.shared ?? 0} from sequences shared by both assemblies, ${s.crossMolecule ?? 0} between nuclear and organelle genomes); sampled identity ${identity}%\n`,
    );
    continue;
  }
  if (/assembly_report\.txt$/i.test(name)) {
    // The sequences of an assembly (length, species) and its names, e.g. an assembly without annotation (GRCh37).
    const text = file.endsWith(".gz") ? gunzipSync(readFileSync(file)).toString("utf8") : readFileSync(file, "utf8");
    const records = assemblyReportSequences(text, basename(file));
    for (const r of records) sink.sequence(r);
    const info = assemblyReportInfo(text);
    for (const [k, v] of Object.entries(info)) meta[k] ??= v;
    meta.aliases ??= JSON.stringify(assemblyReportAliases(text));
    const ucsc = info.assembly && UCSC_DATABASES[assemblyBaseName(info.assembly)];
    if (ucsc) meta.ucsc ??= ucsc;
    process.stderr.write(`${file}: ${records.length} sequences of ${info.assembly ?? "the assembly"}\n`);
    continue;
  }
  if (/\.paf$/i.test(name)) {
    if (!fromNames || !toNames) throw new Error(`${file}: PAF files need --from-report (query assembly) and --to-report (target assembly)`);
    for (const text of chainReports) for (const r of assemblyReportSequences(text, basename(file))) sink.sequence(r);
    const s = await ingestPafFile(file, sink, {
      file: basename(file),
      registry,
      source,
      fromRef: (n) => lookupSeqid(fromNames!, n),
      toRef: (n) => lookupSeqid(toNames!, n),
      shared: sharedBetweenReports(),
      compatible: sameMolecule(),
    });
    const identity = s.sampledBases ? ((100 * s.identicalBases) / s.sampledBases).toFixed(2) : "-";
    process.stderr.write(
      `${file}: ${s.records} records, ${s.alignments} alignments kept, ${s.blocks} blocks, ${s.skipped} skipped ` +
        `(${s.shared ?? 0} shared by both assemblies, ${s.crossMolecule ?? 0} between nuclear and organelle genomes), ` +
        `${s.overlapBases} source bases covered by better alignments; sampled identity ${identity}%\n`,
    );
    continue;
  }
  if (/\.bed$/i.test(name)) {
    const names = seqids ?? reportSeqids;
    if (!names) throw new Error(`${file}: BED files need --assembly-report (chromosome names of the assembly)`);
    const s = await ingestBedFile(file, sink, {
      file: basename(file),
      registry,
      refOf: (c) => lookupSeqid(names, c),
      ...(bedType && { type: bedType }),
      ...(bedColumns && { extraColumns: bedColumns }),
    });
    process.stderr.write(`${file}: ${s.records} regions, ${s.skipped} skipped\n`);
    continue;
  }
  if (/MANE.*summary\.txt$/i.test(name)) {
    const s = await ingestManeSummary(file, sink, { file: basename(file), registry });
    process.stderr.write(`${file}: ${s.genes} MANE genes, ${s.sequences} tagged sequences\n`);
    continue;
  }
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
/** Provenance for reproduction: MD5 of every input and the TogoCoord commit that built the store. */
async function md5(file: string): Promise<string> {
  const hash = createHash("md5");
  for await (const chunk of createReadStream(file)) hash.update(chunk as Buffer);
  return hash.digest("hex");
}
function commit(): string | undefined {
  try {
    const dir = new URL("../..", import.meta.url).pathname; // the repository
    const id = execFileSync("git", ["-C", dir, "rev-parse", "--short", "HEAD"], { encoding: "utf8" }).trim();
    const dirty = execFileSync("git", ["-C", dir, "status", "--porcelain", "--", "."], { encoding: "utf8" }).trim() !== "";
    return `${id}${dirty ? "+local changes" : ""}`;
  } catch {
    return undefined;
  }
}
if (sink instanceof SqliteSink) {
  const checksums = Object.fromEntries(await Promise.all(inputs.map(async (f) => [basename(f), await md5(f)])));
  const built = commit();
  const sequences = Object.fromEntries(await Promise.all(fastaFiles.map(async (f) => [basename(f), await md5(f)])));
  sink.close({
    ...meta,
    inputs: inputs.map((f) => basename(f)).join(","),
    inputs_md5: JSON.stringify(checksums),
    ...(fastaFiles.length && { sequences_md5: JSON.stringify(sequences) }),
    ...(built && { togocoord: built }),
  });
}
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
