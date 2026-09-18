// Data for the TogoCoord stores: download the primary files and build the stores (docs/data.md).
//   node scripts/data.ts list                       datasets, their groups and whether they are built
//   node scripts/data.ts download [NAME...]         the files the named stores (or groups) need; all without names
//   node scripts/data.ts build [NAME...] [--force]  build the stores (missing ones; --force rebuilds)
//   node scripts/data.ts serve [serve options]      start togocoord-serve with the built stores, in scripts/stores.txt order
// NAME: a store (human_rna) or a group (human). Files live under $TOGOCOORD_DATA (default: <repo>/data), which is not
// in git: raw/ downloads (upstream file names), work/ intermediate files (alignments), stores/ SQLite stores, logs/.
import { spawn } from "node:child_process";
import { createWriteStream, existsSync, mkdirSync, readFileSync, renameSync, statSync } from "node:fs";
import { basename, join } from "node:path";
import { Readable } from "node:stream";
import { pipeline } from "node:stream/promises";

const REPO = new URL("..", import.meta.url).pathname;
const DATA = process.env.TOGOCOORD_DATA ?? join(REPO, "data");
const RAW = join(DATA, "raw");
const WORK = join(DATA, "work");
const STORES = join(DATA, "stores");
const LOGS = join(DATA, "logs");

// ---- sources ------------------------------------------------------------------------------------------------------
const NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all";
/** NCBI Datasets FTP prefix of an assembly: GCF_000001405.40 + GRCh38.p14 -> .../GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14 */
const assembly = (acc: string, name: string) => {
  const d = acc.slice(4, 13);
  return `${NCBI}/${acc.slice(0, 3)}/${d.slice(0, 3)}/${d.slice(3, 6)}/${d.slice(6, 9)}/${acc}_${name}/${acc}_${name}`;
};
const GRCH38 = assembly("GCF_000001405.40", "GRCh38.p14");
const GRCH37 = assembly("GCF_000001405.25", "GRCh37.p13");
const GRCM39 = assembly("GCF_000001635.27", "GRCm39");
const GRCM38 = assembly("GCF_000001635.26", "GRCm38.p6");
const TAIR = assembly("GCF_000001735.4", "TAIR10.1");
const MP31 = assembly("GCA_003032435.1", "Marchanta_polymorpha_v1");
const MP71 = assembly("GCA_039105155.1", "MpTak_v7.1");
const UNIPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota";
const ENSEMBL = "https://ftp.ensembl.org/pub/release-116";
const MANE = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5";
const UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath";
const FANTA = "https://data.fanta.bio/cre/v1.2.1";

/** Downloaded files by local name (the upstream file name). */
const RAW_FILES = new Map<string, string>();
const raw = (url: string): string => {
  const name = basename(url);
  RAW_FILES.set(name, url);
  return join(RAW, name);
};
const ncbi = (prefix: string, suffix: string) => raw(`${prefix}_${suffix}`);
const uniprot = (proteome: string, taxon: number, additional = false) =>
  raw(`${UNIPROT}/${proteome}/${proteome}_${taxon}${additional ? "_additional" : ""}.fasta.gz`);
const work = (name: string) => join(WORK, name);

// ---- stores -------------------------------------------------------------------------------------------------------
interface Store {
  name: string;
  group: string;
  /** Arguments of togocoord-ingest (without --db / --overwrite). */
  args: string[];
  /** Files made before ingesting (e.g. an alignment), with the command that makes them. */
  prepare?: Array<{ file: string; command: string[]; inputs: string[] }>;
}

const minimap2 = "minimap2 2.31-r1302 -c --eqx -x asm5 -t 8";
const pafFilter = "PAF filtered by TogoCoord: secondary (tp:A:S) and alignments < 1 kb dropped; one-to-one on the query, best AS first";
const align = (target: string, query: string, out: string) => ({
  file: out,
  command: ["sh", "-c", `minimap2 -c --eqx -x asm5 -t 8 "$0" "$1" > "$2.tmp" && mv "$2.tmp" "$2"`, target, query, out],
  inputs: [target, query],
});

const r38 = ncbi(GRCH38, "assembly_report.txt");
const r37 = ncbi(GRCH37, "assembly_report.txt");
const rm39 = ncbi(GRCM39, "assembly_report.txt");
const rm38 = ncbi(GRCM38, "assembly_report.txt");
const rTair = ncbi(TAIR, "assembly_report.txt");
const rMp31 = ncbi(MP31, "assembly_report.txt");
const rMp71 = ncbi(MP71, "assembly_report.txt");
const g38 = ncbi(GRCH38, "genomic.fna.gz");
const g37 = ncbi(GRCH37, "genomic.fna.gz");
const gm39 = ncbi(GRCM39, "genomic.fna.gz");
const gm38 = ncbi(GRCM38, "genomic.fna.gz");
const gMp31 = ncbi(MP31, "genomic.fna.gz");
const gMp71 = ncbi(MP71, "genomic.fna.gz");
const upHuman = [uniprot("UP000005640", 9606), uniprot("UP000005640", 9606, true)];
const upMouse = [uniprot("UP000000589", 10090), uniprot("UP000000589", 10090, true)];
const upArab = [uniprot("UP000006548", 3702), uniprot("UP000006548", 3702, true)];
const upMarchantia = [uniprot("UP000244005", 3197)];
const paf31to71 = work("MpTak_v3.1_to_v7.1.paf");
const paf71to31 = work("MpTak_v7.1_to_v3.1.paf");
const cre = ["--bed-type", "CRE", "--bed-columns", "Name,attributes", "--id-namespace", "fanta", "--link", "https://fanta.bio/cre/{id}"];

const STORE_LIST: Store[] = [
  // Human, GRCh38 (annotated)
  { name: "human", group: "human", args: ["--label", "Human RefSeq annotation (GRCh38.p14)", "--assembly-report", r38, "--fasta", g38, "--fasta", ncbi(GRCH38, "rna.fna.gz"), "--fasta", ncbi(GRCH38, "protein.faa.gz"), ncbi(GRCH38, "genomic.gff.gz")] },
  { name: "human_rna", group: "human", args: ["--label", "Human RefSeq RNA", ncbi(GRCH38, "rna.gbff.gz")] },
  { name: "human_ensembl", group: "human", args: ["--label", "Human Ensembl release 116 (GRCh38)", "--seqid-map", r38, "--fasta", g38, "--fasta", raw(`${ENSEMBL}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz`), raw(`${ENSEMBL}/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz`)] },
  { name: "human_uniprot", group: "human", args: ["--label", "Human UniProt reference proteome (UP000005640)", ...upHuman] },
  { name: "human_mane", group: "human", args: ["--label", "MANE v1.5 (human)", "--taxon", "9606", "--organism", "Homo sapiens", raw(`${MANE}/MANE.GRCh38.v1.5.summary.txt.gz`), raw(`${MANE}/MANE.GRCh38.v1.5.refseq_rna.fna.gz`), raw(`${MANE}/MANE.GRCh38.v1.5.ensembl_rna.fna.gz`)] },
  { name: "grch38_names", group: "human", args: ["--label", "GRCh38.p14 sequence names (NCBI assembly report)", r38] },
  // Human, GRCh37 and UCSC chains to / from GRCh38
  { name: "grch37", group: "grch37", args: ["--label", "GRCh37.p13 genome (NCBI assembly report and sequences)", "--assembly-report", r37, r37, g37] },
  { name: "chain_hg19ToHg38", group: "grch37", args: ["--label", "UCSC liftOver chains hg19 → hg38 (GRCh37 to GRCh38)", "--from-report", r37, "--to-report", r38, "--fasta", g37, "--fasta", g38, raw(`${UCSC}/hg19/liftOver/hg19ToHg38.over.chain.gz`)] },
  { name: "chain_hg38ToHg19", group: "grch37", args: ["--label", "UCSC liftOver chains hg38 → hg19 (GRCh38 to GRCh37)", "--from-report", r38, "--to-report", r37, "--fasta", g38, "--fasta", g37, raw(`${UCSC}/hg38/liftOver/hg38ToHg19.over.chain.gz`)] },
  // Mouse, GRCm39 (annotated)
  { name: "mouse", group: "mouse", args: ["--label", "Mouse RefSeq annotation (GRCm39)", "--assembly-report", rm39, "--fasta", gm39, "--fasta", ncbi(GRCM39, "rna.fna.gz"), "--fasta", ncbi(GRCM39, "protein.faa.gz"), ncbi(GRCM39, "genomic.gff.gz")] },
  { name: "mouse_rna", group: "mouse", args: ["--label", "Mouse RefSeq RNA", ncbi(GRCM39, "rna.gbff.gz")] },
  { name: "mouse_uniprot", group: "mouse", args: ["--label", "Mouse UniProt reference proteome (UP000000589)", ...upMouse] },
  { name: "grcm39_names", group: "mouse", args: ["--label", "GRCm39 sequence names (NCBI assembly report)", rm39] },
  // Mouse, GRCm38 and UCSC chains to / from GRCm39
  { name: "grcm38", group: "grcm38", args: ["--label", "GRCm38.p6 genome (NCBI assembly report and sequences)", "--assembly-report", rm38, rm38, gm38] },
  { name: "chain_mm10ToMm39", group: "grcm38", args: ["--label", "UCSC liftOver chains mm10 → mm39 (GRCm38 to GRCm39)", "--from-report", rm38, "--to-report", rm39, "--fasta", gm38, "--fasta", gm39, raw(`${UCSC}/mm10/liftOver/mm10ToMm39.over.chain.gz`)] },
  { name: "chain_mm39ToMm10", group: "grcm38", args: ["--label", "UCSC liftOver chains mm39 → mm10 (GRCm39 to GRCm38)", "--from-report", rm39, "--to-report", rm38, "--fasta", gm39, "--fasta", gm38, raw(`${UCSC}/mm39/liftOver/mm39ToMm10.over.chain.gz`)] },
  // Human <-> mouse
  { name: "chain_hg38ToMm39", group: "human_mouse", args: ["--label", "UCSC liftOver chains hg38 → mm39 (human to mouse)", "--from-report", r38, "--to-report", rm39, "--fasta", g38, "--fasta", gm39, raw(`${UCSC}/hg38/liftOver/hg38ToMm39.over.chain.gz`)] },
  { name: "chain_mm39ToHg38", group: "human_mouse", args: ["--label", "UCSC liftOver chains mm39 → hg38 (mouse to human)", "--from-report", rm39, "--to-report", r38, "--fasta", gm39, "--fasta", g38, raw(`${UCSC}/mm39/liftOver/mm39ToHg38.over.chain.gz`)] },
  // Arabidopsis
  { name: "arabidopsis", group: "arabidopsis", args: ["--label", "Arabidopsis RefSeq annotation (TAIR10.1)", "--assembly-report", rTair, "--fasta", ncbi(TAIR, "genomic.fna.gz"), "--fasta", ncbi(TAIR, "rna.fna.gz"), "--fasta", ncbi(TAIR, "protein.faa.gz"), ncbi(TAIR, "genomic.gff.gz")] },
  { name: "arabidopsis_rna", group: "arabidopsis", args: ["--label", "Arabidopsis RefSeq RNA", ncbi(TAIR, "rna.gbff.gz")] },
  { name: "arabidopsis_uniprot", group: "arabidopsis", args: ["--label", "Arabidopsis UniProt reference proteome (UP000006548)", ...upArab] },
  // Marchantia: v7.1 (default) and v3.1 (INSDC only), joined by minimap2 alignments
  { name: "marchantia_v71", group: "marchantia", args: ["--label", "Marchantia polymorpha MpTak_v7.1 INSDC annotation (GCA_039105155.1)", "--assembly-report", rMp71, ncbi(MP71, "genomic.gbff.gz")] },
  { name: "marchantia", group: "marchantia", args: ["--label", "Marchantia polymorpha INSDC annotation (GCA_003032435.1)", "--assembly-report", rMp31, ncbi(MP31, "genomic.gbff.gz")] },
  { name: "marchantia_uniprot", group: "marchantia", args: ["--label", "Marchantia UniProt proteome (UP000244005)", ...upMarchantia] },
  {
    name: "mp_v31_to_v71",
    group: "marchantia",
    prepare: [align(gMp71, gMp31, paf31to71)],
    args: ["--label", "minimap2 alignment: Marchantia MpTak v3.1 → v7.1", "--method", `${minimap2} ${basename(gMp71)} ${basename(gMp31)} (target v7.1, query v3.1). ${pafFilter}`, "--from-report", rMp31, "--to-report", rMp71, "--fasta", gMp31, "--fasta", gMp71, paf31to71],
  },
  {
    name: "mp_v71_to_v31",
    group: "marchantia",
    prepare: [align(gMp31, gMp71, paf71to31)],
    args: ["--label", "minimap2 alignment: Marchantia MpTak v7.1 → v3.1", "--method", `${minimap2} ${basename(gMp31)} ${basename(gMp71)} (target v3.1, query v7.1). ${pafFilter}`, "--from-report", rMp71, "--to-report", rMp31, "--fasta", gMp71, "--fasta", gMp31, paf71to31],
  },
  // Structures of the loaded proteomes
  {
    name: "sifts",
    group: "sifts",
    args: [
      "--label", "SIFTS: UniProt to PDB chains (loaded proteomes)", "--sifts-known-only",
      ...[...upHuman, ...upMouse, ...upArab, ...upMarchantia].flatMap((f) => ["--fasta", f]),
      "--fasta", raw("https://files.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"),
      raw("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/uniprot_segments_observed.tsv.gz"),
    ],
  },
  // Cis-regulatory elements (fanta.bio)
  { name: "fanta_human_hg38", group: "fanta", args: ["--label", "fanta.bio CREs v1.2.1 (human, hg38)", "--assembly-report", r38, ...cre, raw(`${FANTA}/human/human-CREv1.2.1.hg38.cre-peaks.bed.gz`)] },
  { name: "fanta_mouse_mm10", group: "fanta", args: ["--label", "fanta.bio CREs v1.2.1 (mouse, mm10)", "--assembly-report", rm38, ...cre, raw(`${FANTA}/mouse/mouse-CREv1.2.1.mm10.cre-peaks.bed.gz`)] },
];

// ---- commands -----------------------------------------------------------------------------------------------------
const select = (names: string[]): Store[] => {
  if (names.length === 0) return STORE_LIST;
  const out = STORE_LIST.filter((s) => names.includes(s.name) || names.includes(s.group));
  const unknown = names.filter((n) => !STORE_LIST.some((s) => s.name === n || s.group === n));
  if (unknown.length) throw new Error(`unknown store or group: ${unknown.join(", ")} (see: node scripts/data.ts list)`);
  return out;
};
const rawOf = (s: Store) => [...s.args, ...(s.prepare ?? []).flatMap((p) => p.inputs)].filter((a) => a.startsWith(RAW + "/"));
const storePath = (name: string) => join(STORES, `${name}.sqlite`);
const size = (f: string) => (existsSync(f) ? `${(statSync(f).size / 1e6).toFixed(1)} MB` : "-");

async function download(stores: Store[]): Promise<void> {
  mkdirSync(RAW, { recursive: true });
  for (const file of [...new Set(stores.flatMap(rawOf))]) {
    if (existsSync(file)) continue;
    const url = RAW_FILES.get(basename(file))!;
    process.stderr.write(`download ${url}\n`);
    const res = await fetch(url);
    if (!res.ok || !res.body) throw new Error(`${url}: HTTP ${res.status}`);
    await pipeline(Readable.fromWeb(res.body as never), createWriteStream(`${file}.part`));
    renameSync(`${file}.part`, file);
  }
}

function run(command: string[], log: string): Promise<void> {
  return new Promise((resolve, reject) => {
    const out = createWriteStream(log, { flags: "a" });
    const child = spawn(command[0]!, command.slice(1), { stdio: ["ignore", "inherit", "pipe"] });
    child.stderr.on("data", (d: Buffer) => {
      out.write(d);
      process.stderr.write(d);
    });
    child.on("error", reject);
    child.on("close", (code) => {
      out.end();
      code === 0 ? resolve() : reject(new Error(`${command.slice(0, 3).join(" ")} … exited with ${code} (log: ${log})`));
    });
  });
}

async function build(stores: Store[], force: boolean): Promise<void> {
  for (const d of [WORK, STORES, LOGS]) mkdirSync(d, { recursive: true });
  await download(stores);
  for (const s of stores) {
    if (existsSync(storePath(s.name)) && !force) {
      process.stderr.write(`${s.name}: built (use --force to rebuild)\n`);
      continue;
    }
    const log = join(LOGS, `${s.name}.log`);
    for (const p of s.prepare ?? []) {
      if (existsSync(p.file) && !force) continue;
      process.stderr.write(`${s.name}: ${p.command.slice(2, 3).join(" ")}\n`);
      await run(p.command, log);
    }
    process.stderr.write(`${s.name}: building\n`);
    const cli = join(REPO, "ingest/src/cli.ts");
    await run(["node", "--no-warnings", "--max-old-space-size=6144", cli, "--db", storePath(s.name), "--overwrite", ...s.args], log);
  }
}

/** Store names in serving order (scripts/stores.txt). */
function servingOrder(): string[] {
  return readFileSync(join(REPO, "scripts/stores.txt"), "utf8")
    .split("\n")
    .map((l) => l.replace(/#.*/, "").trim())
    .filter(Boolean);
}

const [command, ...rest] = process.argv.slice(2);
const force = rest.includes("--force");
const names = rest.filter((a) => !a.startsWith("--"));
switch (command) {
  case "list":
    for (const s of STORE_LIST) console.log(`${s.group.padEnd(12)} ${s.name.padEnd(22)} ${size(storePath(s.name)).padStart(10)}  ${s.args[1] ?? ""}`);
    console.log(`\ndata: ${DATA}`);
    break;
  case "download":
    await download(select(names));
    break;
  case "build":
    await build(select(names), force);
    break;
  case "serve": {
    const order = servingOrder();
    const unlisted = STORE_LIST.filter((s) => !order.includes(s.name)).map((s) => s.name);
    if (unlisted.length) process.stderr.write(`warning: not in scripts/stores.txt: ${unlisted.join(", ")}\n`);
    await run(["node", "--no-warnings", join(REPO, "service/src/serve.ts"), "--stores", join(REPO, "scripts/stores.txt"), "--store-dir", STORES, "--skip-missing", ...rest], join(LOGS, "serve.log"));
    break;
  }
  default:
    process.stderr.write(readFileSync(new URL(import.meta.url), "utf8").split("\n").slice(0, 8).join("\n") + "\n");
    process.exit(command ? 1 : 0);
}
