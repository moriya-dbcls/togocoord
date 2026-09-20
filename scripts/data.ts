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
const CHM13 = assembly("GCF_009914755.1", "T2T-CHM13v2.0");
const GRCM39 = assembly("GCF_000001635.27", "GRCm39");
const GRCM38 = assembly("GCF_000001635.26", "GRCm38.p6");
const TAIR = assembly("GCF_000001735.4", "TAIR10.1");
const TAIR10 = assembly("GCF_000001735.3", "TAIR10");
const MP31 = assembly("GCA_003032435.1", "Marchanta_polymorpha_v1");
const MP71 = assembly("GCA_039105155.1", "MpTak_v7.1");
const MPTAK2 = assembly("GCA_037833965.1", "MpTak2_v7.1");
const MP51 = assembly("GCA_009936355.2", "ASM993635v2");
const MPCM = assembly("GCA_965642975.2", "cmMarPoly1.2");
const MPV4 = assembly("GCA_001641455.1", "Mp_v4");
const UNIPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota";
const ENSEMBL = "https://ftp.ensembl.org/pub/release-116";
const MANE = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5";
const UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath";
const FANTA = "https://data.fanta.bio/cre/v1.2.1";
const IDMAPPING = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism";

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

const minimap2 = (preset: string) => `minimap2 2.31-r1302 -c --eqx -x ${preset} -t 8`;
const pafFilter = "PAF filtered by TogoCoord: secondary (tp:A:S) and alignments < 1 kb dropped; one-to-one on the query, best AS first";
/** minimap2 presets by sequence divergence: asm5 ~0.1% (versions of one strain), asm10 ~1%, asm20 ~5%. */
const align = (target: string, query: string, out: string, preset = "asm5") => ({
  file: out,
  command: ["sh", "-c", `minimap2 -c --eqx -x "$3" -t 8 "$0" "$1" > "$2.tmp" && mv "$2.tmp" "$2"`, target, query, out, preset],
  inputs: [target, query],
});

const r38 = ncbi(GRCH38, "assembly_report.txt");
const r37 = ncbi(GRCH37, "assembly_report.txt");
const rm39 = ncbi(GRCM39, "assembly_report.txt");
const rm38 = ncbi(GRCM38, "assembly_report.txt");
const rTair = ncbi(TAIR, "assembly_report.txt");
const rTair10 = ncbi(TAIR10, "assembly_report.txt");
const gTair = ncbi(TAIR, "genomic.fna.gz");
const gTair10 = ncbi(TAIR10, "genomic.fna.gz");
const pafTair10 = work("TAIR10_to_TAIR10.1.paf");
const rMp31 = ncbi(MP31, "assembly_report.txt");
const rMp71 = ncbi(MP71, "assembly_report.txt");
const g38 = ncbi(GRCH38, "genomic.fna.gz");
const g37 = ncbi(GRCH37, "genomic.fna.gz");
const rChm13 = ncbi(CHM13, "assembly_report.txt");
const gChm13 = ncbi(CHM13, "genomic.fna.gz");
const gm39 = ncbi(GRCM39, "genomic.fna.gz");
const gm38 = ncbi(GRCM38, "genomic.fna.gz");
const gMp31 = ncbi(MP31, "genomic.fna.gz");
const gMp71 = ncbi(MP71, "genomic.fna.gz");
const upHuman = [uniprot("UP000005640", 9606), uniprot("UP000005640", 9606, true)];
const upMouse = [uniprot("UP000000589", 10090), uniprot("UP000000589", 10090, true)];
const upArab = [uniprot("UP000006548", 3702), uniprot("UP000006548", 3702, true)];
const upMarchantia = [uniprot("UP000244005", 3197)];
const upMarchantiaV4 = [uniprot("UP000077202", 1480154)];
const rMpTak2 = ncbi(MPTAK2, "assembly_report.txt");
const rMp51 = ncbi(MP51, "assembly_report.txt");
const rMpCm = ncbi(MPCM, "assembly_report.txt");
const rMpV4 = ncbi(MPV4, "assembly_report.txt");
const gMpTak2 = ncbi(MPTAK2, "genomic.fna.gz");
const gMp51 = ncbi(MP51, "genomic.fna.gz");
const gMpCm = ncbi(MPCM, "genomic.fna.gz");
const gMpV4 = ncbi(MPV4, "genomic.fna.gz");
const paf31to71 = work("MpTak_v3.1_to_v7.1.paf");
const paf71to31 = work("MpTak_v7.1_to_v3.1.paf");
const cre = ["--bed-type", "CRE", "--bed-columns", "Name,attributes", "--id-namespace", "fanta", "--link", "https://fanta.bio/cre/{id}"];

/**
 * Alignments of other assemblies of a species to its annotated one, in both directions (a star: another assembly is
 * reached through the annotated one, spec-ingest §20).
 */
function starAlignments(
  group: string,
  hubName: string,
  hubReport: string,
  hubGenome: string,
  spokes: Array<{ name: string; label: string; report: string; genome: string; preset?: string }>,
): Store[] {
  return spokes.flatMap((s) => [
    {
      name: `mp_${s.name}_to_hub`,
      group,
      prepare: [align(hubGenome, s.genome, work(`${s.name}_to_${hubName}.paf`), s.preset)],
      args: ["--label", `minimap2 alignment: Marchantia ${s.label} → ${hubName}`, "--method", `${minimap2(s.preset ?? "asm5")} ${basename(hubGenome)} ${basename(s.genome)} (target ${hubName}, query ${s.label}). ${pafFilter}`, "--from-report", s.report, "--to-report", hubReport, "--fasta", s.genome, "--fasta", hubGenome, work(`${s.name}_to_${hubName}.paf`)],
    },
    {
      name: `mp_hub_to_${s.name}`,
      group,
      prepare: [align(s.genome, hubGenome, work(`${hubName}_to_${s.name}.paf`), s.preset)],
      args: ["--label", `minimap2 alignment: Marchantia ${hubName} → ${s.label}`, "--method", `${minimap2(s.preset ?? "asm5")} ${basename(s.genome)} ${basename(hubGenome)} (target ${s.label}, query ${hubName}). ${pafFilter}`, "--from-report", hubReport, "--to-report", s.report, "--fasta", hubGenome, "--fasta", s.genome, work(`${hubName}_to_${s.name}.paf`)],
    },
  ]);
}

const STORE_LIST: Store[] = [
  // Human, GRCh38 (annotated)
  { name: "human", group: "human", args: ["--label", "Human RefSeq annotation (GRCh38.p14)", "--assembly-report", r38, "--fasta", g38, "--fasta", ncbi(GRCH38, "rna.fna.gz"), "--fasta", ncbi(GRCH38, "protein.faa.gz"), ncbi(GRCH38, "genomic.gff.gz")] },
  { name: "human_rna", group: "human", args: ["--label", "Human RefSeq RNA", ncbi(GRCH38, "rna.gbff.gz")] },
  { name: "human_ensembl", group: "human", args: ["--label", "Human Ensembl release 116 (GRCh38)", "--seqid-map", r38, "--fasta", g38, "--fasta", raw(`${ENSEMBL}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz`), raw(`${ENSEMBL}/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz`)] },
  { name: "human_uniprot", group: "human", args: ["--label", "Human UniProt reference proteome (UP000005640)", ...upHuman] },
  { name: "human_mane", group: "human", args: ["--label", "MANE v1.5 (human)", "--taxon", "9606", "--organism", "Homo sapiens", raw(`${MANE}/MANE.GRCh38.v1.5.summary.txt.gz`), raw(`${MANE}/MANE.GRCh38.v1.5.refseq_rna.fna.gz`), raw(`${MANE}/MANE.GRCh38.v1.5.ensembl_rna.fna.gz`)] },
  { name: "grch38_names", group: "human", args: ["--label", "GRCh38.p14 sequence names (NCBI assembly report)", r38] },
  // T2: UniProt entries without an identical annotated protein, aligned to the proteins their ID mapping names
  {
    name: "human_uniprot_alignments",
    group: "human",
    args: [
      "--label", "UniProt to RefSeq / Ensembl protein alignments (human; entries without an identical protein)", "--taxon", "9606", "--organism", "Homo sapiens",
      ...[...upHuman, ncbi(GRCH38, "protein.faa.gz"), raw(`${ENSEMBL}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz`)].flatMap((f) => ["--fasta", f]),
      raw(`${IDMAPPING}/HUMAN_9606_idmapping_selected.tab.gz`),
    ],
  },
  // Human, GRCh37 and UCSC chains to / from GRCh38
  { name: "grch37", group: "grch37", args: ["--label", "GRCh37.p13 genome (NCBI assembly report and sequences)", "--assembly-report", r37, r37, g37] },
  { name: "chain_hg19ToHg38", group: "grch37", args: ["--label", "UCSC liftOver chains hg19 → hg38 (GRCh37 to GRCh38)", "--from-report", r37, "--to-report", r38, "--fasta", g37, "--fasta", g38, raw(`${UCSC}/hg19/liftOver/hg19ToHg38.over.chain.gz`)] },
  { name: "chain_hg38ToHg19", group: "grch37", args: ["--label", "UCSC liftOver chains hg38 → hg19 (GRCh38 to GRCh37)", "--from-report", r38, "--to-report", r37, "--fasta", g38, "--fasta", g37, raw(`${UCSC}/hg38/liftOver/hg38ToHg19.over.chain.gz`)] },
  // Human, T2T-CHM13 (a complete assembly of another individual, annotated) and UCSC chains to / from GRCh38
  { name: "chm13", group: "chm13", args: ["--label", "T2T-CHM13v2.0 RefSeq annotation (GCF_009914755.1)", "--assembly-report", rChm13, "--fasta", gChm13, "--fasta", ncbi(CHM13, "rna.fna.gz"), "--fasta", ncbi(CHM13, "protein.faa.gz"), ncbi(CHM13, "genomic.gff.gz")] },
  { name: "chain_hs1ToHg38", group: "chm13", args: ["--label", "UCSC liftOver chains hs1 → hg38 (T2T-CHM13v2.0 to GRCh38)", "--from-report", rChm13, "--to-report", r38, "--fasta", gChm13, "--fasta", g38, raw(`${UCSC}/hs1/liftOver/hs1ToHg38.over.chain.gz`)] },
  { name: "chain_hg38ToHs1", group: "chm13", args: ["--label", "UCSC liftOver chains hg38 → hs1 (GRCh38 to T2T-CHM13v2.0)", "--from-report", r38, "--to-report", rChm13, "--fasta", g38, "--fasta", gChm13, raw(`${UCSC}/hg38/liftOver/hg38ToHs1.over.chain.gz`)] },
  // Mouse, GRCm39 (annotated)
  { name: "mouse", group: "mouse", args: ["--label", "Mouse RefSeq annotation (GRCm39)", "--assembly-report", rm39, "--fasta", gm39, "--fasta", ncbi(GRCM39, "rna.fna.gz"), "--fasta", ncbi(GRCM39, "protein.faa.gz"), ncbi(GRCM39, "genomic.gff.gz")] },
  { name: "mouse_rna", group: "mouse", args: ["--label", "Mouse RefSeq RNA", ncbi(GRCM39, "rna.gbff.gz")] },
  { name: "mouse_uniprot", group: "mouse", args: ["--label", "Mouse UniProt reference proteome (UP000000589)", ...upMouse] },
  { name: "grcm39_names", group: "mouse", args: ["--label", "GRCm39 sequence names (NCBI assembly report)", rm39] },
  {
    name: "mouse_uniprot_alignments",
    group: "mouse",
    args: [
      "--label", "UniProt to RefSeq protein alignments (mouse; entries without an identical protein)", "--taxon", "10090", "--organism", "Mus musculus",
      ...[...upMouse, ncbi(GRCM39, "protein.faa.gz")].flatMap((f) => ["--fasta", f]),
      raw(`${IDMAPPING}/MOUSE_10090_idmapping_selected.tab.gz`),
    ],
  },
  // Mouse, GRCm38 and UCSC chains to / from GRCm39
  { name: "grcm38", group: "grcm38", args: ["--label", "GRCm38.p6 genome (NCBI assembly report and sequences)", "--assembly-report", rm38, rm38, gm38] },
  { name: "chain_mm10ToMm39", group: "grcm38", args: ["--label", "UCSC liftOver chains mm10 → mm39 (GRCm38 to GRCm39)", "--from-report", rm38, "--to-report", rm39, "--fasta", gm38, "--fasta", gm39, raw(`${UCSC}/mm10/liftOver/mm10ToMm39.over.chain.gz`)] },
  { name: "chain_mm39ToMm10", group: "grcm38", args: ["--label", "UCSC liftOver chains mm39 → mm10 (GRCm39 to GRCm38)", "--from-report", rm39, "--to-report", rm38, "--fasta", gm39, "--fasta", gm38, raw(`${UCSC}/mm39/liftOver/mm39ToMm10.over.chain.gz`)] },
  // Human <-> mouse
  { name: "chain_hg38ToMm39", group: "human_mouse", args: ["--label", "UCSC liftOver chains hg38 → mm39 (human to mouse)", "--from-report", r38, "--to-report", rm39, "--fasta", g38, "--fasta", gm39, raw(`${UCSC}/hg38/liftOver/hg38ToMm39.over.chain.gz`)] },
  { name: "chain_mm39ToHg38", group: "human_mouse", args: ["--label", "UCSC liftOver chains mm39 → hg38 (mouse to human)", "--from-report", rm39, "--to-report", r38, "--fasta", gm39, "--fasta", g38, raw(`${UCSC}/mm39/liftOver/mm39ToHg38.over.chain.gz`)] },
  // Arabidopsis
  { name: "arabidopsis", group: "arabidopsis", args: ["--label", "Arabidopsis RefSeq annotation (TAIR10.1)", "--assembly-report", rTair, "--fasta", gTair, "--fasta", ncbi(TAIR, "rna.fna.gz"), "--fasta", ncbi(TAIR, "protein.faa.gz"), ncbi(TAIR, "genomic.gff.gz")] },
  { name: "arabidopsis_rna", group: "arabidopsis", args: ["--label", "Arabidopsis RefSeq RNA", ncbi(TAIR, "rna.gbff.gz")] },
  { name: "arabidopsis_uniprot", group: "arabidopsis", args: ["--label", "Arabidopsis UniProt reference proteome (UP000006548)", ...upArab] },
  {
    name: "arabidopsis_uniprot_alignments",
    group: "arabidopsis",
    args: [
      "--label", "UniProt to RefSeq protein alignments (Arabidopsis; entries without an identical protein)", "--taxon", "3702", "--organism", "Arabidopsis thaliana",
      ...[...upArab, ncbi(TAIR, "protein.faa.gz")].flatMap((f) => ["--fasta", f]),
      raw(`${IDMAPPING}/ARATH_3702_idmapping_selected.tab.gz`),
    ],
  },
  // Arabidopsis TAIR10 (the previous RefSeq version): the same nuclear chromosomes and chloroplast; another
  // mitochondrial genome. UCSC GenArk has a chain TAIR10.1 -> TAIR10 only; the other direction is aligned here.
  { name: "tair10", group: "tair10", args: ["--label", "TAIR10 genome (GCF_000001735.3; NCBI assembly report and sequences)", "--assembly-report", rTair10, rTair10, gTair10] },
  { name: "chain_tair10.1ToTair10", group: "tair10", args: ["--label", "UCSC GenArk liftOver chains TAIR10.1 → TAIR10", "--from-report", rTair, "--to-report", rTair10, "--fasta", gTair, "--fasta", gTair10, raw("https://hgdownload.soe.ucsc.edu/hubs/GCF/000/001/735/GCF_000001735.4/liftOver/GCF_000001735.4_TAIR10.1ToGCF_000001735.3_TAIR10.over.chain.gz")] },
  {
    name: "tair10_to_tair10.1",
    group: "tair10",
    prepare: [align(gTair, gTair10, pafTair10)],
    args: ["--label", "minimap2 alignment: TAIR10 → TAIR10.1", "--method", `${minimap2("asm5")} ${basename(gTair)} ${basename(gTair10)} (target TAIR10.1, query TAIR10). ${pafFilter}; sequences of both assemblies skipped (identity)`, "--from-report", rTair10, "--to-report", rTair, "--fasta", gTair10, "--fasta", gTair, pafTair10],
  },
  // Marchantia: v7.1 (default) and v3.1 (INSDC only), joined by minimap2 alignments
  { name: "marchantia_v71", group: "marchantia", args: ["--label", "Marchantia polymorpha MpTak_v7.1 INSDC annotation (GCA_039105155.1)", "--assembly-report", rMp71, ncbi(MP71, "genomic.gbff.gz")] },
  { name: "marchantia", group: "marchantia", args: ["--label", "Marchantia polymorpha INSDC annotation (GCA_003032435.1)", "--assembly-report", rMp31, ncbi(MP31, "genomic.gbff.gz")] },
  { name: "marchantia_uniprot", group: "marchantia", args: ["--label", "Marchantia UniProt proteome (UP000244005)", ...upMarchantia] },
  {
    name: "mp_v31_to_v71",
    group: "marchantia",
    prepare: [align(gMp71, gMp31, paf31to71)],
    args: ["--label", "minimap2 alignment: Marchantia MpTak v3.1 → v7.1", "--method", `${minimap2("asm5")} ${basename(gMp71)} ${basename(gMp31)} (target v7.1, query v3.1). ${pafFilter}`, "--from-report", rMp31, "--to-report", rMp71, "--fasta", gMp31, "--fasta", gMp71, paf31to71],
  },
  {
    name: "mp_v71_to_v31",
    group: "marchantia",
    prepare: [align(gMp31, gMp71, paf71to31)],
    args: ["--label", "minimap2 alignment: Marchantia MpTak v7.1 → v3.1", "--method", `${minimap2("asm5")} ${basename(gMp31)} ${basename(gMp71)} (target v3.1, query v7.1). ${pafFilter}`, "--from-report", rMp71, "--to-report", rMp31, "--fasta", gMp71, "--fasta", gMp31, paf71to31],
  },
  // Marchantia: the other assemblies are aligned to MpTak_v7.1, the annotated one (a star, spec-ingest §20)
  { name: "marchantia_tak2", group: "marchantia", args: ["--label", "Marchantia polymorpha MpTak2_v7.1 (Tak-2) INSDC annotation (GCA_037833965.1)", "--assembly-report", rMpTak2, ncbi(MPTAK2, "genomic.gbff.gz")] },
  { name: "marchantia_v51", group: "marchantia", args: ["--label", "Marchantia polymorpha v5.1 INSDC annotation (GCA_009936355.2)", "--assembly-report", rMp51, ncbi(MP51, "genomic.gbff.gz")] },
  { name: "marchantia_cmv12", group: "marchantia", args: ["--label", "Marchantia polymorpha cmMarPoly1.2 genome, no annotation (GCA_965642975.2)", "--assembly-report", rMpCm, rMpCm, gMpCm] },
  { name: "marchantia_mpv4", group: "marchantia", args: ["--label", "Marchantia polymorpha Mp_v4 INSDC annotation (GCA_001641455.1)", "--assembly-report", rMpV4, ncbi(MPV4, "genomic.gbff.gz")] },
  { name: "marchantia_uniprot_v4", group: "marchantia", args: ["--label", "Marchantia UniProt proteome of Mp_v4 (UP000077202)", ...upMarchantiaV4] },
  ...starAlignments("marchantia", "MpTak_v7.1", rMp71, gMp71, [
    { name: "tak2", label: "MpTak2_v7.1 (Tak-2)", report: rMpTak2, genome: gMpTak2 },
    { name: "v51", label: "v5.1", report: rMp51, genome: gMp51 },
    // Another accession, 2-3% divergent and structurally different: with asm5 (~0.1%) only 40% of it is lifted,
    // with asm10 (~1%) 55%; asm20 (~5%) is the fitting preset (the alignment covers 60% / 71% of the genome).
    { name: "cmv12", label: "cmMarPoly1.2", report: rMpCm, genome: gMpCm, preset: "asm20" },
    // Scaffold-level assembly of pooled Tak-1 and Tak-2 (Oxford 2016); UniProt's UP000077202 is built on it.
    { name: "v4", label: "Mp_v4", report: rMpV4, genome: gMpV4 },
  ]),
  // Structures of the loaded proteomes
  {
    name: "sifts",
    group: "sifts",
    args: [
      "--label", "SIFTS: UniProt to PDB chains (loaded proteomes)", "--sifts-known-only",
      ...[...upHuman, ...upMouse, ...upArab, ...upMarchantia, ...upMarchantiaV4].flatMap((f) => ["--fasta", f]),
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
