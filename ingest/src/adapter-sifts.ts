// SIFTS adapter (enrichment, tier T1): UniProt <-> PDB chain residue correspondences.
// Input: uniprot_segments_observed.tsv(.gz) from https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/
//   PDB  CHAIN  SP_PRIMARY  RES_BEG  RES_END  PDB_BEG  PDB_END  SP_BEG  SP_END
// RES_* index the chain's SEQRES (label_seq_id, 1-based), our coordinate for PDB chains (design §3.7); PDB_* are author
// numbers (with insertion codes, kept as attributes for display); SP_* are UniProt residues. CHAIN is the author chain ID.
import { NamespaceRegistry, residueBlock, type Block } from "@togocoord/core";
import { ownString } from "./common.ts";
import type { Edge, Provenance, Sink } from "./model.ts";
import type { SequenceSource } from "./sequence.ts";
import { readLines } from "./stream.ts";

export interface SiftsOptions {
  file?: string;
  registry?: NamespaceRegistry;
  /** UniProt and PDB SEQRES residues (e.g. UniProt FASTA and pdb_seqres.txt) for validation. */
  source?: SequenceSource;
  /** Keep only rows whose UniProt accession passes (e.g. accessions of one proteome). */
  accept?: (uniprot: string) => boolean;
  /** Identity below which an edge is reported as a mismatch (misaligned rather than mutated); default 0.5. */
  minIdentity?: number;
}

interface Group {
  uniprot: string;
  pdb: string;
  blocks: Block[];
  author: string[];
}

export async function ingestSiftsFile(path: string, sink: Sink, options: SiftsOptions = {}): Promise<{ rows: number; edges: number; skipped: number }> {
  const registry = options.registry ?? new NamespaceRegistry();
  const minIdentity = options.minIdentity ?? 0.5;
  const provenance: Provenance = { adapter: "sifts", ...(options.file && { file: options.file }) };
  const seen = new Set<string>();
  let rows = 0;
  let edges = 0;
  let skipped = 0;
  let group: Group | undefined;

  const key = (ns: string, acc: string) => {
    try {
      return registry.refKey(ns, acc);
    } catch {
      return undefined;
    }
  };
  const sequence = (ref: string) => {
    if (seen.has(ref)) return;
    seen.add(ref);
    const length = options.source?.length?.(ref);
    sink.sequence({ ref, moltype: "protein", unit: "aa", length: length ?? 0, provenance: { ...provenance, record: ref } });
  };
  const flush = () => {
    if (!group) return;
    sequence(group.uniprot);
    sequence(group.pdb);
    const edge: Edge = {
      kind: "alignment",
      from: group.uniprot,
      to: group.pdb,
      blocks: group.blocks,
      attributes: { type: "SIFTS", authorNumbering: group.author.join(",") },
      provenance: { ...provenance, record: group.pdb },
      validation: validate(group.blocks, options.source, minIdentity),
    };
    sink.edge(edge);
    edges++;
    group = undefined;
  };

  let header: string[] | undefined;
  for await (const line of readLines(path)) {
    if (line.startsWith("#") || line === "") continue;
    const cols = line.split("\t");
    if (!header) {
      header = cols;
      continue;
    }
    const row = Object.fromEntries(header.map((h, i) => [h, cols[i] ?? ""]));
    rows++;
    const accession = row.SP_PRIMARY!;
    if (options.accept && !options.accept(accession)) continue;
    const uniprot = key("uniprot", accession);
    const pdb = key("pdb", `${row.PDB}.${row.CHAIN}`);
    const [resBeg, resEnd, spBeg, spEnd] = [row.RES_BEG, row.RES_END, row.SP_BEG, row.SP_END].map(Number) as [number, number, number, number];
    const length = resEnd - resBeg + 1;
    if (!uniprot || !pdb || !(length > 0) || spEnd - spBeg + 1 !== length || resBeg < 1 || spBeg < 1) {
      skipped++;
      if (skipped <= 5) sink.warning(`${options.file ?? path}: unusable SIFTS row: ${line}`);
      continue;
    }
    if (!group || group.uniprot !== uniprot || group.pdb !== pdb) {
      flush();
      group = { uniprot: ownString(uniprot), pdb: ownString(pdb), blocks: [], author: [] };
    }
    group.blocks.push(residueBlock({ srcRef: group.uniprot, srcBegin: spBeg, tgtRef: group.pdb, tgtBegin: resBeg, length }));
    group.author.push(`${row.PDB_BEG}-${row.PDB_END}`);
  }
  flush();
  if (skipped > 5) sink.warning(`${options.file ?? path}: ${skipped} unusable SIFTS rows`);
  return { rows, edges, skipped };
}

/** Residue identity over the aligned segments (engineered mutations are expected; misalignment gives ~5%). */
function validate(blocks: Block[], source: SequenceSource | undefined, minIdentity: number): Edge["validation"] {
  if (!source) return { status: "skipped", detail: "sequences not available" };
  let same = 0;
  let total = 0;
  for (const b of blocks) {
    const a = source.get(b.srcRef, b.src / 3, (b.src + b.len) / 3);
    const t = source.get(b.tgtRef, b.tgt / 3, (b.tgt + b.len) / 3);
    if (a === undefined || t === undefined) return { status: "skipped", detail: "sequences not available" };
    for (let i = 0; i < a.length; i++) if (a[i] === t[i]) same++;
    total += a.length;
  }
  const identity = total ? same / total : 0;
  return {
    status: identity >= minIdentity ? "ok" : "mismatch",
    basis: "full",
    detail: `${same}/${total} aligned residues identical`,
  };
}
