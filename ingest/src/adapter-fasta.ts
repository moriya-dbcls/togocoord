// FASTA adapter: sequence records with checksums (identity by refget digest; design §6, spec-service §3).
import { NamespaceRegistry, type Unit } from "@togocoord/core";
import { accessionRef, ownString } from "./common.ts";
import type { Provenance, SequenceRecord, Sink } from "./model.ts";
import { checksums } from "./sequence.ts";
import { readLines } from "./stream.ts";

export interface FastaOptions {
  file?: string;
  registry?: NamespaceRegistry;
  /** Header -> sequence key; the default understands UniProt, Ensembl, RefSeq and INSDC headers. */
  refOf?: (header: string) => string | undefined;
  /** Residue unit when the namespace does not imply one (default: from the key, else "aa"). */
  unit?: Unit;
}

/**
 * Sequence key from a FASTA header:
 *   `sp|P07203|GPX1_HUMAN ...` / `tr|...`       -> uniprot:P07203 (isoforms keep their suffix, e.g. P07203-2)
 *   `ENSP00000451042.1 pep ...`                  -> ensembl:ENSP00000451042.1
 *   `NP_000572.2 glutathione ...`, `AB000001.1`  -> refseq:/insdc:
 *   `101m_A mol:protein length:154  MYOGLOBIN`   -> pdb:101M.A (wwPDB pdb_seqres.txt; `mol:na` chains are skipped)
 */
export function defaultFastaRef(header: string, registry: NamespaceRegistry): string | undefined {
  const first = header.trim().split(/\s+/)[0] ?? "";
  const up = /^(?:sp|tr)\|([^|]+)\|/.exec(first);
  const tryKey = (ns: string, acc: string) => {
    try {
      return registry.refKey(ns, acc);
    } catch {
      return undefined;
    }
  };
  if (up) return tryKey("uniprot", up[1]!);
  const pdb = /^([0-9][A-Za-z0-9]{3})_(\S+)$/.exec(first);
  if (pdb) return /\bmol:protein\b/.test(header) ? tryKey("pdb", `${pdb[1]}.${pdb[2]}`) : undefined;
  if (/^ENS[A-Z]*[GTPER]\d{11}/.test(first)) return tryKey("ensembl", first);
  return accessionRef(first, registry);
}

/** Stream a (gzipped) FASTA file into sequence records. Residues are not stored, only length and checksums. */
export async function ingestFastaFile(path: string, sink: Sink, options: FastaOptions = {}): Promise<{ records: number; skipped: number }> {
  const registry = options.registry ?? new NamespaceRegistry();
  const refOf = options.refOf ?? ((h: string) => defaultFastaRef(h, registry));
  const provenance: Provenance = { adapter: "fasta", ...(options.file && { file: options.file }) };
  let header: string | undefined;
  let chunks: string[] = [];
  let records = 0;
  let skipped = 0;
  const flush = () => {
    if (header === undefined) return;
    const ref = refOf(header);
    if (!ref) {
      skipped++;
      if (skipped <= 5) sink.warning(`${options.file ?? path}: FASTA header '${header.slice(0, 60)}' has no recognised accession; skipped`);
    } else {
      const residues = chunks.join("").toUpperCase();
      const unit = options.unit ?? registry.defaultUnit(ref) ?? "aa";
      const transcript = /^(?:refseq:[NX][MR]_|ensembl:ENS[A-Z]*T\d)/.test(ref);
      const record: SequenceRecord = {
        ref: ownString(ref),
        moltype: unit === "aa" ? "protein" : transcript ? "RNA" : "DNA",
        unit,
        length: residues.length,
        provenance: { ...provenance, record: ownString(ref) },
        ...checksums(residues),
      };
      sink.sequence(record);
      records++;
    }
    chunks = [];
  };
  for await (const line of readLines(path)) {
    if (line.startsWith(">")) {
      flush();
      header = line.slice(1);
    } else if (header !== undefined) {
      chunks.push(line.trim());
    }
  }
  flush();
  if (skipped > 5) sink.warning(`${options.file ?? path}: ${skipped} FASTA records without a recognised accession`);
  return { records, skipped };
}
