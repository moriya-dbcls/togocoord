// MANE summary adapter (enrichment): tags the RefSeq and Ensembl transcripts and proteins of MANE Select and
// MANE Plus Clinical (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.*.summary.txt.gz).
// Only sequence records are written (tags, gene symbol); lengths stay unknown (0) and come from other stores.
import { NamespaceRegistry } from "@togocoord/core";
import { accessionRef, ownString } from "./common.ts";
import type { Provenance, SequenceRecord, Sink } from "./model.ts";
import { readLines } from "./stream.ts";

export async function ingestManeSummary(path: string, sink: Sink, options: { file?: string; registry?: NamespaceRegistry } = {}): Promise<{ genes: number; sequences: number }> {
  const registry = options.registry ?? new NamespaceRegistry();
  const provenance: Provenance = { adapter: "mane", ...(options.file && { file: options.file }) };
  let header: string[] | undefined;
  let genes = 0;
  let sequences = 0;
  for await (const line of readLines(path)) {
    if (!line) continue;
    const cols = line.split("\t");
    if (line.startsWith("#")) {
      header = cols.map((c) => c.replace(/^#/, ""));
      continue;
    }
    if (!header) throw new Error(`${path}: missing header line`);
    const row = Object.fromEntries(header.map((h, i) => [h, cols[i] ?? ""]));
    const status = row.MANE_status;
    if (!status) continue;
    genes++;
    const kinds: Array<[string, SequenceRecord["moltype"], "nt" | "aa"]> = [
      [row.RefSeq_nuc!, "RNA", "nt"],
      [row.RefSeq_prot!, "protein", "aa"],
      [row.Ensembl_nuc!, "RNA", "nt"],
      [row.Ensembl_prot!, "protein", "aa"],
    ];
    for (const [accession, moltype, unit] of kinds) {
      const ref = accession && accessionRef(accession, registry);
      if (!ref) continue;
      sink.sequence({
        ref: ownString(ref),
        moltype,
        unit,
        length: 0,
        tags: [status],
        ...(row.symbol && { gene: row.symbol }),
        provenance: { ...provenance, record: ownString(ref) },
      });
      sequences++;
    }
  }
  return { genes, sequences };
}
