import { readFileSync } from "node:fs";
import { parseFasta } from "../src/fasta.ts";
import { parseGenBank } from "../src/gbff.ts";
import { MemorySequenceSource } from "../src/sequence.ts";

export const fixture = (name: string) => readFileSync(new URL(`./fixtures/${name}`, import.meta.url), "utf8");

/** Residues of every GenBank fixture record, keyed by refseq:/insdc: accession. */
export function genbankSource(...names: string[]): MemorySequenceSource {
  const source = new MemorySequenceSource();
  for (const name of names) {
    for (const r of parseGenBank(fixture(name))) {
      if (r.sequence && r.version) source.add(`${/^[A-Z]{2}_/.test(r.version) ? "refseq" : "insdc"}:${r.version}`, r.sequence);
    }
  }
  return source;
}

/** Sequences for the chr3 cDNA_match fixture: transcripts plus two genome slices at their offsets. */
export function chr3Source(): MemorySequenceSource {
  const source = new MemorySequenceSource();
  for (const name of ["NM_012234.7.fa", "NR_047574.1.fa"]) {
    for (const [id, s] of parseFasta(fixture(name))) source.add(`refseq:${id}`, s);
  }
  for (const name of ["NC_000003.12_72374597-72446623.fa", "NC_000003.12_194584268-194590832.fa"]) {
    const offset = Number(/_(\d+)-/.exec(name)![1]) - 1;
    for (const [, s] of parseFasta(fixture(name))) source.add("refseq:NC_000003.12", s, offset);
  }
  return source;
}
