// BED adapter (enrichment): genome regions such as cis-regulatory elements (e.g. fanta.bio CREs), stored as
// annotations of the genome sequences. Format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
//   chrom start end [name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts] [extra...]
// Coordinates are 0-based half-open. Chromosome names (chr1, 1, ...) are read through the assembly report.
import { createContext, formatLocationId, NamespaceRegistry, type Location, type Segment } from "@togocoord/core";
import { ownString } from "./common.ts";
import type { Provenance, Sink } from "./model.ts";
import { readLines } from "./stream.ts";

export interface BedOptions {
  file?: string;
  registry?: NamespaceRegistry;
  /** Chromosome name -> sequence key (see assemblyReportSeqids). */
  refOf: (chrom: string) => string | undefined;
  /** Annotation type (default "region"), e.g. "CRE". */
  type?: string;
  /**
   * Names of the columns after the standard ones (after column 12, or after the last standard column the file uses:
   * BED9+2 means `--bed-columns Name,attributes`). "attributes" parses `key:value|key:value`.
   */
  extraColumns?: string[];
  /** Number of standard columns (3-12); default: 12, or fewer when extra columns follow (BED9+2: 9). */
  standardColumns?: number;
}

export interface BedStats {
  records: number;
  skipped: number;
}

export async function ingestBedFile(path: string, sink: Sink, options: BedOptions): Promise<BedStats> {
  const registry = options.registry ?? new NamespaceRegistry();
  const ctx = createContext({ registry, units: () => "nt" });
  const provenance: Provenance = { adapter: "bed", ...(options.file && { file: options.file }) };
  const extra = options.extraColumns ?? [];
  const stats: BedStats = { records: 0, skipped: 0 };
  const missing = new Set<string>();
  for await (const line of readLines(path)) {
    if (!line || line.startsWith("#") || line.startsWith("track") || line.startsWith("browser")) continue;
    const f = line.split("\t");
    const standard = Math.min(options.standardColumns ?? (extra.length ? f.length - extra.length : 12), 12);
    const [chrom, s, e, name, score, strand, thickStart, thickEnd, , blockCount, blockSizes, blockStarts] = f.slice(0, standard);
    const ref = chrom ? options.refOf(chrom) : undefined;
    if (!ref) {
      stats.skipped++;
      if (chrom && !missing.has(chrom)) {
        missing.add(chrom);
        sink.warning(`${options.file ?? path}: no sequence for ${chrom}; its regions are skipped`);
      }
      continue;
    }
    const start = Number(s);
    const end = Number(e);
    const minus = strand === "-";
    // BED12 blocks become the segments of a join, in the order of the strand.
    let spans: Array<[number, number]> = [[start, end]];
    if (blockCount && blockSizes && blockStarts && Number(blockCount) > 1) {
      const sizes = blockSizes.split(",").filter(Boolean).map(Number);
      const starts = blockStarts.split(",").filter(Boolean).map(Number);
      spans = starts.map((b, i) => [start + b, start + b + sizes[i]!]);
    }
    const segments: Segment[] = (minus ? [...spans].reverse() : spans).map(([a, b]) => ({ ref, start: a, end: b, strand: minus ? -1 : 1 }));
    const location: Location = { outer: ref, kind: "join", segments };
    const attributes: Record<string, string[]> = {};
    if (name) attributes.ID = [ownString(name)];
    if (score && score !== "0" && score !== ".") attributes.score = [score];
    if (thickStart && thickEnd && (Number(thickStart) !== start || Number(thickEnd) !== end) && Number(thickEnd) > Number(thickStart)) {
      attributes.thick = [formatLocationId({ outer: ref, kind: "join", segments: [{ ref, start: Number(thickStart), end: Number(thickEnd), strand: minus ? -1 : 1 }] }, ctx)];
    }
    extra.forEach((column, i) => {
      const v = f[standard + i];
      if (v === undefined || v === "" || v === ".") return;
      if (column === "attributes") {
        for (const kv of v.split("|")) {
          const at = kv.indexOf(":");
          if (at > 0) attributes[kv.slice(0, at)] = [ownString(kv.slice(at + 1))];
        }
      } else attributes[column] = [ownString(v)];
    });
    sink.annotation({
      location: formatLocationId(location, ctx),
      type: options.type ?? "region",
      attributes,
      extent: { ref, start, end },
      provenance: { ...provenance, ...(name && { record: ownString(name) }) },
    });
    stats.records++;
  }
  return stats;
}
