// Streaming readers for large inputs (scaling §3): line by line, gzip-aware, constant memory per feature/record.
import { createReadStream, writeSync } from "node:fs";
import { StringDecoder } from "node:string_decoder";
import { createGunzip } from "node:zlib";
import { GenBankIngestor } from "./adapter-gbff.ts";
import { Gff3Ingestor, type Gff3Options } from "./adapter-gff3.ts";
import type { AdapterOptions } from "./common.ts";
import { parseGenBankRecord } from "./gbff.ts";
import { FeatureGrouper, parseGffLine, sequenceRegion } from "./gff3.ts";
import type { Annotation, Edge, SequenceRecord, Sink } from "./model.ts";

/**
 * Lines of a (gzipped) text file. Chunks are pulled only after the previous chunk's lines were consumed, so memory
 * stays bounded when parsing is slower than decompression (readline's async iterator buffers without limit there).
 */
export async function* readLines(path: string): AsyncGenerator<string> {
  const raw = createReadStream(path, { highWaterMark: 1 << 20 });
  const input = path.endsWith(".gz") ? raw.pipe(createGunzip({ chunkSize: 1 << 20 })) : raw;
  const decoder = new StringDecoder("utf8");
  let rest = "";
  for await (const chunk of input as AsyncIterable<Buffer>) {
    const lines = (rest + decoder.write(chunk)).split("\n");
    rest = lines.pop()!;
    for (const line of lines) yield line.endsWith("\r") ? line.slice(0, -1) : line;
  }
  rest += decoder.end();
  if (rest) yield rest.endsWith("\r") ? rest.slice(0, -1) : rest;
}

export interface StreamStats {
  lines: number;
  features: number;
}

/**
 * Stream a GFF3 file into a sink. Rows of one feature must be within `window` features of each other
 * (NCBI RefSeq GFF3 keeps them adjacent); violations are reported as warnings. A `##FASTA` section is not read:
 * pass residues through `options.source` (e.g. a FaiSequenceSource).
 */
export async function ingestGff3File(path: string, sink: Sink, options: Gff3Options & { window?: number } = {}): Promise<StreamStats> {
  const ingestor = new Gff3Ingestor(sink, options);
  const grouper = new FeatureGrouper(options.window ?? 1000);
  let lines = 0;
  let features = 0;
  const emit = (fs: ReturnType<FeatureGrouper["push"]>) => {
    for (const f of fs) {
      ingestor.feature(f);
      features++;
    }
  };
  for await (const line of readLines(path)) {
    lines++;
    if (line.startsWith("#")) {
      const region = sequenceRegion(line);
      if (region) ingestor.sequenceRegion(region[0], region[1]);
      if (line.startsWith("##FASTA")) {
        sink.warning(`${path}: ##FASTA section ignored when streaming; pass sequences with --fasta`);
        break;
      }
      continue;
    }
    const row = parseGffLine(line, lines);
    if (row) emit(grouper.push(row));
  }
  emit(grouper.end());
  for (const v of grouper.violations) sink.warning(`${path}: ${v}`);
  ingestor.finish();
  return { lines, features };
}

/** Stream a GenBank/GenPept file record by record. */
export async function ingestGenBankFile(path: string, sink: Sink, options: AdapterOptions = {}): Promise<StreamStats> {
  const ingestor = new GenBankIngestor(sink, options);
  let lines = 0;
  let features = 0;
  let record: string[] = [];
  const flush = () => {
    if (record.some((l) => l.startsWith("LOCUS"))) {
      const r = parseGenBankRecord(record);
      features += r.features.length;
      ingestor.record(r);
    }
    record = [];
  };
  for await (const line of readLines(path)) {
    lines++;
    if (line.startsWith("//")) flush();
    else record.push(line);
  }
  flush();
  return { lines, features };
}

/** Writes JSON Lines to a file descriptor with buffered synchronous writes (no unbounded stream buffering). */
export class JsonlSink implements Sink {
  readonly #fd: number;
  #buffer: string[] = [];
  #size = 0;
  readonly counts = { sequence: 0, edge: 0, annotation: 0, warning: 0 };
  readonly validation: Record<string, number> = {};
  readonly mismatches: string[] = [];

  constructor(fd = 1) {
    this.#fd = fd;
  }

  #write(obj: object): void {
    const line = JSON.stringify(obj) + "\n";
    this.#buffer.push(line);
    this.#size += line.length;
    if (this.#size > 1 << 20) this.flush();
  }

  flush(): void {
    if (this.#buffer.length) writeSync(this.#fd, this.#buffer.join(""));
    this.#buffer = [];
    this.#size = 0;
  }

  sequence(record: SequenceRecord): void {
    this.counts.sequence++;
    this.#write({ record: "sequence", ...record });
  }
  edge(edge: Edge): void {
    this.counts.edge++;
    this.validation[edge.validation.status] = (this.validation[edge.validation.status] ?? 0) + 1;
    if (edge.validation.status === "mismatch") this.mismatches.push(`${edge.from} -> ${edge.to}: ${edge.validation.detail}`);
    this.#write({ record: "edge", ...edge });
  }
  annotation(annotation: Annotation): void {
    this.counts.annotation++;
    this.#write({ record: "annotation", ...annotation });
  }
  warning(message: string): void {
    this.counts.warning++;
    this.#write({ record: "warning", message });
  }
}
