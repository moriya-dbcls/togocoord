// Block-list mappings and their algebra (spec-core §4).
import { MappingError } from "./errors.ts";
import type { Location } from "./location.ts";

/** Ungapped one-to-one correspondence of `len` units. */
export interface Block {
  readonly srcRef: string;
  readonly src: number;
  readonly tgtRef: string;
  readonly tgt: number;
  readonly len: number;
  readonly rev: boolean;
}

export interface BlockHit {
  block: Block;
  /** Registration order within the mapping; used as a deterministic tie-breaker. */
  index: number;
}

interface RefIndex {
  /** Block indices sorted by `src`. */
  order: number[];
  maxLen: number;
}

/** Immutable set of blocks, indexed by source reference. Overlaps are allowed on both sides. */
export class Mapping {
  readonly blocks: readonly Block[];
  readonly #index = new Map<string, RefIndex>();

  constructor(blocks: Iterable<Block>) {
    this.blocks = Object.freeze([...blocks].map(validateBlock));
    this.blocks.forEach((b, i) => {
      let entry = this.#index.get(b.srcRef);
      if (!entry) this.#index.set(b.srcRef, (entry = { order: [], maxLen: 0 }));
      entry.order.push(i);
      entry.maxLen = Math.max(entry.maxLen, b.len);
    });
    for (const entry of this.#index.values()) {
      entry.order.sort((x, y) => this.blocks[x]!.src - this.blocks[y]!.src || x - y);
    }
  }

  /** Blocks whose source interval intersects `[start, end)` (or contains unit `start` when empty). */
  overlapping(ref: string, start: number, end: number): BlockHit[] {
    const entry = this.#index.get(ref);
    if (!entry) return [];
    const stop = Math.max(end, start + 1);
    // Any overlapping block starts after `start - maxLen`.
    const floor = start - entry.maxLen;
    let lo = 0;
    let hi = entry.order.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (this.blocks[entry.order[mid]!]!.src <= floor) lo = mid + 1;
      else hi = mid;
    }
    const hits: BlockHit[] = [];
    for (let k = lo; k < entry.order.length; k++) {
      const index = entry.order[k]!;
      const block = this.blocks[index]!;
      if (block.src >= stop) break;
      if (block.src + block.len > start) hits.push({ block, index });
    }
    return hits.sort((a, b) => a.index - b.index);
  }

  inverse(): Mapping {
    return invert(this);
  }

  /** `this` (a->b) followed by `next` (b->c). */
  then(next: Mapping): Mapping {
    return compose(this, next);
  }

  concat(...others: Mapping[]): Mapping {
    return new Mapping([this, ...others].flatMap((m) => m.blocks));
  }
}

function validateBlock(b: Block): Block {
  const ints = [b.src, b.tgt, b.len];
  if (!ints.every(Number.isSafeInteger) || b.src < 0 || b.tgt < 0 || b.len <= 0) {
    throw new MappingError(`invalid block ${JSON.stringify(b)}`);
  }
  return Object.freeze({ srcRef: b.srcRef, src: b.src, tgtRef: b.tgtRef, tgt: b.tgt, len: b.len, rev: b.rev });
}

/** Target interval of the source sub-interval `[start, end)` of `b` (must lie inside the block). */
export function projectInterval(b: Block, start: number, end: number): [number, number] {
  return b.rev
    ? [b.tgt + b.len - (end - b.src), b.tgt + b.len - (start - b.src)]
    : [b.tgt + (start - b.src), b.tgt + (end - b.src)];
}

export function invert(m: Mapping): Mapping {
  return new Mapping(m.blocks.map((b) => ({ srcRef: b.tgtRef, src: b.tgt, tgtRef: b.srcRef, tgt: b.src, len: b.len, rev: b.rev })));
}

export function compose(ab: Mapping, bc: Mapping): Mapping {
  const out: Block[] = [];
  for (const a of ab.blocks) {
    for (const { block: b } of bc.overlapping(a.tgtRef, a.tgt, a.tgt + a.len)) {
      const lo = Math.max(a.tgt, b.src);
      const hi = Math.min(a.tgt + a.len, b.src + b.len);
      out.push({
        srcRef: a.srcRef,
        src: a.rev ? a.src + (a.tgt + a.len - hi) : a.src + (lo - a.tgt),
        tgtRef: b.tgtRef,
        tgt: b.rev ? b.tgt + (b.src + b.len - hi) : b.tgt + (lo - b.src),
        len: hi - lo,
        rev: a.rev !== b.rev,
      });
    }
  }
  return new Mapping(out);
}

/** Mapping from the feature sequence `featureRef` (1..L along `loc`) onto the referenced sequences. */
export function mappingFromLocation(featureRef: string, loc: Location): Mapping {
  const blocks: Block[] = [];
  let offset = 0;
  for (const seg of loc.segments) {
    const len = seg.end - seg.start;
    if (len === 0) throw new MappingError("a feature location cannot contain between-positions");
    blocks.push({ srcRef: featureRef, src: offset, tgtRef: seg.ref, tgt: seg.start, len, rev: seg.strand === -1 });
    offset += len;
  }
  return new Mapping(blocks);
}

export interface CdsOptions {
  /** Protein sequence key (unit: aa). */
  protein: string;
  /** CDS location on nucleotide sequences, including the stop codon if present. */
  cds: Location;
  /** INSDC /codon_start (1..3). */
  codonStart?: number;
  /** Protein length in residues; supplied by the caller (spec-core §4.2). */
  aaLength: number;
  /**
   * The protein's first residue stands for the incomplete codon formed by the `codonStart - 1` bases before the first
   * complete codon (Ensembl writes it as `X`); INSDC /codon_start has no such residue. Requires codonStart > 1.
   */
  leadingPartialCodon?: boolean;
}

/** Protein (codon units) -> nucleotide mapping of a CDS. */
export function cdsMapping({ protein, cds, codonStart = 1, aaLength, leadingPartialCodon = false }: CdsOptions): Mapping {
  if (![1, 2, 3].includes(codonStart)) throw new MappingError(`codon_start must be 1, 2 or 3 (got ${codonStart})`);
  if (!Number.isSafeInteger(aaLength) || aaLength < 1) throw new MappingError(`invalid aaLength ${aaLength}`);
  if (leadingPartialCodon && codonStart === 1) throw new MappingError("leadingPartialCodon needs codon_start 2 or 3");
  const cdsLength = cds.segments.reduce((n, s) => n + s.end - s.start, 0);
  // Protein units [skip, 3·aaLength) map onto CDS bases starting at `tgt` (the partial first residue keeps only its
  // last codonStart-1 units).
  const phase = codonStart - 1;
  const skip = leadingPartialCodon ? 3 - phase : 0;
  const tgt = leadingPartialCodon ? 0 : phase;
  if (tgt + 3 * aaLength - skip > cdsLength) {
    throw new MappingError(`${aaLength} aa (codon_start ${codonStart}) do not fit in a ${cdsLength} nt CDS`);
  }
  const virtual = `${protein}#cds`;
  const scale = new Mapping([{ srcRef: protein, src: skip, tgtRef: virtual, tgt, len: 3 * aaLength - skip, rev: false }]);
  return compose(scale, mappingFromLocation(virtual, cds));
}

export interface ResidueBlockSpec {
  srcRef: string;
  /** 1-based first residue on the source protein. */
  srcBegin: number;
  tgtRef: string;
  /** 1-based first residue on the target protein. */
  tgtBegin: number;
  /** Number of residues. */
  length: number;
}

/** aa->aa block from residue numbering (converted to codon units). */
export function residueBlock(spec: ResidueBlockSpec): Block {
  return {
    srcRef: spec.srcRef,
    src: 3 * (spec.srcBegin - 1),
    tgtRef: spec.tgtRef,
    tgt: 3 * (spec.tgtBegin - 1),
    len: 3 * spec.length,
    rev: false,
  };
}
