// Naive per-unit reference implementation (design §10). Deliberately shares no logic with src/.
import type { Location, Mapping, MapResult, Segment } from "../src/index.ts";

export interface UnitImage {
  ref: string;
  pos: number;
  rev: boolean;
}

/** All images of a single unit, by linear scan over every block. */
export function mapUnit(m: Mapping, ref: string, u: number): UnitImage[] {
  const out: UnitImage[] = [];
  for (const b of m.blocks) {
    if (b.srcRef !== ref || u < b.src || u >= b.src + b.len) continue;
    const off = u - b.src;
    out.push({ ref: b.tgtRef, pos: b.rev ? b.tgt + b.len - 1 - off : b.tgt + off, rev: b.rev });
  }
  return out;
}

/** Units of a segment in traversal order. */
export function traversal(seg: Segment): number[] {
  const units: number[] = [];
  for (let u = seg.start; u < seg.end; u++) units.push(u);
  return seg.strand === 1 ? units : units.reverse();
}

const key = (srcRef: string, u: number, tgtRef: string, v: number, strand: number) => `${srcRef}|${u}|${tgtRef}|${v}|${strand}`;

/** Expected (source unit -> target unit, strand) pairs and unmapped source units for a non-between location. */
export function expectedPairs(loc: Location, m: Mapping): { pairs: string[]; unmapped: string[] } {
  const pairs: string[] = [];
  const unmapped: string[] = [];
  for (const seg of loc.segments) {
    for (const u of traversal(seg)) {
      const images = mapUnit(m, seg.ref, u);
      if (images.length === 0) unmapped.push(`${seg.ref}|${u}`);
      for (const img of images) pairs.push(key(seg.ref, u, img.ref, img.pos, img.rev ? -seg.strand : seg.strand));
    }
  }
  return { pairs: pairs.sort(), unmapped: unmapped.sort() };
}

/** The same pairs, reconstructed from a MapResult's pieces and unmapped location. */
export function actualPairs(result: MapResult): { pairs: string[]; unmapped: string[] } {
  const pairs: string[] = [];
  for (const { source, target } of result.pieces) {
    const s = traversal(source);
    const t = traversal(target);
    if (s.length !== t.length) throw new Error("piece length mismatch");
    s.forEach((u, i) => pairs.push(key(source.ref, u, target.ref, t[i]!, target.strand)));
  }
  const unmapped = (result.unmapped?.segments ?? []).flatMap((seg) => traversal(seg).map((u) => `${seg.ref}|${u}`));
  return { pairs: pairs.sort(), unmapped: unmapped.sort() };
}
