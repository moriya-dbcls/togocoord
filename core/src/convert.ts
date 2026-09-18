// Mapping a Location through a Mapping (spec-core §5).
import type { CoordContext, Location, Segment } from "./location.ts";
import { projectInterval, type Mapping } from "./mapping.ts";

export interface Piece {
  /** Index of the input segment this piece came from. */
  inputIndex: number;
  source: Segment;
  target: Segment;
}

export interface TargetLocation {
  location: Location;
  orientation: "forward" | "reverse" | "mixed";
}

export interface MapResult {
  /** Unmerged one-to-one pieces in traversal order (for correspondence tables). */
  pieces: Piece[];
  /** One location per target sequence, in order of first appearance. */
  targets: TargetLocation[];
  /** Source parts without a counterpart, or null. */
  unmapped: Location | null;
  /** A one-of (`a.b`) input was split into several target intervals. */
  uncertain: boolean;
}

interface Hit {
  a: number;
  b: number;
  ta: number;
  tb: number;
  tgtRef: string;
  strand: 1 | -1;
  index: number;
  /** Input segment index. */
  seg: number;
}

export function mapLocation(loc: Location, m: Mapping, ctx: CoordContext): MapResult {
  const pieces: Piece[] = [];
  const outputs: Segment[] = [];
  const unmapped: Segment[] = [];
  /** Per input segment and target sequence: is a source unit mapped onto that target? */
  const coverage: Array<(u: number, target: string) => boolean> = [];
  let run: Hit[] = [];
  let uncertain = false;

  const flush = () => {
    for (const g of mergeHits(run, loc)) {
      outputs.push(finishGroup(g, loc, coverage));
      const s = loc.segments[g.first.seg]!;
      if (s.uncertain && g.first.seg === g.last.seg && g.only) outputs.at(-1)!.uncertain = true;
    }
    // A one-of input split into several target intervals can no longer be written as `a.b`.
    for (const i of new Set(run.map((h) => h.seg))) {
      if (loc.segments[i]!.uncertain && mergeHits(run.filter((h) => h.seg === i), loc).length > 1) uncertain = true;
    }
    run = [];
  };

  loc.segments.forEach((seg, inputIndex) => {
    if (seg.start === seg.end) {
      flush();
      coverage[inputIndex] = () => false;
      const targets = mapBetween(seg, m);
      for (const target of targets) {
        pieces.push({ inputIndex, source: seg, target });
        outputs.push(target);
      }
      if (targets.length === 0) unmapped.push(seg);
      return;
    }

    const hits: Hit[] = m.overlapping(seg.ref, seg.start, seg.end).map(({ block, index }) => {
      const a = Math.max(seg.start, block.src);
      const b = Math.min(seg.end, block.src + block.len);
      const [ta, tb] = projectInterval(block, a, b);
      return { a, b, ta, tb, tgtRef: block.tgtRef, strand: block.rev ? flip(seg.strand) : seg.strand, index, seg: inputIndex };
    });
    hits.sort(seg.strand === 1 ? (x, y) => x.a - y.a || x.index - y.index : (x, y) => y.b - x.b || x.index - y.index);

    for (const h of hits) {
      pieces.push({
        inputIndex,
        source: { ref: seg.ref, start: h.a, end: h.b, strand: seg.strand },
        target: { ref: h.tgtRef, start: h.ta, end: h.tb, strand: h.strand },
      });
    }

    const covered = union(hits.map((h) => [h.a, h.b] as const));
    const byTarget = new Map<string, Array<[number, number]>>();
    for (const h of hits) {
      const list = byTarget.get(h.tgtRef);
      if (list) list.push([h.a, h.b]);
      else byTarget.set(h.tgtRef, [[h.a, h.b]]);
    }
    coverage[inputIndex] = (u, target) => (byTarget.get(target) ?? []).some(([x, y]) => x <= u && u < y);
    unmapped.push(...uncovered(seg, covered));
    run.push(...hits);
  });
  flush();

  return {
    pieces,
    targets: groupTargets(outputs, loc.kind, ctx),
    unmapped: unmapped.length ? { outer: loc.outer, kind: loc.kind, segments: unmapped } : null,
    uncertain,
  };
}

function flip(strand: 1 | -1): 1 | -1 {
  return strand === 1 ? -1 : 1;
}

interface Merged {
  tgtRef: string;
  strand: 1 | -1;
  ta: number;
  tb: number;
  /** First and last hit (traversal order); they carry the group's traversal start and end. */
  first: Hit;
  last: Hit;
  /** Source extent within the last segment, for the within-segment progression test. */
  srcLo: number;
  srcHi: number;
  /** The group holds every hit of its (single) segment. */
  only: boolean;
}

/**
 * Join hits that are consecutive for their target sequence and contiguous on it (spec-core §5.4). Within one input segment the source must
 * advance (a base read twice, as in ribosomal slippage, stays split); across segments of a `join` the hits may come
 * from different exons, which is how a CDS location maps to one protein interval. One-of segments are not merged
 * with others.
 */
function mergeHits(hits: Hit[], loc: Location): Merged[] {
  const out: Merged[] = [];
  // Hits for several target sequences interleave (e.g. overlapping transcripts); merge per target sequence.
  const open = new Map<string, Merged>();
  for (const h of hits) {
    const cur = open.get(h.tgtRef);
    const seg = loc.segments[h.seg]!;
    let joinable = false;
    if (cur && cur.tgtRef === h.tgtRef && cur.strand === h.strand && (h.strand === 1 ? h.ta === cur.tb : h.tb === cur.ta)) {
      if (cur.last.seg === h.seg) joinable = seg.strand === 1 ? h.a >= cur.srcHi : h.b <= cur.srcLo;
      else joinable = loc.kind === "join" && !seg.uncertain && !loc.segments[cur.last.seg]!.uncertain;
    }
    if (cur && joinable) {
      cur.ta = Math.min(cur.ta, h.ta);
      cur.tb = Math.max(cur.tb, h.tb);
      if (cur.last.seg === h.seg) {
        cur.srcLo = Math.min(cur.srcLo, h.a);
        cur.srcHi = Math.max(cur.srcHi, h.b);
      } else {
        cur.srcLo = h.a;
        cur.srcHi = h.b;
        cur.only = false;
      }
      cur.last = h;
    } else {
      const g: Merged = { tgtRef: h.tgtRef, strand: h.strand, ta: h.ta, tb: h.tb, first: h, last: h, srcLo: h.a, srcHi: h.b, only: true };
      out.push(g);
      open.set(h.tgtRef, g);
    }
  }
  // `only` also requires that no other group starts in the same segment.
  const perSegment = new Map<number, number>();
  for (const g of out) for (const s of new Set([g.first.seg, g.last.seg])) perSegment.set(s, (perSegment.get(s) ?? 0) + 1);
  for (const g of out) if ((perSegment.get(g.first.seg) ?? 0) > 1) g.only = false;
  return out;
}

/** Output segment of a merged group, with truncation marks (spec-core §5.5). */
function finishGroup(g: Merged, loc: Location, coverage: Array<(u: number, target: string) => boolean>): Segment {
  const startSeg = loc.segments[g.first.seg]!;
  const endSeg = loc.segments[g.last.seg]!;
  const covStart = (u: number) => coverage[g.first.seg]!(u, g.tgtRef);
  const covEnd = (u: number) => coverage[g.last.seg]!(u, g.tgtRef);
  // Truncated where the unit just outside the group (in traversal order) is part of the input but not mapped onto
  // this target sequence; other targets do not matter (results must not depend on unrelated edges).
  const startCut =
    startSeg.strand === 1
      ? g.first.a === startSeg.start ? !!startSeg.fuzzyLow : !covStart(g.first.a - 1)
      : g.first.b === startSeg.end ? !!startSeg.fuzzyHigh : !covStart(g.first.b);
  const endCut =
    endSeg.strand === 1
      ? g.last.b === endSeg.end ? !!endSeg.fuzzyHigh : !covEnd(g.last.b)
      : g.last.a === endSeg.start ? !!endSeg.fuzzyLow : !covEnd(g.last.a - 1);
  const [fuzzyLow, fuzzyHigh] = g.strand === 1 ? [startCut, endCut] : [endCut, startCut];
  const out: Segment = { ref: g.tgtRef, start: g.ta, end: g.tb, strand: g.strand };
  if (fuzzyLow) out.fuzzyLow = true;
  if (fuzzyHigh) out.fuzzyHigh = true;
  return out;
}

function union(intervals: ReadonlyArray<readonly [number, number]>): Array<[number, number]> {
  const sorted = [...intervals].sort((x, y) => x[0] - y[0]);
  const out: Array<[number, number]> = [];
  for (const [a, b] of sorted) {
    const last = out.at(-1);
    if (last && a <= last[1]) last[1] = Math.max(last[1], b);
    else out.push([a, b]);
  }
  return out;
}

/** Parts of `seg` outside `covered` (sorted, disjoint), in traversal order. */
function uncovered(seg: Segment, covered: Array<[number, number]>): Segment[] {
  const out: Segment[] = [];
  let cursor = seg.start;
  for (const [a, b] of [...covered, [seg.end, seg.end] as [number, number]]) {
    if (a > cursor) {
      const part: Segment = { ref: seg.ref, start: cursor, end: Math.min(a, seg.end), strand: seg.strand };
      if (part.start === seg.start && seg.fuzzyLow) part.fuzzyLow = true;
      if (part.end === seg.end && seg.fuzzyHigh) part.fuzzyHigh = true;
      out.push(part);
    }
    cursor = Math.max(cursor, b);
  }
  return seg.strand === 1 ? out : out.reverse();
}

/**
 * Between-position: for each target sequence, both flanking units must map uniquely onto adjacent units of that
 * sequence (spec-core §5.6). Overlapping features on other sequences (e.g. genes in different frames) do not interfere.
 */
function mapBetween(seg: Segment, m: Mapping): Segment[] {
  const k = seg.start;
  if (k === 0) return [];
  const unit = (u: number) =>
    m.overlapping(seg.ref, u, u + 1).map(({ block }) => ({ ref: block.tgtRef, pos: projectInterval(block, u, u + 1)[0] }));
  const left = unit(k - 1);
  const right = unit(k);
  const out: Segment[] = [];
  for (const ref of new Set(left.map((x) => x.ref))) {
    const l = left.filter((x) => x.ref === ref);
    const r = right.filter((x) => x.ref === ref);
    if (l.length !== 1 || r.length !== 1) continue;
    const d = r[0]!.pos - l[0]!.pos;
    if (Math.abs(d) !== 1) continue;
    const boundary = Math.max(l[0]!.pos, r[0]!.pos);
    out.push({ ref, start: boundary, end: boundary, strand: d === 1 ? seg.strand : flip(seg.strand) });
  }
  return out;
}

function groupTargets(outputs: Segment[], kind: Location["kind"], ctx: CoordContext): TargetLocation[] {
  const byRef = new Map<string, Segment[]>();
  for (const s of outputs) {
    const list = byRef.get(s.ref);
    if (list) list.push(s);
    else byRef.set(s.ref, [s]);
  }
  return [...byRef].map(([ref, segs]) => {
    const minus = segs.filter((s) => s.strand === -1).length;
    const orientation = minus === 0 ? "forward" : minus === segs.length ? "reverse" : "mixed";
    if (ctx.unitOf(ref) === "aa" && minus > 0) {
      // Proteins have no strand: report the orientation and keep numeric order N->C.
      segs = (orientation === "reverse" ? [...segs].reverse() : segs).map((s) => ({ ...s, strand: 1 as const }));
    }
    return { location: { outer: ref, kind, segments: segs }, orientation };
  });
}
