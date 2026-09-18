// Location IDs <-> semantic segment lists (spec-core §2, §3.3, §3.5).
import { LocationSemanticError, LocationSyntaxError } from "./errors.ts";
import { NamespaceRegistry, splitRef, type Unit } from "./namespace.ts";
import { parseLocationText, type AstNode, type AstPos } from "./parser.ts";

/** Oriented interval in internal units: 0-based half-open; `start === end` is a between-position. */
export interface Segment {
  ref: string;
  start: number;
  end: number;
  strand: 1 | -1;
  /** '<' on the numerically lower end. */
  fuzzyLow?: boolean;
  /** '>' on the numerically higher end. */
  fuzzyHigh?: boolean;
  /** INSDC one-of (`a.b`). */
  uncertain?: boolean;
}

/** Segments in traversal order (5'->3' / N->C of the described feature). */
export interface Location {
  outer: string;
  kind: "join" | "order";
  segments: Segment[];
}

export type CodonMode = "auto" | "never";

export interface CoordContext {
  registry: NamespaceRegistry;
  unitOf(ref: string): Unit;
  lengthOf(ref: string): number | undefined;
}

export interface ContextOptions {
  registry?: NamespaceRegistry;
  /** Per-sequence unit overrides (e.g. for refget keys or test sequences). */
  units?: Record<string, Unit> | ((ref: string) => Unit | undefined);
  /** Sequence lengths in residues (not units); enables range checks. */
  lengths?: Record<string, number> | ((ref: string) => number | undefined);
}

export function createContext(options: ContextOptions = {}): CoordContext {
  const registry = options.registry ?? new NamespaceRegistry();
  const lookup = <T>(src: Record<string, T> | ((ref: string) => T | undefined) | undefined, ref: string) =>
    typeof src === "function" ? src(ref) : src?.[ref];
  return {
    registry,
    unitOf(ref) {
      const unit = lookup(options.units, ref) ?? registry.defaultUnit(ref) ?? (ref.includes("#") ? "nt" : undefined);
      if (!unit) throw new LocationSemanticError(`unit of '${ref}' is unknown; supply it in the context`);
      return unit;
    },
    lengthOf(ref) {
      return lookup(options.lengths, ref);
    },
  };
}

/** Split `namespace:accession:location` at its first two colons. */
export function splitLocationId(text: string): { namespace: string; accession: string; locationText: string } {
  const s = text.replace(/\s+/g, "");
  const i = s.indexOf(":");
  const j = i < 0 ? -1 : s.indexOf(":", i + 1);
  if (i <= 0 || j <= i + 1 || j === s.length - 1) {
    throw new LocationSyntaxError("expected namespace:accession:location", Math.max(i, 0));
  }
  return { namespace: s.slice(0, i), accession: s.slice(i + 1, j), locationText: s.slice(j + 1) };
}

export interface ParseOptions {
  /**
   * Accept `namespace:accession` without a location as the whole sequence (1..length), using ctx.lengthOf.
   * Input shorthand only: formatting always writes the explicit range (spec-core §3.1).
   */
  wholeSequence?: boolean;
}

export function parseLocationId(text: string, ctx: CoordContext, options: ParseOptions = {}): Location {
  if (options.wholeSequence) {
    const s = text.replace(/\s+/g, "");
    const i = s.indexOf(":");
    if (i > 0 && s.indexOf(":", i + 1) < 0) return wholeSequence(s.slice(0, i), s.slice(i + 1), ctx);
  }
  const { namespace, accession, locationText } = splitLocationId(text);
  const outer = ctx.registry.refKey(namespace, accession);
  const ast = parseLocationText(locationText);
  const state = { order: false };
  const segments = build(ast, outer, splitRef(outer).namespace, ctx, true, state);
  for (const seg of segments) {
    if (seg.strand === -1 && ctx.unitOf(seg.ref) === "aa") {
      throw new LocationSemanticError(`complement is not allowed on protein sequence '${seg.ref}'`);
    }
  }
  return { outer, kind: state.order ? "order" : "join", segments };
}

function wholeSequence(namespace: string, accession: string, ctx: CoordContext): Location {
  const outer = ctx.registry.refKey(namespace, accession);
  const length = ctx.lengthOf(outer);
  if (!length || length < 1) {
    throw new LocationSemanticError(`the length of '${outer}' is unknown; give an explicit location (namespace:accession:location)`);
  }
  const units = ctx.unitOf(outer) === "aa" ? 3 : 1;
  return { outer, kind: "join", segments: [{ ref: outer, start: 0, end: length * units, strand: 1 }] };
}

function build(
  node: AstNode,
  ref: string,
  namespace: string,
  ctx: CoordContext,
  orderAllowed: boolean,
  state: { order: boolean },
): Segment[] {
  const here = node.ref !== undefined ? ctx.registry.refKey(namespace, node.ref) : ref;
  switch (node.type) {
    case "complement":
      return build(node.child, here, namespace, ctx, orderAllowed, state)
        .reverse()
        .map((s) => ({ ...s, strand: s.strand === 1 ? -1 : 1 }));
    case "order":
      if (!orderAllowed) throw new LocationSemanticError("order() is only allowed at the top level");
      state.order = true;
      return node.children.flatMap((c) => build(c, here, namespace, ctx, false, state));
    case "join":
      return node.children.flatMap((c) => build(c, here, namespace, ctx, false, state));
    default:
      return [spanToSegment(node, here, ctx)];
  }
}

function spanToSegment(node: Exclude<AstNode, { children: unknown } | { child: unknown }>, ref: string, ctx: CoordContext): Segment {
  const aa = ctx.unitOf(ref) === "aa";
  const first = (p: AstPos) => (aa ? 3 * (p.value - 1) + (p.codon ?? 1) - 1 : p.value - 1);
  const last = (p: AstPos) => (aa ? 3 * (p.value - 1) + (p.codon ?? 3) : p.value);
  const positions =
    node.type === "point" ? [node.pos] : node.type === "range" ? [node.begin, node.end] : node.type === "between" ? [node.left, node.right] : [];
  if (!aa && positions.some((p) => p.codon !== undefined)) {
    throw new LocationSemanticError(`codon extension is only allowed on protein sequences ('${ref}')`);
  }

  let seg: Segment;
  switch (node.type) {
    case "point":
      seg = { ref, start: first(node.pos), end: last(node.pos), strand: 1 };
      break;
    case "range":
      seg = { ref, start: first(node.begin), end: last(node.end), strand: 1 };
      if (node.fuzzyLow) seg.fuzzyLow = true;
      if (node.fuzzyHigh) seg.fuzzyHigh = true;
      break;
    case "oneof":
      seg = { ref, start: first({ value: node.begin }), end: last({ value: node.end }), strand: 1, uncertain: true };
      break;
    case "between": {
      // Right unit of the boundary; the left unit must be immediately before it.
      const leftUnit = aa ? 3 * (node.left.value - 1) + (node.left.codon ?? 3) - 1 : node.left.value - 1;
      const rightUnit = aa ? 3 * (node.right.value - 1) + (node.right.codon ?? 1) - 1 : node.right.value - 1;
      if (rightUnit !== leftUnit + 1) {
        throw new LocationSemanticError(
          "'^' needs adjacent positions (the circular n^1 form is not supported in v0.1)",
        );
      }
      seg = { ref, start: rightUnit, end: rightUnit, strand: 1 };
      break;
    }
  }
  const length = ctx.lengthOf(ref);
  if (length !== undefined && seg.end > length * (aa ? 3 : 1)) {
    throw new LocationSemanticError(`location exceeds the length of '${ref}' (${length})`);
  }
  return seg;
}

/** Canonical Location ID text (spec-core §3.5). */
export function formatLocationId(loc: Location, ctx: CoordContext, codon: CodonMode = "auto"): string {
  return `${loc.outer}:${formatLocation(loc, ctx, codon)}`;
}

export function formatLocation(loc: Location, ctx: CoordContext, codon: CodonMode = "auto"): string {
  const segs = loc.segments;
  if (segs.length === 0) throw new LocationSemanticError("cannot format an empty location");
  const outerNs = splitRef(loc.outer).namespace;
  const prefix = (ref: string) => {
    if (ref === loc.outer) return "";
    const { namespace, accession } = splitRef(ref);
    if (namespace !== outerNs) {
      throw new LocationSemanticError(`remote reference '${ref}' is outside namespace '${outerNs}'`);
    }
    return `${accession}:`;
  };
  const span = (s: Segment) => prefix(s.ref) + formatSpan(s, ctx.unitOf(s.ref), codon);
  const oriented = (s: Segment) => (s.strand === -1 ? `complement(${span(s)})` : span(s));

  if (segs.length === 1) return oriented(segs[0]!);
  const first = segs[0]!;
  if (segs.every((s) => s.strand === -1 && s.ref === first.ref)) {
    return `complement(${loc.kind}(${[...segs].reverse().map(span).join(",")}))`;
  }
  return `${loc.kind}(${segs.map(oriented).join(",")})`;
}

function formatSpan(s: Segment, unit: Unit, codon: CodonMode): string {
  const lo = s.fuzzyLow ? "<" : "";
  const hi = s.fuzzyHigh ? ">" : "";
  if (unit === "nt") {
    if (s.start === s.end) return `${s.start}^${s.start + 1}`;
    if (s.uncertain) return `${s.start + 1}.${s.end}`;
    if (s.end - s.start === 1 && !lo && !hi) return `${s.end}`;
    return `${lo}${s.start + 1}..${hi}${s.end}`;
  }

  const residue = (u: number) => Math.floor(u / 3) + 1;
  const phase = (u: number) => (u % 3) + 1;
  if (s.start === s.end) {
    const k = s.start;
    if (k % 3 === 0) return `${k / 3}^${k / 3 + 1}`;
    if (codon === "never") return `${residue(k)}`;
    return `${residue(k - 1)}c${phase(k - 1)}^${residue(k)}c${phase(k)}`;
  }
  const bRes = residue(s.start);
  const eRes = residue(s.end - 1);
  if (s.uncertain) return `${bRes}.${eRes}`;
  const bC = codon === "never" ? 1 : phase(s.start);
  const eC = codon === "never" ? 3 : phase(s.end - 1);
  if (!lo && !hi && bRes === eRes) {
    if (bC === 1 && eC === 3) return `${bRes}`;
    if (bC === eC) return `${bRes}c${bC}`;
  }
  return `${lo}${bRes}${bC === 1 ? "" : `c${bC}`}..${hi}${eRes}${eC === 3 ? "" : `c${eC}`}`;
}

/** Parse and re-emit in canonical form. */
export function canonicalize(text: string, ctx: CoordContext, codon: CodonMode = "auto"): string {
  return formatLocationId(parseLocationId(text, ctx), ctx, codon);
}
