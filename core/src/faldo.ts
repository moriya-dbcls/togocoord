// Location -> FALDO JSON-LD and Location ID -> IRI (design §3.5, §3.6).
// Vocabulary: http://biohackathon.org/resource/faldo (checked against faldo.ttl, OBF/FALDO master).
import { formatLocationId, type CodonMode, type CoordContext, type Location, type Segment } from "./location.ts";

export const FALDO = "http://biohackathon.org/resource/faldo#";
export const RDF = "http://www.w3.org/1999/02/22-rdf-syntax-ns#";

export interface FaldoOptions {
  /** Base of location IRIs (identifiers.org style: `<base><namespace>:<accession>:<location>`). */
  base?: string;
  /** IRI of a sequence; default identifiers.org (`https://identifiers.org/refseq:NC_000001.11`). */
  referenceIri?: (ref: string) => string;
  /** Namespace of TogoCoord terms (codonPosition), until FALDO adopts them. */
  vocabulary?: string;
}

export const DEFAULT_BASE = "https://togocoord.dbcls.jp/";
export const DEFAULT_VOCABULARY = "https://togocoord.dbcls.jp/ontology#";

/** Percent-encode the characters of a Location ID that may not appear in an IRI path (`<`, `>`, `^`). */
export function encodeLocationId(id: string): string {
  return id.replace(/</g, "%3C").replace(/>/g, "%3E").replace(/\^/g, "%5E");
}

export function decodeLocationId(path: string): string {
  return decodeURIComponent(path);
}

export function locationIri(loc: Location, ctx: CoordContext, base = DEFAULT_BASE): string {
  return base + encodeLocationId(formatLocationId(loc, ctx));
}

type Node = Record<string, unknown>;

/**
 * FALDO JSON-LD of a location. The top-level node is the location itself (`@id` = its IRI), so annotations can point
 * to it with `faldo:location`.
 *   single residue/base   -> faldo:ExactPosition
 *   range                 -> faldo:Region with faldo:begin / faldo:end (on the reverse strand begin > end)
 *   a^b                   -> faldo:InBetweenPosition with faldo:after / faldo:before
 *   a.b                   -> faldo:InRangePosition with faldo:begin / faldo:end
 *   join(...) / order(...) -> faldo:ListOfRegions (rdf:Seq) / faldo:BagOfRegions (rdf:Bag), members rdf:_1, rdf:_2, ...
 *   < / >                 -> the begin / end is typed faldo:FuzzyPosition (keeping faldo:position)
 *   codon extension       -> tgc:codonPosition 1..3 on the position (aa sequences)
 */
export function toFaldo(loc: Location, ctx: CoordContext, options: FaldoOptions = {}, codon: CodonMode = "auto"): Node {
  const refIri = options.referenceIri ?? ((ref: string) => `https://identifiers.org/${ref}`);
  const context = {
    faldo: FALDO,
    rdf: RDF,
    tgc: options.vocabulary ?? DEFAULT_VOCABULARY,
  };
  const segmentNode = (s: Segment): Node => {
    const aa = ctx.unitOf(s.ref) === "aa";
    const strandType = aa ? [] : [s.strand === 1 ? "faldo:ForwardStrandPosition" : "faldo:ReverseStrandPosition"];
    const position = (unit: number, fuzzy: boolean, side: "first" | "last"): Node => {
      const node: Node = {
        "@type": [fuzzy ? "faldo:FuzzyPosition" : "faldo:ExactPosition", ...strandType],
        "faldo:position": aa ? Math.floor(unit / 3) + 1 : unit + 1,
        "faldo:reference": { "@id": refIri(s.ref) },
      };
      const k = (unit % 3) + 1;
      if (aa && codon !== "never" && k !== (side === "first" ? 1 : 3)) node["tgc:codonPosition"] = k;
      return node;
    };

    if (s.start === s.end) {
      // Between units s.start - 1 and s.start.
      const [left, right] = s.strand === 1 ? [s.start - 1, s.start] : [s.start, s.start - 1];
      const exact = (u: number) => {
        const n = position(u, false, "first");
        if (aa && codon !== "never") n["tgc:codonPosition"] = (u % 3) + 1;
        if (aa && s.start % 3 === 0) delete n["tgc:codonPosition"];
        return n;
      };
      return { "@type": "faldo:InBetweenPosition", "faldo:after": exact(left), "faldo:before": exact(right) };
    }

    const last = s.end - 1;
    if (s.uncertain) {
      return {
        "@type": "faldo:InRangePosition",
        "faldo:begin": position(s.start, false, "first"),
        "faldo:end": position(last, false, "last"),
      };
    }
    const singleNt = !aa && s.end - s.start === 1;
    const singleResidue = aa && s.start % 3 === 0 && s.end - s.start === 3;
    const singleUnit = aa && s.end - s.start === 1 && codon !== "never";
    if ((singleNt || singleResidue || singleUnit) && !s.fuzzyLow && !s.fuzzyHigh) {
      const n = position(s.start, false, "first");
      if (singleUnit) n["tgc:codonPosition"] = (s.start % 3) + 1;
      return n;
    }
    // Biological begin/end: on the reverse strand the begin is the numerically higher end.
    const low = position(s.start, !!s.fuzzyLow, "first");
    const high = position(last, !!s.fuzzyHigh, "last");
    return { "@type": "faldo:Region", ...(s.strand === 1 ? { "faldo:begin": low, "faldo:end": high } : { "faldo:begin": high, "faldo:end": low }) };
  };

  const id = locationIri(loc, ctx, options.base);
  if (loc.segments.length === 1) return { "@context": context, "@id": id, ...segmentNode(loc.segments[0]!) };
  const collection: Node = { "@context": context, "@id": id, "@type": loc.kind === "order" ? "faldo:BagOfRegions" : "faldo:ListOfRegions" };
  loc.segments.forEach((s, i) => {
    collection[`rdf:_${i + 1}`] = segmentNode(s);
  });
  return collection;
}
