// REST API (design §9) on node:http. JSON in, JSON out; CORS open for browser clients (UI, TogoStanza).
import { readFileSync } from "node:fs";
import { createServer, type IncomingMessage, type Server, type ServerResponse } from "node:http";
import {
  decodeLocationId,
  DEFAULT_BASE,
  formatLocationId,
  LocationSemanticError,
  LocationSyntaxError,
  locationIri,
  parseLocationId,
  toFaldo,
  type CodonMode,
  type CoordContext,
  type Location,
  type Segment,
} from "@togocoord/core";
import { CATEGORIES, type Category } from "./category.ts";
import { convert, type Conversion, type Target } from "./search.ts";
import type { StoreSet } from "./stores.ts";

export interface ApiOptions {
  /** Base of location IRIs (default https://togocoord.example.org/). */
  base?: string;
  /** Maximum number of locations in one POST /v1/convert (default 1000). */
  maxBatch?: number;
  /** Maximum total length of an input location, in bases or residues (default 5,000,000). */
  maxInputLength?: number;
  /** Maximum number of results per conversion (default 1000); `truncated` reports a cut. */
  maxResults?: number;
  /** Tags of preferred sequences for equal-cost choices (default: MANE Select, MANE Plus Clinical). */
  prefer?: string[];
}

export const DEFAULT_PREFER = ["MANE Select", "MANE Plus Clinical"];

/** Web UI files (service/public), served at `/` and `/ui/<file>`. */
const UI_FILES: Record<string, { type: string; path: URL }> = {
  "index.html": { type: "text/html; charset=utf-8", path: new URL("../public/index.html", import.meta.url) },
  "app.js": { type: "text/javascript; charset=utf-8", path: new URL("../public/app.js", import.meta.url) },
  "style.css": { type: "text/css; charset=utf-8", path: new URL("../public/style.css", import.meta.url) },
};

class HttpError extends Error {
  readonly status: number;
  readonly extra: Record<string, unknown>;
  constructor(status: number, message: string, extra: Record<string, unknown> = {}) {
    super(message);
    this.status = status;
    this.extra = extra;
  }
}

export function createApi(stores: StoreSet, options: ApiOptions = {}): Server {
  const base = options.base ?? DEFAULT_BASE;
  const maxBatch = options.maxBatch ?? 1000;
  const maxInputLength = options.maxInputLength ?? 5_000_000;
  const maxResults = options.maxResults ?? 1000;
  const prefer = options.prefer ?? DEFAULT_PREFER;
  const ctx = stores.context();

  /** Location IDs; `namespace:accession` alone means the whole sequence (1..length). */
  const parse = (text: string | null | undefined): Location => {
    if (!text) throw new HttpError(400, "missing parameter 'loc'");
    try {
      return parseLocationId(text, ctx, { wholeSequence: true });
    } catch (e) {
      if (e instanceof LocationSyntaxError) throw new HttpError(400, e.message, { position: e.position });
      if (e instanceof LocationSemanticError) throw new HttpError(400, e.message);
      throw e;
    }
  };

  const targets = (values: string[]): Target[] | undefined => {
    if (values.length === 0) return undefined;
    return values.map((v) => {
      if ((CATEGORIES as readonly string[]).includes(v)) return { category: v as Category };
      if (v.includes(":")) return { ref: parse(`${v}:1`).outer };
      if (stores.registry.get(v)) return { namespace: v.toLowerCase() };
      throw new HttpError(400, `unknown target '${v}' (a category ${CATEGORIES.join("/")}, a namespace, or namespace:accession)`);
    });
  };

  /** `taxon`: an NCBI taxon ID (`10090`, `taxon:10090`) or the name of a loaded species (`Mus musculus`, `mouse`). */
  const taxonParam = (v: unknown): number | undefined => {
    if (v === undefined || v === null || v === "") return undefined;
    const text = String(v).trim();
    const id = /^(?:taxon:|ncbitaxon:)?(\d+)$/i.exec(text);
    if (id) return Number(id[1]);
    const name = text.toLowerCase();
    // Scientific name, or the common name in parentheses ("Mus musculus (house mouse)": "house mouse", "mouse").
    const named = (n: string) => {
      const [, scientific, common] = /^(.*?)(?:\s*\((.*)\))?$/.exec(n.toLowerCase())!;
      return n.toLowerCase() === name || scientific === name || (common !== undefined && (common === name || common.endsWith(` ${name}`)));
    };
    const hit = stores.species().find((s) => s.names.some(named));
    if (!hit) throw new HttpError(400, `unknown species '${text}' (an NCBI taxon ID or a loaded organism name)`);
    return hit.taxon;
  };
  const assemblyParam = (v: unknown): string | undefined => {
    if (v === undefined || v === null || v === "") return undefined;
    const name = String(v).trim().toLowerCase();
    const hit = stores.species().flatMap((s) => s.assemblies).find((a) => a.toLowerCase() === name);
    if (!hit) throw new HttpError(400, `unknown assembly '${String(v)}' (loaded: ${stores.species().flatMap((s) => s.assemblies).join(", ")})`);
    return hit;
  };

  interface Scope {
    taxon?: number;
    assembly?: string;
  }

  const convertOne = (text: string, to: string[], maxHops: number | undefined, codon: CodonMode, tags: string[] = [], scope: Scope = {}) => {
    const loc = parse(text);
    const length = loc.segments.reduce((n, s) => n + (s.end - s.start) / (ctx.unitOf(s.ref) === "aa" ? 3 : 1), 0);
    if (length > maxInputLength) {
      throw new HttpError(413, `input location spans ${length} bases/residues; at most ${maxInputLength} are converted per request`);
    }
    const t = targets(to);
    const found = convert(
      stores,
      loc,
      { ...(t && { to: t }), ...(maxHops !== undefined && { maxHops }), prefer, ...(tags.length === 0 && { maxResults: maxResults + 1 }), ...scope },
      ctx,
    );
    // `tag` keeps only targets carrying one of the tags (e.g. tag=MANE Select).
    const results = tags.length ? found.filter((r) => r.tags.some((x) => tags.includes(x))) : found;
    return {
      input: formatLocationId(loc, ctx, codon),
      ...(stores.taxonOf(loc.outer) !== undefined && { inputTaxon: stores.taxonOf(loc.outer) }),
      results: results.slice(0, maxResults).map((r) => conversionJson(r, ctx, base, codon, stores)),
      ...(results.length > maxResults && { truncated: true }),
    };
  };

  const routes: Array<[string, RegExp, (m: RegExpMatchArray, q: URLSearchParams, body: unknown, req: IncomingMessage) => unknown]> = [
    ["GET", /^\/v1\/meta$/, () => ({ stores: stores.meta(), species: stores.species(), base })],
    [
      "GET",
      /^\/v1\/convert$/,
      (_m, q) =>
        convertOne(q.get("loc") ?? "", q.getAll("to"), intParam(q, "maxHops"), codonParam(q), q.getAll("tag"), {
          taxon: taxonParam(q.get("taxon")),
          assembly: assemblyParam(q.get("assembly")),
        }),
    ],
    [
      "POST",
      /^\/v1\/convert$/,
      (_m, q, body) => {
        const b = (body ?? {}) as { locations?: unknown; to?: unknown; maxHops?: unknown; codon?: unknown; tag?: unknown; taxon?: unknown; assembly?: unknown };
        if (!Array.isArray(b.locations) || !b.locations.every((x) => typeof x === "string")) {
          throw new HttpError(400, "body must be {\"locations\": [\"<Location ID>\", ...], \"to\"?: string | string[]}");
        }
        if (b.locations.length > maxBatch) throw new HttpError(413, `at most ${maxBatch} locations per request`);
        const to = b.to === undefined ? [] : Array.isArray(b.to) ? b.to.map(String) : [String(b.to)];
        const codon: CodonMode = b.codon === "never" ? "never" : codonParam(q);
        const hops = typeof b.maxHops === "number" ? b.maxHops : intParam(q, "maxHops");
        const scope = { taxon: taxonParam(b.taxon ?? q.get("taxon")), assembly: assemblyParam(b.assembly ?? q.get("assembly")) };
        return {
          results: (b.locations as string[]).map((text) => {
            try {
              return convertOne(text, to, hops, codon, b.tag === undefined ? [] : [b.tag].flat().map(String), scope);
            } catch (e) {
              if (e instanceof HttpError) return { input: text, error: e.message, ...e.extra };
              throw e;
            }
          }),
        };
      },
    ],
    [
      "GET",
      /^\/v1\/location$/,
      (_m, q) => {
        const loc = parse(q.get("loc"));
        const codon = codonParam(q);
        return {
          id: formatLocationId(loc, ctx, codon),
          iri: locationIri(loc, ctx, base),
          sequence: loc.outer,
          unit: ctx.unitOf(loc.outer),
          kind: loc.kind,
          segments: loc.segments.map((s) => segmentJson(s, ctx)),
        };
      },
    ],
    ["GET", /^\/v1\/location\/faldo$/, (_m, q) => new JsonLd(toFaldo(parse(q.get("loc")), ctx, { base }, codonParam(q)))],
    [
      "GET",
      /^\/v1\/sequences\/([^/]+)$/,
      (m) => {
        const ref = decodeURIComponent(m[1]!);
        const seq = stores.sequence(ref);
        if (!seq) throw new HttpError(404, `unknown sequence '${ref}'`);
        return { ...seq, category: stores.category(ref), identical: stores.identical(ref) };
      },
    ],
    [
      "GET",
      /^\/v1\/sequences\/([^/]+)\/edges$/,
      (m) => {
        const ref = decodeURIComponent(m[1]!);
        return { sequence: ref, edges: stores.edges(ref).map(({ id: _id, ...e }) => e) };
      },
    ],
    [
      "GET",
      /^\/v1\/annotations$/,
      (_m, q) => {
        const loc = parse(q.get("loc"));
        const seen = new Set<string>();
        const annotations = loc.segments.flatMap((s) =>
          stores.annotations(s.ref, s.start === s.end ? s.start - 1 : s.start, s.start === s.end ? s.start + 1 : s.end),
        );
        // The index holds each feature's bounding interval; keep features whose own segments overlap the input
        // (a position in an intron is not "in" the transcript).
        const overlaps = (text: string) => {
          let feature: Location;
          try {
            feature = parseLocationId(text, ctx);
          } catch {
            return true;
          }
          return feature.segments.some((f) =>
            loc.segments.some((s) => f.ref === s.ref && f.start < Math.max(s.end, s.start + 1) && Math.max(s.start, s.end === s.start ? s.start - 1 : s.start) < f.end),
          );
        };
        return {
          input: formatLocationId(loc, ctx),
          annotations: annotations.filter((a) => {
            const key = `${a.type}\t${a.location}`;
            if (seen.has(key)) return false;
            seen.add(key);
            return overlaps(a.location);
          }),
        };
      },
    ],
  ];

  return createServer(async (req, res) => {
    const url = new URL(req.url ?? "/", "http://localhost");
    res.setHeader("Access-Control-Allow-Origin", "*");
    res.setHeader("Access-Control-Allow-Headers", "Content-Type");
    try {
      if (req.method === "OPTIONS") return send(res, 204, undefined);
      if (req.method === "GET" && (url.pathname === "/" || url.pathname.startsWith("/ui/"))) {
        const file = UI_FILES[url.pathname === "/" ? "index.html" : url.pathname.slice(4)];
        if (!file) throw new HttpError(404, `no such UI file ${url.pathname}`);
        res.writeHead(200, { "Content-Type": file.type, "Cache-Control": "no-cache" });
        return void res.end(readFileSync(file.path));
      }
      // identifiers.org-style IRI of a sequence (no location): the sequence itself, not a range on it.
      if (req.method === "GET" && /^\/[A-Za-z][\w.-]*:[^/:]+$/.test(url.pathname)) {
        const text = decodeLocationId(url.pathname.slice(1));
        const accept = req.headers.accept ?? "";
        let ref: string;
        try {
          const i = text.indexOf(":");
          ref = ctx.registry.refKey(text.slice(0, i), text.slice(i + 1));
        } catch (e) {
          throw new HttpError(400, (e as Error).message);
        }
        res.setHeader("Location", accept.includes("text/html") ? `/?loc=${encodeURIComponent(ref)}` : `/v1/sequences/${encodeURIComponent(ref)}`);
        return send(res, 303, undefined);
      }
      // identifiers.org-style resolution: /<namespace>:<accession>:<location>
      if (req.method === "GET" && /^\/[A-Za-z][\w.-]*:[^/]+:/.test(url.pathname)) {
        const loc = parse(decodeLocationId(url.pathname.slice(1)));
        const accept = req.headers.accept ?? "";
        if (accept.includes("text/html")) {
          // Browsers get the web UI for this location.
          res.setHeader("Location", `/?loc=${encodeURIComponent(formatLocationId(loc, ctx))}`);
          return send(res, 303, undefined);
        }
        if (accept.includes("application/json") && !accept.includes("ld+json")) {
          res.setHeader("Location", `/v1/location?loc=${encodeURIComponent(formatLocationId(loc, ctx))}`);
          return send(res, 303, undefined);
        }
        return send(res, 200, new JsonLd(toFaldo(loc, ctx, { base })));
      }
      for (const [method, pattern, handler] of routes) {
        const m = url.pathname.match(pattern);
        if (!m) continue;
        if (req.method !== method) continue;
        const body = method === "POST" ? await readJson(req) : undefined;
        return send(res, 200, await handler(m, url.searchParams, body, req));
      }
      throw new HttpError(routes.some(([, p]) => p.test(url.pathname)) ? 405 : 404, `no route for ${req.method} ${url.pathname}`);
    } catch (e) {
      if (e instanceof HttpError) return send(res, e.status, { error: e.message, ...e.extra });
      console.error(e);
      return send(res, 500, { error: "internal error" });
    }
  });
}

class JsonLd {
  readonly value: unknown;
  constructor(value: unknown) {
    this.value = value;
  }
}

function send(res: ServerResponse, status: number, body: unknown): void {
  if (body === undefined) {
    res.writeHead(status).end();
    return;
  }
  const ld = body instanceof JsonLd;
  res.writeHead(status, { "Content-Type": ld ? "application/ld+json" : "application/json; charset=utf-8" });
  res.end(JSON.stringify(ld ? body.value : body));
}


async function readJson(req: IncomingMessage): Promise<unknown> {
  const chunks: Buffer[] = [];
  let size = 0;
  for await (const chunk of req as AsyncIterable<Buffer>) {
    size += chunk.length;
    if (size > 10 * 1024 * 1024) throw new HttpError(413, "request body too large");
    chunks.push(chunk);
  }
  if (size === 0) return undefined;
  try {
    return JSON.parse(Buffer.concat(chunks).toString("utf8"));
  } catch {
    throw new HttpError(400, "request body is not valid JSON");
  }
}

function intParam(q: URLSearchParams, name: string): number | undefined {
  const v = q.get(name);
  if (v === null) return undefined;
  const n = Number(v);
  if (!Number.isInteger(n) || n < 0 || n > 8) throw new HttpError(400, `${name} must be an integer 0..8`);
  return n;
}

function codonParam(q: URLSearchParams): CodonMode {
  return q.get("codon") === "never" ? "never" : "auto";
}

/** Segment in 1-based closed coordinates of its sequence (residues for proteins, with codon positions). */
function segmentJson(s: Segment, ctx: CoordContext): Record<string, unknown> {
  const aa = ctx.unitOf(s.ref) === "aa";
  const at = (u: number) => (aa ? { residue: Math.floor(u / 3) + 1, codonPosition: (u % 3) + 1 } : { position: u + 1 });
  const out: Record<string, unknown> = { sequence: s.ref, strand: s.strand === 1 ? "+" : "-" };
  if (s.start === s.end) Object.assign(out, { between: [at(s.start - 1), at(s.start)] });
  else Object.assign(out, { begin: at(s.start), end: at(s.end - 1) });
  if (s.fuzzyLow) out.fuzzyLow = true;
  if (s.fuzzyHigh) out.fuzzyHigh = true;
  if (s.uncertain) out.uncertain = true;
  return out;
}

function conversionJson(r: Conversion, ctx: CoordContext, base: string, codon: CodonMode, stores: StoreSet): Record<string, unknown> {
  const organism = r.taxon !== undefined ? stores.organismName(r.taxon) : undefined;
  return {
    location: formatLocationId(r.location, ctx, codon),
    iri: locationIri(r.location, ctx, base),
    sequence: r.location.outer,
    category: r.category,
    ...(r.tags.length && { tags: r.tags }),
    ...(r.taxon !== undefined && { taxon: r.taxon }),
    ...(organism && { organism }),
    cost: r.cost,
    approximate: r.approximate,
    orientation: r.orientation,
    path: r.path.map((s) => ({
      kind: s.kind,
      from: s.from,
      to: s.to,
      direction: s.direction,
      input: s.input,
      unmapped: s.unmapped,
      ...(s.edgeLocation !== undefined && { edgeLocation: s.edgeLocation }),
      attributes: s.attributes,
      validation: s.validation,
      ...(s.provenance && { provenance: s.provenance }),
    })),
  };
}
