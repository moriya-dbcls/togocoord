// REST API (design §9) on node:http. JSON in, JSON out; CORS open for browser clients (UI, TogoStanza).
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
}

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
  const ctx = stores.context();

  const parse = (text: string | null | undefined): Location => {
    if (!text) throw new HttpError(400, "missing parameter 'loc'");
    try {
      return parseLocationId(text, ctx);
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

  const convertOne = (text: string, to: string[], maxHops: number | undefined, codon: CodonMode) => {
    const loc = parse(text);
    const t = targets(to);
    const results = convert(stores, loc, { ...(t && { to: t }), ...(maxHops !== undefined && { maxHops }) }, ctx);
    return { input: formatLocationId(loc, ctx, codon), results: results.map((r) => conversionJson(r, ctx, base, codon)) };
  };

  const routes: Array<[string, RegExp, (m: RegExpMatchArray, q: URLSearchParams, body: unknown, req: IncomingMessage) => unknown]> = [
    ["GET", /^\/v1\/meta$/, () => ({ stores: stores.meta(), base })],
    [
      "GET",
      /^\/v1\/convert$/,
      (_m, q) => convertOne(q.get("loc") ?? "", q.getAll("to"), intParam(q, "maxHops"), codonParam(q)),
    ],
    [
      "POST",
      /^\/v1\/convert$/,
      (_m, q, body) => {
        const b = (body ?? {}) as { locations?: unknown; to?: unknown; maxHops?: unknown; codon?: unknown };
        if (!Array.isArray(b.locations) || !b.locations.every((x) => typeof x === "string")) {
          throw new HttpError(400, "body must be {\"locations\": [\"<Location ID>\", ...], \"to\"?: string | string[]}");
        }
        if (b.locations.length > maxBatch) throw new HttpError(413, `at most ${maxBatch} locations per request`);
        const to = b.to === undefined ? [] : Array.isArray(b.to) ? b.to.map(String) : [String(b.to)];
        const codon: CodonMode = b.codon === "never" ? "never" : codonParam(q);
        const hops = typeof b.maxHops === "number" ? b.maxHops : intParam(q, "maxHops");
        return {
          results: (b.locations as string[]).map((text) => {
            try {
              return convertOne(text, to, hops, codon);
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
        return {
          input: formatLocationId(loc, ctx),
          annotations: annotations.filter((a) => {
            const key = `${a.type}\t${a.location}`;
            if (seen.has(key)) return false;
            seen.add(key);
            return true;
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
      // identifiers.org-style resolution: /<namespace>:<accession>:<location>
      if (req.method === "GET" && /^\/[A-Za-z][\w.-]*:[^/]+:/.test(url.pathname)) {
        const loc = parse(decodeLocationId(url.pathname.slice(1)));
        const accept = req.headers.accept ?? "";
        if (accept.includes("text/html")) return sendHtml(res, locationPage(loc, ctx, base));
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

function sendHtml(res: ServerResponse, html: string): void {
  res.writeHead(200, { "Content-Type": "text/html; charset=utf-8" }).end(html);
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

function conversionJson(r: Conversion, ctx: CoordContext, base: string, codon: CodonMode): Record<string, unknown> {
  return {
    location: formatLocationId(r.location, ctx, codon),
    iri: locationIri(r.location, ctx, base),
    sequence: r.location.outer,
    category: r.category,
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

function escapeHtml(s: string): string {
  return s.replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]!);
}

function locationPage(loc: Location, ctx: CoordContext, base: string): string {
  const id = formatLocationId(loc, ctx);
  const q = encodeURIComponent(id);
  return `<!doctype html><meta charset="utf-8"><title>${escapeHtml(id)}</title>
<h1>${escapeHtml(id)}</h1>
<ul>
<li><a href="/v1/location?loc=${q}">location (JSON)</a></li>
<li><a href="/v1/location/faldo?loc=${q}">FALDO JSON-LD</a></li>
<li><a href="/v1/convert?loc=${q}">directly connected sequences</a></li>
<li>convert to: ${["genome", "transcript", "protein"].map((c) => `<a href="/v1/convert?loc=${q}&amp;to=${c}">${c}</a>`).join(" · ")}</li>
</ul>
<pre>${escapeHtml(JSON.stringify(toFaldo(loc, ctx, { base }), null, 2))}</pre>`;
}
