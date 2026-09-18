// TogoCoord web UI: a thin client of the REST API (spec-service §6). State lives in the URL (?loc=&to=&codon=).
const $ = (sel, root = document) => root.querySelector(sel);

// Without `to`, an example keeps the current selectors; with it, it sets them all (unset ones to their default).
const EXAMPLES = [
  { label: "UniProt residue", loc: "uniprot:P07203:49" },
  { label: "genome codon", loc: "refseq:NC_000003.12:complement(49358132..49358134)" },
  { label: "structure residue", loc: "pdb:2F8A.A:59" },
  { label: "Ensembl protein range", loc: "ensembl:ENSP00000407375.1:40..60" },
  { label: "RefSeq protein", loc: "refseq:NP_036366.3:20" },
  { label: "whole protein", loc: "refseq:NP_000572.2" },
  { label: "overlapping genes (mtDNA)", loc: "refseq:NC_012920.1:8527..8529" },
  { label: "hg19 → GRCh38 (BRAF V600E)", loc: "hg19:chr7:140453136", to: "genome", assembly: "GRCh38", needs: "GRCh37" },
  { label: "hg19 → protein", loc: "hg19:chr7:complement(140453135..140453137)", to: "protein", db: "refseq", needs: "GRCh37" },
  { label: "protein → hg19", loc: "uniprot:P15056:600", to: "genome", assembly: "GRCh37", needs: "GRCh37" },
  { label: "human → mouse UniProt", loc: "uniprot:P07203:49", to: "protein", taxon: "10090", db: "uniprot" },
  { label: "mouse → human genome", loc: "refseq:NC_000075.7:106312500..106312550", to: "genome", taxon: "9606" },
];

const form = $("#query");
const locInput = $("#loc");
const toSelect = $("#to");
const dbSelect = $("#db");
const taxonSelect = $("#taxon");
const assemblySelect = $("#assembly");
const codonBox = $("#codon");
const maneBox = $("#mane");

function el(tag, props = {}, ...children) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === "class") node.className = v;
    else if (k === "dataset") Object.assign(node.dataset, v);
    else if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) node.append(c.nodeType ? c : String(c));
  return node;
}

/**
 * Root of the service, from where this script is served (`<root>/ui/app.js`): the UI also works when a reverse proxy
 * mounts the service under a subdirectory (e.g. https://example.org/togocoord/).
 */
const ROOT = new URL("..", import.meta.url);
const endpoint = (path) => new URL(path.replace(/^\//, ""), ROOT);

async function api(path) {
  const res = await fetch(endpoint(path), { headers: { accept: "application/json" } });
  const body = await res.json();
  if (!res.ok) throw Object.assign(new Error(body.error ?? res.statusText), { body });
  return body;
}

const qs = (params) => new URLSearchParams(Object.entries(params).filter(([, v]) => v !== "" && v !== undefined)).toString();

/** A Location ID that continues from itself: shows what it connects to directly (the target is reset). */
function locationLink(id) {
  return el("a", { class: "location", href: `?${qs({ loc: id })}`, title: "continue from here", onclick: (e) => { e.preventDefault(); run(id, ""); } }, el("code", {}, id));
}

/** Features describing a whole molecule are not useful as annotations of a position. */
const WHOLE_MOLECULE_TYPES = new Set(["region", "chromosome", "scaffold", "source", "supercontig", "contig"]);

function showError(err, text) {
  const box = $("#error");
  let msg = err.message;
  if (typeof err.body?.position === "number" && text) {
    // Positions count characters of the location part after `namespace:accession:`.
    const offset = text.replace(/\s+/g, "").split(":", 2).join(":").length + 1 + err.body.position;
    msg += `\n${text.replace(/\s+/g, "")}\n${" ".repeat(offset)}^`;
  }
  box.textContent = msg;
  box.hidden = false;
}

function segmentRows(segments) {
  // Codon positions are shown only where a range does not start (c1) or end (c3) on a residue boundary.
  const pos = (p, boundary) => (p.residue !== undefined ? `${p.residue}${p.codonPosition !== boundary ? ` (codon position ${p.codonPosition})` : ""}` : p.position);
  return [
    el("tr", {}, el("th", {}, "sequence"), el("th", {}, "strand"), el("th", {}, "begin"), el("th", {}, "end"), el("th", {}, "")),
    ...segments.map((s) =>
      el("tr", {},
        el("td", {}, el("code", {}, s.sequence)),
        el("td", {}, s.strand),
        el("td", {}, s.between ? `after ${pos(s.between[0], 3)}` : `${s.fuzzyLow ? "<" : ""}${pos(s.begin, 1)}`),
        el("td", {}, s.between ? `before ${pos(s.between[1], 1)}` : `${s.fuzzyHigh ? ">" : ""}${pos(s.end, 3)}`),
        el("td", { class: "muted" }, s.uncertain ? "one of" : ""),
      ),
    ),
  ];
}

/** Toggle an extra panel (FALDO or annotations) below a card. */
async function toggleExtra(card, kind, id) {
  const box = $(".extra", card);
  if (!box.hidden && box.dataset.kind === kind) {
    box.hidden = true;
    return;
  }
  box.dataset.kind = kind;
  box.replaceChildren(el("span", { class: "muted" }, "loading…"));
  box.hidden = false;
  try {
    if (kind === "faldo") {
      const res = await fetch(endpoint(`/v1/location/faldo?${qs({ loc: id })}`));
      box.replaceChildren(el("pre", {}, JSON.stringify(await res.json(), null, 2)));
    } else {
      const all = (await api(`/v1/annotations?${qs({ loc: id })}`)).annotations;
      // Narrowest first: the features that actually describe this position come before gene-long ones.
      const span = (a) => { const m = a.location.match(/\d+/g); return m ? Math.max(...m.map(Number)) - Math.min(...m.map(Number)) : 0; };
      const annotations = all.filter((a) => !WHOLE_MOLECULE_TYPES.has(a.type)).sort((a, b) => span(a) - span(b));
      if (!annotations.length) {
        box.replaceChildren(el("span", { class: "muted" }, "no annotations overlap this location"));
        return;
      }
      const label = (a) => {
        const at = a.attributes ?? {};
        const pick = (k) => (at[k] ? [].concat(at[k])[0] : undefined);
        return [pick("gene"), pick("product"), pick("Name"), pick("note")].filter(Boolean).slice(0, 2).join(" · ");
      };
      box.replaceChildren(
        el("ul", { class: "annotations" },
          annotations.slice(0, 100).map((a) => el("li", {}, el("span", { class: "badge" }, a.type), " ", locationLink(a.location), " ", el("span", { class: "muted" }, label(a)))),
          annotations.length > 100 ? el("li", { class: "muted" }, `… ${annotations.length - 100} more`) : null,
        ),
      );
    }
  } catch (err) {
    box.replaceChildren(el("span", { class: "unmapped" }, err.message));
  }
}

function renderPath(result) {
  const nodes = [];
  for (const [i, s] of result.path.entries()) {
    const from = s.direction === "forward" ? s.from : s.to;
    const to = s.direction === "forward" ? s.to : s.from;
    if (i === 0) nodes.push(el("span", { class: "node" }, from));
    const bad = s.validation?.status === "mismatch" || (s.attributes?.exception && s.validation?.basis !== "full");
    const title = [
      `${s.kind}${s.direction === "inverse" ? " (inverse)" : ""}`,
      s.edgeLocation && `edge: ${s.edgeLocation}`,
      s.validation && `validation: ${s.validation.status}${s.validation.detail ? ` — ${s.validation.detail}` : ""}`,
      s.attributes?.exception && `exception: ${s.attributes.exception}`,
      s.provenance?.file && `source: ${s.provenance.file}`,
    ].filter(Boolean).join("\n");
    nodes.push(el("span", { class: `step ${s.kind}${bad ? " bad" : ""}`, title }, s.kind === "identity" ? "identical" : s.kind));
    nodes.push(el("span", { class: "node" }, to));
  }
  return el("div", { class: "path" }, nodes);
}

/** Taxon of the input sequence, to mark results from other species. */
let inputTaxon;

function renderResult(r) {
  const card = el("li", { class: "card" });
  const unmapped = r.path.filter((s) => s.unmapped).map((s) => `${s.unmapped} (at ${s.kind} step)`);
  card.append(
    ...[el("div", { class: "id-row" },
      locationLink(r.location),
      el("span", { class: "badge" }, r.category),
      ...(r.tags ?? []).map((t) => el("span", { class: "badge tag", title: "curated tag of the target sequence" }, t)),
      r.taxon !== undefined && r.taxon !== inputTaxon
        ? el("span", { class: "badge species", title: `taxon ${r.taxon}` }, r.organism ?? `taxon ${r.taxon}`)
        : null,
      r.assembly && multiAssembly.has(r.taxon) ? el("span", { class: "badge species", title: "genome assembly" }, r.assembly) : null,
      el("span", { class: "badge", title: "path cost" }, `cost ${r.cost}`),
      r.approximate ? el("span", { class: "badge warn", title: "the path uses an edge not verified against the sequences; positions may be shifted" }, "approximate") : null,
      r.orientation !== "forward" ? el("span", { class: "badge" }, r.orientation) : null,
      el("button", { type: "button", class: "small", onclick: () => toggleExtra(card, "faldo", r.location) }, "FALDO"),
      el("button", { type: "button", class: "small", onclick: () => toggleExtra(card, "annotations", r.location) }, "Annotations"),
    ),
    renderPath(r),
    unmapped.length ? el("div", { class: "unmapped" }, `not mapped: ${unmapped.join("; ")}`) : null,
    el("div", { class: "extra", hidden: "" })].filter(Boolean),
  );
  return card;
}

/** Species with several assemblies (their genome results carry an assembly badge). */
const multiAssembly = new Set();
/** Loaded assemblies by name, UCSC name and name without patch level (for examples and URLs: hg19, GRCh37). */
const assemblyNames = new Map();

/** Loaded species (and assemblies, shown only when a species has several) for the scope selectors. */
const speciesReady = api("/v1/meta").then(({ species = [], assemblies = [] }) => {
  taxonSelect.append(...species.map((s) => el("option", { value: String(s.taxon) }, s.organism ?? `taxon ${s.taxon}`)));
  for (const s of species) if (s.assemblies.length > 1) multiAssembly.add(s.taxon);
  for (const a of assemblies) for (const n of [a.name, a.name.replace(/\.p\d+$/, ""), a.ucsc]) if (n) assemblyNames.set(n.toLowerCase(), a.name);
  // Only species with a choice; the UCSC name helps those who know hg19 / hg38.
  assemblySelect.append(
    ...assemblies
      .filter((a) => multiAssembly.has(a.taxon))
      .map((a) => {
        const s = species.find((x) => x.taxon === a.taxon);
        const def = s?.defaultAssembly === a.name ? " · annotated" : "";
        const title = `${s?.organism ?? `taxon ${a.taxon}`}${a.accession ? ` · ${a.accession}` : ""}${def}`;
        return el("option", { value: a.name, title }, `${a.name}${a.ucsc ? ` / ${a.ucsc}` : ""}`);
      }),
  );
  assemblySelect.hidden = multiAssembly.size === 0;
  // Examples that need data not loaded here (e.g. GRCh37) are left out.
  for (const b of $("#examples").querySelectorAll("button[data-needs]")) b.hidden = !assemblyNames.has(b.dataset.needs.toLowerCase());
}).catch(() => {});

async function run(loc, to, push = true, taxon = "", assembly = "", db = "") {
  loc = loc.trim();
  locInput.value = loc;
  if (![...toSelect.options].some((o) => o.value === (to ?? ""))) {
    toSelect.append(el("option", { value: to }, to)); // a namespace or sequence target, e.g. from an example or a URL
  }
  toSelect.value = to ?? "";
  taxonSelect.value = taxon;
  assembly = assemblyNames.get((assembly ?? "").toLowerCase()) ?? assembly;
  assemblySelect.value = assembly;
  dbSelect.value = db;
  $("#error").hidden = true;
  if (!loc) return;
  const codon = codonBox.checked ? "" : "never";
  const tag = maneBox.checked ? "MANE Select" : "";
  if (push) history.pushState(null, "", `?${qs({ loc, to, db, taxon, assembly, codon, tag })}`);
  try {
    const [info, conv] = await Promise.all([
      api(`/v1/location?${qs({ loc, codon })}`),
      api(`/v1/convert?${qs({ loc, to, db, taxon, assembly, codon, tag })}`),
    ]);
    $("#input-id").textContent = info.id;
    $("#input-kind").textContent = `${info.unit === "aa" ? "protein" : "nucleotide"}${info.kind === "order" ? " · order" : ""}`;
    // Species and assembly of the input; the name as written (hg19:chr7:...) when it was given that way.
    const scope = [info.organism, info.assembly, info.written && `written as ${info.written.name}`].filter(Boolean).join(" · ");
    $("#input-scope").textContent = scope;
    $("#input-scope").hidden = !scope;
    $("#segments").replaceChildren(...segmentRows(info.segments));
    const inputCard = $("#input .card");
    $(".extra", inputCard).hidden = true;
    $("#input").hidden = false;
    inputTaxon = conv.inputTaxon;
    $("#results").replaceChildren(...conv.results.map(renderResult));
    $("#count").textContent = conv.results.length
      ? `(${conv.results.length}${conv.truncated ? ", truncated — narrow the target or the input" : ""})`
      : tag || db ? "— none matching the filters" : taxon || assembly ? "— none reachable in the selected species / assembly" : "— none reachable";
    $("#output").hidden = false;
  } catch (err) {
    $("#input").hidden = true;
    $("#output").hidden = true;
    showError(err, loc);
  }
}

form.addEventListener("submit", (e) => {
  e.preventDefault();
  run(locInput.value, toSelect.value, true, taxonSelect.value, assemblySelect.value, dbSelect.value);
});
const rerun = () => run(locInput.value, toSelect.value, true, taxonSelect.value, assemblySelect.value, dbSelect.value);
for (const control of [toSelect, dbSelect, taxonSelect, assemblySelect, codonBox, maneBox]) {
  control.addEventListener("change", () => locInput.value && rerun());
}
$("#input .card").addEventListener("click", (e) => {
  const action = e.target.dataset?.action;
  if (action) toggleExtra($("#input .card"), action, $("#input-id").textContent);
});
$("#examples").append(
  "Examples: ",
  ...EXAMPLES.map((x) => {
    const b = el("button", {
      type: "button",
      class: "small",
      title: x.loc,
      onclick: () =>
        x.to
          ? run(x.loc, x.to, true, x.taxon ?? "", x.assembly ?? "", x.db ?? "")
          : run(x.loc, toSelect.value, true, taxonSelect.value, assemblySelect.value, dbSelect.value),
    }, x.label);
    if (x.needs) {
      b.dataset.needs = x.needs;
      b.hidden = true;
    }
    return b;
  }),
);

// ---- Loaded data view ------------------------------------------------------------------------------------------

const fmt = (n) => Number(n).toLocaleString("en-US");
// Sequences without a record of their own (only the end of an edge, e.g. GFF3 transcripts) have no molecule type.
const counts = (o) =>
  Object.entries(o ?? {}).filter(([, n]) => n > 0).map(([k, n]) => `${fmt(n)} ${k === "unknown" ? "referenced only" : k}`).join(", ") || "—";

/** Organism a store belongs to: its recorded organism, else the organism most of its sequences carry. */
function storeOrganism(s) {
  if (s.organism) return { name: s.organism, taxon: s.taxon };
  const t = s.summary?.taxa?.[0];
  return t ? { name: t.organism ?? `taxon ${t.taxon}`, taxon: t.taxon } : { name: "Cross-species / structures", taxon: undefined };
}

async function showData() {
  $("#convert-view").hidden = true;
  $("#data").hidden = false;
  document.querySelectorAll("nav a").forEach((a) => a.classList.toggle("current", a.dataset.view === "data"));
  const { stores } = await api("/v1/meta");
  const groups = new Map();
  for (const s of stores) {
    const o = storeOrganism(s);
    const key = o.taxon !== undefined ? String(o.taxon) : o.name; // "Homo sapiens" and "Homo sapiens (human)" are one group
    if (!groups.has(key)) groups.set(key, { ...o, stores: [] });
    groups.get(key).stores.push(s);
  }
  $("#data-count").textContent = `(${stores.length} stores, ${groups.size} groups)`;
  $("#species").replaceChildren(
    ...[...groups.values()].map((g) =>
      el("div", { class: "species" },
        el("h3", {}, g.name, g.taxon ? el("span", { class: "muted" }, ` · taxon ${g.taxon}`) : null),
        ...g.stores.map((s) =>
          el("div", { class: "card store" },
            el("div", { class: "label" }, s.label ?? s.file),
            el("dl", {},
              s.assembly ? [el("dt", {}, "assembly"), el("dd", {}, `${s.assembly}${s.accession ? ` (${s.accession})` : ""}`)] : null,
              el("dt", {}, "inputs"), el("dd", {}, el("code", {}, (s.inputs ?? "").split(",").join(", "))),
              el("dt", {}, "sequences"), el("dd", {}, counts(s.summary?.sequences)),
              el("dt", {}, "edges"), el("dd", {}, `${counts(s.summary?.edges)}${s.summary?.blocks ? ` (${fmt(s.summary.blocks)} blocks)` : ""}`),
              s.summary?.annotations ? [el("dt", {}, "annotations"), el("dd", {}, fmt(s.summary.annotations))] : null,
              el("dt", {}, "built"), el("dd", {}, `${(s.created ?? "").slice(0, 10)} · ${s.file} · schema ${s.schema}`),
              s.summary?.examples?.length
                ? [el("dt", {}, "try"), el("dd", {}, ...s.summary.examples.map((id) =>
                    el("button", { type: "button", class: "small", title: "convert this location", onclick: () => { showConvert(); run(id, ""); } }, id)))]
                : null,
            ),
          ),
        ),
      ),
    ),
  );
}

function showConvert() {
  $("#data").hidden = true;
  $("#convert-view").hidden = false;
  document.querySelectorAll("nav a").forEach((a) => a.classList.toggle("current", a.dataset.view === "convert"));
}

document.querySelectorAll("nav a").forEach((a) =>
  a.addEventListener("click", (e) => {
    e.preventDefault();
    if (a.dataset.view === "data") {
      history.pushState(null, "", "?view=data");
      showData().catch((err) => showError(err));
    } else {
      history.pushState(null, "", "?");
      showConvert();
    }
  }),
);

async function fromUrl(push = false) {
  const p = new URLSearchParams(location.search);
  if (p.get("view") === "data") {
    showData().catch((err) => showError(err));
    return;
  }
  showConvert();
  codonBox.checked = p.get("codon") !== "never";
  maneBox.checked = p.get("tag") === "MANE Select";
  await speciesReady;
  if (p.get("loc")) run(p.get("loc"), p.get("to") ?? "", push, p.get("taxon") ?? "", p.get("assembly") ?? "", p.get("db") ?? "");
}
window.addEventListener("popstate", () => fromUrl(false));
fromUrl(false);
