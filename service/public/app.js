// TogoCoord web UI: a thin client of the REST API (spec-service §6). State lives in the URL (?loc=&to=&codon=).
const $ = (sel, root = document) => root.querySelector(sel);

const EXAMPLES = [
  ["UniProt residue", "uniprot:P07203:49"],
  ["genome codon", "refseq:NC_000003.12:complement(49358132..49358134)"],
  ["structure residue", "pdb:2F8A.A:59"],
  ["Ensembl protein range", "ensembl:ENSP00000407375.1:40..60"],
  ["RefSeq protein", "refseq:NP_036366.3:20"],
  ["whole protein", "refseq:NP_000572.2"],
  ["overlapping genes (mtDNA)", "refseq:NC_012920.1:8527..8529"],
  ["human → mouse protein", "refseq:NP_000572.2:49", "protein"],
  ["mouse → human genome", "refseq:NC_000075.7:106312500..106312550", "genome"],
];

const form = $("#query");
const locInput = $("#loc");
const toSelect = $("#to");
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

async function api(path) {
  const res = await fetch(path, { headers: { accept: "application/json" } });
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
      const res = await fetch(`/v1/location/faldo?${qs({ loc: id })}`);
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

async function run(loc, to, push = true) {
  loc = loc.trim();
  locInput.value = loc;
  toSelect.value = to ?? "";
  $("#error").hidden = true;
  if (!loc) return;
  const codon = codonBox.checked ? "" : "never";
  const tag = maneBox.checked ? "MANE Select" : "";
  if (push) history.pushState(null, "", `?${qs({ loc, to, codon, tag })}`);
  try {
    const [info, conv] = await Promise.all([
      api(`/v1/location?${qs({ loc, codon })}`),
      api(`/v1/convert?${qs({ loc, to, codon, tag })}`),
    ]);
    $("#input-id").textContent = info.id;
    $("#input-kind").textContent = `${info.unit === "aa" ? "protein" : "nucleotide"}${info.kind === "order" ? " · order" : ""}`;
    $("#segments").replaceChildren(...segmentRows(info.segments));
    const inputCard = $("#input .card");
    $(".extra", inputCard).hidden = true;
    $("#input").hidden = false;
    inputTaxon = conv.inputTaxon;
    $("#results").replaceChildren(...conv.results.map(renderResult));
    $("#count").textContent = conv.results.length
      ? `(${conv.results.length}${conv.truncated ? ", truncated — narrow the target or the input" : ""})`
      : tag ? "— none with this tag" : "— none reachable";
    $("#output").hidden = false;
  } catch (err) {
    $("#input").hidden = true;
    $("#output").hidden = true;
    showError(err, loc);
  }
}

form.addEventListener("submit", (e) => {
  e.preventDefault();
  run(locInput.value, toSelect.value);
});
toSelect.addEventListener("change", () => locInput.value && run(locInput.value, toSelect.value));
codonBox.addEventListener("change", () => locInput.value && run(locInput.value, toSelect.value));
maneBox.addEventListener("change", () => locInput.value && run(locInput.value, toSelect.value));
$("#input .card").addEventListener("click", (e) => {
  const action = e.target.dataset?.action;
  if (action) toggleExtra($("#input .card"), action, $("#input-id").textContent);
});
$("#examples").append(
  "Examples: ",
  ...EXAMPLES.map(([label, id, to]) => el("button", { type: "button", class: "small", title: id, onclick: () => run(id, to ?? toSelect.value) }, label)),
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

function fromUrl(push = false) {
  const p = new URLSearchParams(location.search);
  if (p.get("view") === "data") {
    showData().catch((err) => showError(err));
    return;
  }
  showConvert();
  codonBox.checked = p.get("codon") !== "never";
  maneBox.checked = p.get("tag") === "MANE Select";
  if (p.get("loc")) run(p.get("loc"), p.get("to") ?? "", push);
}
window.addEventListener("popstate", () => fromUrl(false));
fromUrl(false);
