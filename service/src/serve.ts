#!/usr/bin/env node
// togocoord-serve [--port 8080] [--host 127.0.0.1] [--base URL] [--workers N] [--log] [--budget MS]
//                  [--max-input-length N] [--max-results N] [--max-batch N] [--max-annotations N]
//                  [--stores LIST [--store-dir DIR] [--skip-missing]] [STORE.sqlite...]
// Stores are read in the order given (the first store's value wins where several describe one sequence). A LIST file
// names them one per line (`#` comments; `name` means name.sqlite), relative to --store-dir (default: the LIST's
// directory). See docs/data.md.
import cluster from "node:cluster";
import { existsSync, readFileSync } from "node:fs";
import { availableParallelism } from "node:os";
import { dirname, isAbsolute, join } from "node:path";
import { createApi } from "./api.ts";
import { StoreSet } from "./stores.ts";

const args = process.argv.slice(2);
const USAGE =
  "usage: togocoord-serve [--port 8080] [--host 127.0.0.1] [--base URL] [--workers N] [--log] [--budget MS]\n" +
  "                       [--max-input-length N] [--max-results N] [--max-batch N] [--max-annotations N]\n" +
  "                       [--stores LIST [--store-dir DIR] [--skip-missing]] [STORE.sqlite...]\n";
let port = 8080;
let host = "127.0.0.1";
let base: string | undefined;
let list: string | undefined;
let storeDir: string | undefined;
let skipMissing = false;
let workers = 1;
let log = false;
const limits: { searchBudgetMs?: number; maxInputLength?: number; maxResults?: number; maxBatch?: number; maxAnnotations?: number } = {};
const paths: string[] = [];
for (let i = 0; i < args.length; i++) {
  const a = args[i]!;
  if (a === "--port") port = Number(args[++i]);
  else if (a === "--host") host = args[++i]!;
  else if (a === "--base") base = args[++i];
  else if (a === "--stores") list = args[++i];
  else if (a === "--store-dir") storeDir = args[++i];
  else if (a === "--skip-missing") skipMissing = true;
  // SQLite reads are synchronous, so one process answers one request at a time: a slow conversion delays every other
  // request behind it (spec-service §15). Workers are separate processes sharing the (read-only) stores.
  // `auto` takes half the cores (8 workers on 8 cores measured slower than 4, and the host is usually shared).
  else if (a === "--workers") workers = args[++i] === "auto" ? Math.max(1, Math.min(8, availableParallelism() >> 1)) : Number(args[i]);
  else if (a === "--log") log = true;
  else if (a === "--budget") limits.searchBudgetMs = Number(args[++i]);
  else if (a === "--max-input-length") limits.maxInputLength = Number(args[++i]);
  else if (a === "--max-results") limits.maxResults = Number(args[++i]);
  else if (a === "--max-batch") limits.maxBatch = Number(args[++i]);
  else if (a === "--max-annotations") limits.maxAnnotations = Number(args[++i]);
  else if (a === "-h" || a === "--help") {
    process.stdout.write(USAGE);
    process.exit(0);
  } else if (a.startsWith("--")) {
    process.stderr.write(`unknown option ${a}\n${USAGE}`);
    process.exit(1);
  } else paths.push(a);
}
if (list) {
  const dir = storeDir ?? dirname(list);
  const listed = readFileSync(list, "utf8")
    .split("\n")
    .map((l) => l.replace(/#.*/, "").trim())
    .filter(Boolean)
    .map((n) => (isAbsolute(n) ? n : join(dir, n.endsWith(".sqlite") ? n : `${n}.sqlite`)));
  const missing = listed.filter((p) => !existsSync(p));
  if (missing.length && !skipMissing) {
    process.stderr.write(`missing stores (build them, or pass --skip-missing):\n${missing.map((p) => `  ${p}\n`).join("")}`);
    process.exit(1);
  }
  if (missing.length) process.stderr.write(`skipping ${missing.length} missing store(s): ${missing.map((p) => p.split("/").pop()).join(", ")}\n`);
  paths.unshift(...listed.filter((p) => existsSync(p)));
}
if (paths.length === 0) {
  process.stderr.write(USAGE);
  process.exit(1);
}

if (workers > 1 && cluster.isPrimary) {
  process.stderr.write(`togocoord: ${workers} workers on http://${host}:${port}/\n`);
  for (let i = 0; i < workers; i++) cluster.fork();
  // A worker that dies (out of memory, a bug) is replaced, so the service keeps answering.
  cluster.on("exit", (w, code, signal) => {
    if (!stopping) {
      process.stderr.write(`togocoord: worker ${w.process.pid} exited (${signal ?? code}); starting another\n`);
      cluster.fork();
    }
  });
  let stopping = false;
  for (const sig of ["SIGINT", "SIGTERM"] as const) {
    process.on(sig, () => {
      stopping = true;
      for (const w of Object.values(cluster.workers ?? {})) w?.kill(sig);
      process.exit(0);
    });
  }
} else {
  const stores = new StoreSet(paths);
  const server = createApi(stores, { ...(base && { base }), ...limits });
  if (log) {
    // One line per request with how long it took: a request over the budget is the one delaying the others.
    // Before the API's own handler, which runs synchronously: registered after it, the clock would start once the
    // work was already done and every request would look instant.
    server.prependListener("request", (req, res) => {
      const started = performance.now();
      res.on("finish", () =>
        process.stderr.write(`${new Date().toISOString()} ${res.statusCode} ${Math.round(performance.now() - started)}ms ${req.method} ${req.url}\n`),
      );
    });
  }
  // Warm the caches that /v1/meta and the scope rules build (species, assemblies, crossings). Cold, that work takes
  // seconds, and it would be paid by the first request each worker answers - right after a restart, when several
  // arrive at once.
  const warmed = performance.now();
  stores.species();
  stores.assemblies();
  stores.crossings();
  stores.tagSpecies();
  server.listen(port, host, () => {
    const who = cluster.isWorker ? `worker ${process.pid}` : `${paths.length} store(s)`;
    process.stderr.write(`togocoord: ${who} on http://${host}:${port}/ (warmed in ${Math.round(performance.now() - warmed)} ms)\n`);
  });
  const shutdown = () => server.close(() => (stores.close(), process.exit(0)));
  process.on("SIGINT", shutdown);
  process.on("SIGTERM", shutdown);
}
