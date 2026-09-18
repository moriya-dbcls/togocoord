#!/usr/bin/env node
// togocoord-serve [--port 8080] [--host 127.0.0.1] [--base URL] [--stores LIST [--store-dir DIR] [--skip-missing]] [STORE.sqlite...]
// Stores are read in the order given (the first store's value wins where several describe one sequence). A LIST file
// names them one per line (`#` comments; `name` means name.sqlite), relative to --store-dir (default: the LIST's
// directory). See docs/data.md.
import { existsSync, readFileSync } from "node:fs";
import { dirname, isAbsolute, join } from "node:path";
import { createApi } from "./api.ts";
import { StoreSet } from "./stores.ts";

const args = process.argv.slice(2);
const USAGE =
  "usage: togocoord-serve [--port 8080] [--host 127.0.0.1] [--base URL] [--stores LIST [--store-dir DIR] [--skip-missing]] [STORE.sqlite...]\n";
let port = 8080;
let host = "127.0.0.1";
let base: string | undefined;
let list: string | undefined;
let storeDir: string | undefined;
let skipMissing = false;
const paths: string[] = [];
for (let i = 0; i < args.length; i++) {
  const a = args[i]!;
  if (a === "--port") port = Number(args[++i]);
  else if (a === "--host") host = args[++i]!;
  else if (a === "--base") base = args[++i];
  else if (a === "--stores") list = args[++i];
  else if (a === "--store-dir") storeDir = args[++i];
  else if (a === "--skip-missing") skipMissing = true;
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

const stores = new StoreSet(paths);
const server = createApi(stores, { ...(base && { base }) });
server.listen(port, host, () => {
  process.stderr.write(`togocoord: ${paths.length} store(s) on http://${host}:${port}/\n`);
});
const shutdown = () => server.close(() => (stores.close(), process.exit(0)));
process.on("SIGINT", shutdown);
process.on("SIGTERM", shutdown);
