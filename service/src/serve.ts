#!/usr/bin/env node
// togocoord-serve [--port 8080] [--host 127.0.0.1] [--base URL] STORE.sqlite...
import { createApi } from "./api.ts";
import { StoreSet } from "./stores.ts";

const args = process.argv.slice(2);
const USAGE = "usage: togocoord-serve [--port 8080] [--host 127.0.0.1] [--base URL] STORE.sqlite...\n";
let port = 8080;
let host = "127.0.0.1";
let base: string | undefined;
const paths: string[] = [];
for (let i = 0; i < args.length; i++) {
  const a = args[i]!;
  if (a === "--port") port = Number(args[++i]);
  else if (a === "--host") host = args[++i]!;
  else if (a === "--base") base = args[++i];
  else if (a === "-h" || a === "--help") {
    process.stdout.write(USAGE);
    process.exit(0);
  } else if (a.startsWith("--")) {
    process.stderr.write(`unknown option ${a}\n${USAGE}`);
    process.exit(1);
  } else paths.push(a);
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
