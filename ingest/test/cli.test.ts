import assert from "node:assert/strict";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import { it } from "node:test";

it("togocoord-ingest writes JSON Lines and a validation summary", () => {
  const cli = fileURLToPath(new URL("../src/cli.ts", import.meta.url));
  const gb = fileURLToPath(new URL("./fixtures/NC_045512.2.gb", import.meta.url));
  const run = spawnSync(process.execPath, ["--no-warnings", cli, gb], { encoding: "utf8" });
  assert.equal(run.status, 0, run.stderr);
  const records = run.stdout.trim().split("\n").map((l) => JSON.parse(l) as { record: string; type?: string });
  const count = (r: string) => records.filter((x) => x.record === r).length;
  assert.equal(count("sequence"), 13);
  assert.equal(count("edge"), 12);
  assert.equal(count("annotation"), 56);
  assert.ok(records.filter((x) => x.record === "annotation").every((x) => typeof x.type === "string" && x.type !== "annotation"));
  assert.match(run.stderr, /edges 12 \{"ok":12\}/);
});
