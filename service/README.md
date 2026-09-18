# @togocoord/service

```sh
node service/src/serve.ts --port 8080 human_genome.sqlite human_rna.sqlite human_ensembl.sqlite human_uniprot.sqlite sifts_human.sqlite
# Web UI: http://127.0.0.1:8080/   REST API: /v1/convert?loc=uniprot:P07203:49&to=genome
```

TogoCoord のサービス層（v0.1）。複数の保存先（`@togocoord/ingest` で作る SQLite）をまたいで経路を探索し、Location ID を変換する。

- 規則と検証結果: [../docs/spec-service.md](../docs/spec-service.md)

```ts
import { parseLocationId } from "@togocoord/core";
import { convert, StoreSet } from "@togocoord/service";

const stores = new StoreSet(["human_genome.sqlite", "human_rna.sqlite"]);
const ctx = stores.context();
for (const r of convert(stores, parseLocationId("refseq:NP_036366.3:20", ctx), { to: { category: "genome" } }, ctx)) {
  console.log(r.id, r.cost, r.approximate, r.path.map((s) => `${s.kind}:${s.to}`).join(" > "));
}
// refseq:NW_011332691.1:complement(121036..121038)  1 false  annotation:refseq:NW_011332691.1
// refseq:NC_000003.12:complement(72446594..72446596) 2 false  annotation:refseq:NM_012234.7 > alignment:refseq:NC_000003.12
```

結果は、変換先の配列ごとに1件ずつ、コストの安い順に並ぶ。上の例では、代替配列（NW_）には検証済みの CDS を直接たどって到達し、GRCh38 の主染色体（NC_）には RefSeq 転写産物と `cDNA_match` を経由して到達している。

検証: `node service/bench/verify-search.ts GENOME.fna PROTEINS.faa.gz SAMPLES SEED [--exception-only] STORE.sqlite...`
