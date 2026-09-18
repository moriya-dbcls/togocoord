# @togocoord/ingest

TogoCoord のアダプタ（v0.2）。GenBank/GenPept のフラットファイルと GFF3 を読み、配列、写像の edge、アノテーションを出力する。取り込み時に自己検証（翻訳の照合、アライメントの同一塩基数の照合）を行う。

- 規則: [../docs/spec-ingest.md](../docs/spec-ingest.md)
- コア: [../core](../core)

## CLI

```sh
# 小さなファイル: JSON Lines を標準出力へ
node ingest/src/cli.ts NC_012920.1.gb > mt.jsonl

# ゲノム全体: SQLite へ（.gz をそのまま読む。FASTA は .fai で必要な区間だけ読む）
node --max-old-space-size=512 ingest/src/cli.ts \
  --db human.sqlite --fasta GRCh38.fna --fasta GRCh38_protein.faa.gz GRCh38_genomic.gff.gz
# stderr: sequences 137512, edges 505277 {"skipped":359831,"ok":142549,"mismatch":2897}, ... (73 s)
#         mismatches: 0 unexplained, 2897 with an INSDC /exception
```

| オプション | 内容 |
|---|---|
| `--db FILE` | SQLite に保存する（なければ JSON Lines を出力する）。既存のファイルは `--overwrite` で置き換える |
| `--fasta FILE` | 自己検証に使う配列（ゲノム、転写産物、タンパク質）。何度でも指定できる。64MB を超えるか `.fai` があるファイルは、ランダムアクセスで読む |
| `--all-annotations` | exon も annotation として保存する |

出力する JSON Lines では、各行の `record` が `sequence` / `edge` / `annotation` / `warning` のいずれかになる。

## 保存先の検証

```sh
node ingest/bench/verify-store.ts human.sqlite GRCh38.fna GRCh38_protein.faa.gz 10000
# {"residues":10000,"genomicHits":10616,"agree":10594,"disagree":0,"disagreeOnExceptionEdges":22,...,"roundTrip":10616,"roundTripMissing":0,...}
# protein -> genome: p50 0.06 ms, p95 0.62 ms, p99 0.95 ms
```

ランダムに選んだタンパク質の残基を、保存先を使ってゲノム上に変換し、そのコドンを翻訳して公開配列と比べる。さらにゲノムから逆に変換して、元の残基に戻るかも調べる。結果の読み方は [../docs/scaling.md](../docs/scaling.md) §7 を参照。

## ライブラリ

```ts
import { ingestGenBank, edgeMapping } from "@togocoord/ingest";
import { createContext, mapLocation, parseLocationId, formatLocationId } from "@togocoord/core";

const result = ingestGenBank(text);
const nd6 = result.edges.find((e) => e.from === "refseq:YP_003024037.1")!;
const units = new Map(result.sequences.map((s) => [s.ref, s.unit]));
const ctx = createContext({ units: (ref) => units.get(ref) });
mapLocation(parseLocationId("refseq:YP_003024037.1:174", ctx), edgeMapping(nd6), ctx)
  .targets.map((t) => formatLocationId(t.location, ctx));
// => ["refseq:NC_012920.1:complement(14152..14154)"]
```

## 構成

| ファイル | 内容 |
|---|---|
| `src/gbff.ts`、`src/gff3.ts`、`src/fasta.ts` | パーサ |
| `src/adapter-gbff.ts`、`src/adapter-gff3.ts` | アダプタ（spec-ingest §3, §4） |
| `src/validate.ts` | 自己検証（§6）とタンパク質長の推定（§5） |
| `src/sequence.ts` | 配列の取得、相補鎖、翻訳、refget ダイジェスト |
| `src/stream.ts` | ストリーミング読み込み（`.gz` 対応）、`JsonlSink` |
| `src/fasta-index.ts` | `.fai` の作成と、FASTA へのランダムアクセス |
| `src/store.ts` | SQLite への書き込み（`SqliteSink`）と問い合わせ（`TogoCoordStore`） |
| `bench/verify-store.ts` | 公開タンパク質配列による照合と、問い合わせ時間の計測 |
| `src/genetic-codes.ts` | NCBI の `gc.prt` から生成した翻訳表 |
| `test/fixtures/` | NCBI の実データ（2026-09-18 取得） |
