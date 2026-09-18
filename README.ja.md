# TogoCoord

[English](README.md) | 日本語

TogoCoord は、ゲノム・転写産物・タンパク質・立体構造など、生命科学の異なるレイヤーの配列座標を相互に変換する仕組みです。位置は INSDC の location 記法に基づく **Location ID**（例: `refseq:NC_000003.12:complement(49358132..49358134)`、`uniprot:P07203:49`）で表し、FALDO JSON-LD としても出力します。

- 一次リポジトリ（GenBank / GFF3 + FASTA）を基盤にして、全生物種に一般化する。ヒトやマウスの豊富なリソース（Ensembl、UniProt、SIFTS、MANE など）は拡張として加える。
- 座標の対応は「location = 写像」として、ブロック列の演算（写像・反転・合成）で扱う。
- 配列が完全に一致するもの（refget ダイジェスト）は同じものとみなし、完全一致の経路を優先する。
- 同じ種のアセンブリ間（GRCh37 ↔ GRCh38 など）は必要に応じて自動でまたぎ、種をまたぐ変換は指定したときだけ行う。
- すべての対応は、取り込み時に実際の配列で自己検証する。

## 構成

| ディレクトリ | 内容 |
|---|---|
| `core/` | Location ID のパース・正規形・FALDO JSON-LD、ブロック列による写像（依存なし） |
| `ingest/` | GenBank / GFF3 / FASTA / SIFTS / MANE / UCSC chain / PAF / BED / NCBI assembly report のアダプタ、自己検証、SQLite（R*Tree）の保存先、CLI `togocoord-ingest` |
| `service/` | 複数の保存先をまたぐ経路探索、REST API、Web UI（`togocoord-serve`） |
| `docs/` | 設計書と仕様（[design.md](docs/design.md)、[spec-core.md](docs/spec-core.md)、[spec-ingest.md](docs/spec-ingest.md)、[spec-service.md](docs/spec-service.md)、[scaling.md](docs/scaling.md)） |
| `poc/` | コンセプト検証版の実装（Web UI と SPARQList。参考） |

## 使い方

Node.js 23.6 以上（TypeScript を型の除去だけで直接実行します）。

```sh
npm install
npm test            # core / ingest / service のテスト

# 保存先を作る（データはリポジトリに含めていません。NCBI などから取得してください）
node ingest/src/cli.ts --db human.sqlite --assembly-report GRCh38.p14_assembly_report.txt \
  --fasta GRCh38.p14_genomic.fna --fasta GRCh38.p14_protein.faa.gz GRCh38.p14_genomic.gff.gz
node ingest/src/cli.ts --db human_uniprot.sqlite UP000005640_9606.fasta.gz

# REST API と Web UI
node service/src/serve.ts --port 8080 human.sqlite human_uniprot.sqlite
# http://127.0.0.1:8080/  ·  /v1/convert?loc=uniprot:P07203:49&to=genome
# 種をまたぐ（UCSC chain を読み込んだとき）: /v1/convert?loc=uniprot:P07203:49&to=protein&db=uniprot&taxon=10090
# アセンブリの配列名で入力（GRCh37 と chain を読み込んだとき）: /v1/convert?loc=hg19:chr7:140453136&to=genome&assembly=GRCh38
```

リバースプロキシでサブディレクトリに置く場合は、[spec-service §6](docs/spec-service.md) を参照してください。

取り込める入力と CLI のオプションは [ingest/README.md](ingest/README.md)、API は [docs/spec-service.md](docs/spec-service.md) を参照してください。

## ライセンス

MIT（[LICENSE](LICENSE)）。
