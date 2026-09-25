# @togocoord/core

[English](README.md) | 日本語

TogoCoord のコアライブラリ（v0.1）。Location ID のパース・正規化と、ブロック列による座標写像を提供する。I/O は持たない。

- 設計: [../docs/design.md](../docs/design.ja.md)
- 規則: [../docs/spec-core.md](../docs/spec-core.ja.md)

## 使い方

```ts
import { cdsMapping, createContext, formatLocationId, mapLocation, parseLocationId } from "@togocoord/core";

const ctx = createContext();
// SARS-CoV-2 ORF1ab (NC_045512.2): CDS with a -1 ribosomal slippage
const toGenome = cdsMapping({
  protein: "refseq:YP_009724389.1",
  cds: parseLocationId("refseq:NC_045512.2:join(266..13468,13468..21555)", ctx),
  aaLength: 7096,
});

mapLocation(parseLocationId("refseq:YP_009724389.1:4401..4402", ctx), toGenome, ctx)
  .targets.map((t) => formatLocationId(t.location, ctx));
// => ["refseq:NC_045512.2:join(13466..13468,13468..13470)"]

mapLocation(parseLocationId("refseq:NC_045512.2:13468", ctx), toGenome.inverse(), ctx)
  .targets.map((t) => formatLocationId(t.location, ctx));
// => ["refseq:YP_009724389.1:join(4401c3,4402c1)"]   ("never" モードでは join(4401,4402))
```

## 構成

| ファイル | 内容 |
|---|---|
| `src/parser.ts` | location の文法（spec-core §3.2） |
| `src/namespace.ts` | 名前空間の登録・accession の検査・既定の単位（§3.4） |
| `src/location.ts` | 意味値（セグメント列）への変換、正規形での出力（§3.3, §3.5） |
| `src/mapping.ts` | ブロック、`invert`、`compose`、`mappingFromLocation`、`cdsMapping`（§4） |
| `src/convert.ts` | `mapLocation`（§5） |
| `test/corpus/` | フェーズ0のテストデータ（実データ4件、合成データ3件） |
| `test/oracle.ts` | 残基単位の素朴な実装（テストの正解） |

## 開発

Node.js 24 以上（TypeScript を型の除去だけで直接実行する）。

```sh
npm install
npm test          # node --test（テストデータ + property テスト + 単体テスト）
npm run typecheck # tsc --noEmit
```
