# いまの状況

[English](status.md) | [日本語](status.ja.md)

再開のための記録。何が動いていて、何が判断待ちで、次に何から手を付けるか。理由は各 spec にあるので、ここには決まって
いないことだけを書く。最終更新 2026-10-19、`6d90ec6` 時点。

## 状態

- `main` は `6d90ec6`。作業ツリーはクリーンで、`github.com/moriya-dbcls/togocoord` に push 済み。
- テストは通る: core 119、ingest 71、service 62。Node 24（`npm test`、`npm run typecheck`）。
- `data/stores` に51の保存先（10GB）。一覧は `scripts/stores.txt`。作り直しは `node scripts/data.ts build`
  （docs/data.md）。データは commit しない。

## デプロイがリポジトリより古い

**504 の修正が入っていない。** nginx の後ろの環境で 504 が出ていた件は、先頭の要求が後続を塞ぐ現象と判明し
（spec-service §15）、`da47887` / `6d90ec6` で直したが、反映には pull と再起動が要る。

```sh
git pull
node --no-warnings service/src/serve.ts --host 127.0.0.1 --port 8080 --base https://<host>/<path>/ \
  --workers auto --log --stores scripts/stores.txt --store-dir /path/to/stores
```

遅い変換1本でサービス全体が止まるのを防ぐのが `--workers`、次に同じことが起きたときに原因を特定する手段が `--log`。
nginx の `proxy_read_timeout` は上げないこと。1分を超える要求はこちら側の不具合である。

同じ pull で入る未デプロイの変更: `/v1/meta` の修正（`c8c0a96`。10.2秒 → 10ミリ秒未満。これだけでもページを開く
たびに 504 の原因になりうる）、名前空間の推定、JSON リンク、favicon、`coverage` と `minMatch`。

前回データをコピーした後に作った、サーバに要る保存先4つ（95MB）: `marchantia_mpv4`、`marchantia_uniprot_v4`、
`mp_v4_to_hub`、`mp_hub_to_v4`。

**未確認**: デプロイ先がどの URL で応答するか。2026-10-08 にこちらから叩いたところ、
`https://sparql-support.dbcls.jp/togocoord/v1/meta` は 301、`/togocoord_v1/` は 302 の先が 404 で、どちらも測れな
かった。デプロイを試すときは URL の確認から。

## 判断待ち

- **オーソログ。** 実測では、OMA の Arabidopsis↔Marchantia は 9,668対あり、`sequence_md5` を持つのでこちらの配列の
  鍵とそのまま突き合わせられる。標本の同一性は50%前後。ユーザの立場（2026-09-25）は、座標サービスの価値は変異や
  ドメイン、二次構造を辿れることにあり、BLAST や OMA が既に答えることを言い直すことではない、というもの。したがって
  オーソログ対を雑にアラインメントするのは筋が悪い。オーソログの辺を作るなら、共有の座標系（InterPro/Pfam の HMM の
  match state、または構造の重ね合わせ）を経由し、その根拠が結果に残る形にすべき。未着手。design.md §12 では出典自体
  が未決のまま。
- **残基レベルの注釈（こちらの方が有望）。** UniProt は配列としてしか取り込んでいないので、活性部位、修飾残基、
  ドメイン、二次構造、変異がどれも無い。提案した順序は、UniProt の feature → 構造由来（DSSP、AlphaFold の pLDDT）
  → 変異（TogoVar / ClinVar）。いずれも既存の配列に対する注釈として入る。あわせて `/v1/annotations` を、到達できる
  他の層の注釈まで集めるよう広げる。未着手。
- **`idmapping_selected.tab` が無い種（ゼニゴケ）の T2。** UniProt の REST から UniRef90 のクラスタを取れば、同一
  配列では到達できない UP000077202 のエントリのうち 6,047件に候補を作れる。未着手。Mp_v4 を入れたことで対象は
  約400件（cDNA 由来、属内の他種、ミトコンドリアの45件）まで縮んでいる。

## 決着済み（蒸し返さない）

- 1つの位置に対する複数の対応は返さない。`-multiple` 相当は設けない（spec-service §2）。
- 保存先は残基も配列も持たない。差分の位置だけを記録する（spec-ingest §19）。
- `MpTak1_v7.1` は入れない。配列も CDS のタンパク質も中心と同じ accession（spec-ingest §20）。
- `minMatch` に既定値は置かない。UCSC liftOver と違い、ここでは部分的な結果が通常のことだから（spec-service §4）。
- BRCA1 の CRE の件で chain score による順位付けを検討したが、効かない（誤った座位の方が score が高い）。

## 細かい未決事項

- `togocoord-serve` には流量制限も認証も同時実行数の上限も無い。リバースプロキシ側で持つべき。
- Col-CEN v1.2（シロイヌナズナ）は INSDC の accession が無く、配列の鍵を決められない（spec-ingest §20）。
- design.md §12: `faldo:codonPosition` の FALDO 開発者への提案、ドメイン名と identifiers.org の prefix、どこで
  運用するか、経路が見つからなかった問い合わせを記録するかどうか。
