# データの取得と構築

[English](data.md) | 日本語

デモの保存先（SQLite）を、一次データの取得から作り直す手順。定義と実行は `scripts/data.ts`、配信する順番は `scripts/stores.txt` にある。データそのものは GitHub に置かない（`data/` は `.gitignore` で除外）。

## 必要なもの

- Node.js 23.6 以上（リポジトリの `npm install` 済み）
- ゼニゴケの v3.1 ↔ v7.1 のアライメントだけ、[minimap2](https://github.com/lh3/minimap2)（2.31 で作成。`brew install minimap2`）
- ディスク約 30GB（取得したファイルと展開したゲノムで約 18GB、保存先で約 9GB）。メモリ 8GB 程度（ヒトの構築で最大約 6GB）

## 使い方

```sh
node scripts/data.ts list                  # データセットの一覧、グループ、構築済みかどうか（大きさ）
node scripts/data.ts download [名前...]     # 必要なファイルを取得（名前を省くと全部）
node scripts/data.ts build [名前...]        # 保存先を作る（ないものだけ。--force で作り直す）
node scripts/data.ts serve [--port 8080] [--base URL] [--host 0.0.0.0]
                                           # scripts/stores.txt の順で、作った保存先を配信する
```

- 名前は、保存先（`human_rna`）かグループ（`human`）。`build` は、必要なファイルがなければ先に取得する。
- 置き場所は `$TOGOCOORD_DATA`（既定はリポジトリの `data/`）。`raw/` に取得したファイル（配布元のファイル名のまま）、`work/` に途中のファイル（アライメントの PAF）、`stores/` に保存先、`logs/` に構築のログ。
- 大きな gzip の FASTA は、取り込み時に `raw/` の中に展開したファイルと `.fai` を作る（spec-ingest §8）。
- すべてを作ると、取得済みのファイルからで約13分（8コアの Mac、最大メモリ約7GB）。取得は回線次第で、合計約6GB。
- 例: ヒトだけを作って配信する。`node scripts/data.ts build human && node scripts/data.ts serve`（`stores.txt` にあってまだ作っていない保存先は飛ばす）。

## データセット

| グループ | 保存先 | 入力（配布元） | 大きさ |
|---|---|---|---|
| human | `human` | RefSeq GRCh38.p14（GCF_000001405.40）の GFF3、ゲノム・RNA・タンパク質の配列、assembly report（NCBI） | 1.6GB |
| | `human_rna` | RefSeq GRCh38.p14 の RNA の GenBank | 650MB |
| | `human_ensembl` | Ensembl release 116 の GFF3 とタンパク質配列 | 3.1GB |
| | `human_uniprot` | UniProt 参照プロテオーム UP000005640（current_release） | 57MB |
| | `human_mane` | MANE v1.5 の summary と RNA 配列 | 23MB |
| | `grch38_names` | GRCh38.p14 の assembly report（配列名） | 0.3MB |
| | `human_uniprot_alignments` | 同一配列のない UniProt のエントリと、RefSeq・Ensembl のタンパク質のアライメント（T2。UniProt の `HUMAN_9606_idmapping_selected.tab.gz`） | 5MB |
| grch37 | `grch37` | GRCh37.p13（GCF_000001405.25）の assembly report とゲノム配列 | 0.3MB |
| | `chain_hg19ToHg38`、`chain_hg38ToHg19` | UCSC liftOver chain | 1MB、10MB |
| mouse | `mouse`、`mouse_rna`、`mouse_uniprot`、`grcm39_names`、`mouse_uniprot_alignments` | RefSeq GRCm39（GCF_000001635.27）、UniProt UP000000589 | 930MB、390MB、21MB、0.1MB |
| grcm38 | `grcm38`、`chain_mm10ToMm39`、`chain_mm39ToMm10` | GRCm38.p6（GCF_000001635.26）、UCSC liftOver chain | 0.2MB、0.3MB、0.5MB |
| human_mouse | `chain_hg38ToMm39`、`chain_mm39ToHg38` | UCSC liftOver chain | 158MB、156MB |
| arabidopsis | `arabidopsis`、`arabidopsis_rna`、`arabidopsis_uniprot`、`arabidopsis_uniprot_alignments` | RefSeq TAIR10.1（GCF_000001735.4）、UniProt UP000006548 | 300MB、174MB、14MB |
| tair10 | `tair10`、`chain_tair10.1ToTair10`、`tair10_to_tair10.1` | 一つ前の RefSeq の版 TAIR10（GCF_000001735.3）の assembly report とゲノム配列、UCSC GenArk の chain（TAIR10.1 → TAIR10 だけ）、逆向きは minimap2 のアライメント | 0.1MB、0.1MB、0.1MB |
| marchantia | `marchantia_v71` | INSDC MpTak_v7.1（GCA_039105155.1）の GenBank | 72MB |
| | `marchantia` | INSDC MpTak v3.1（GCA_003032435.1）の GenBank | 80MB |
| | `marchantia_uniprot` | UniProt UP000244005 | 7MB |
| | `mp_v31_to_v71`、`mp_v71_to_v31` | 両方のゲノム配列を minimap2 でアライメントした PAF（`work/` に作る） | 3MB、3MB |
| sifts | `sifts` | SIFTS `uniprot_segments_observed.tsv.gz`（PDBe）、`pdb_seqres.txt.gz`（wwPDB）、上の UniProt すべて | 180MB |
| fanta | `fanta_human_hg38`、`fanta_mouse_mm10` | fanta.bio CRE v1.2.1 の BED | 197MB、116MB |

各グループの結果（件数、自己検証、変換の検証）は spec-ingest と spec-service にある。

**2026-09-18 の作り直しの確認**: このスクリプトで全29の保存先を作り直し、以前の手作業の保存先と、配列・edge・注釈の数と自己検証の結果が一致することを確かめた。違いは chain の4つだけで、その後の改良による（両側のアセンブリの配列を記録する、パッチの UCSC 名を読み替えて取り込む chain が 5〜277 本増えた）。増えた chain のうち、RefSeq に収録されていない未配置 scaffold（`KI270752.1` など）に着くものは、RefSeq のゲノム配列で照合できないので自己検証が `skipped` になる。

**TAIR10 と TAIR10.1**: 核の染色体と葉緑体は同じ accession・同じ配列で、違うのはミトコンドリアのゲノム（TAIR10 `NC_001284.2`、TAIR10.1 `NC_037304.1`）だけ。核の位置は、どちらのアセンブリでもそのまま答えになる（spec-service §2.2）。

**入れていないもの**: Col-CEN v1.2（セントロメアまでつながった Col-0 のアセンブリ。UCSC GenArk に TAIR10.1 との chain がある）は、INSDC の accession がなく GitHub（schatzlab/Col-CEN）でだけ配布されているので、配列の鍵を決められず、入れていない。

**版について**: NCBI のアセンブリ、Ensembl（release-116）、MANE（release_1.5）、fanta.bio（v1.2.1）は版を固定した URL から取る。UniProt（current_release）、SIFTS、PDB の配列は、配布元が同じ URL で更新するので、取得した時期によって中身が変わる。

## 再現のための記録

保存先の `meta` には、次を記録する（Web UI の Loaded data に表示）。

- 入力ファイルの MD5（`inputs_md5`）と、`--fasta` で与えた配列ファイルの MD5（`sequences_md5`）
- 構築した TogoCoord の commit（`togocoord`。未 commit の変更があれば `+local changes`）
- 入力の作り方（`method`。アライメントなら minimap2 の版と引数、PAF の絞り方）

同じ入力（MD5 が同じ）と同じ commit で `build --force` すれば、同じ保存先ができる。

## 配信とデプロイ

`serve` は `togocoord-serve --stores scripts/stores.txt --store-dir data/stores --skip-missing` を実行する。`stores.txt` の順番は結果に影響する（同じ配列を複数の保存先が記録していれば、先の保存先の値を使う。アノテーションのあるアセンブリを先に、アライメントを最後に置く）。

別のサーバに置くときは、リポジトリと `data/stores/`（保存先だけでよい。`raw/` と `work/` は不要）をコピーし、次のように起動する。

```sh
node --no-warnings service/src/serve.ts --host 127.0.0.1 --port 8080 --base https://example.org/togocoord/ \
  --stores scripts/stores.txt --store-dir /path/to/stores
```

サブディレクトリに置く場合は spec-service §6 を参照。`*.sqlite` をまとめて渡すと、余計なファイルも読み込み、順番もアルファベット順になるので、`--stores` を使う。

## 新しいデータセットを加える

`scripts/data.ts` の `STORE_LIST` に、名前、グループ、`togocoord-ingest` の引数を加える。取得するファイルは `raw(URL)`（NCBI のアセンブリは `ncbi(assembly(accession, name), suffix)`）で書くと、`download` の対象になる。配信するなら `scripts/stores.txt` の適切な位置に名前を加える。
