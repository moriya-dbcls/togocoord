# TogoCoord 設計書（v0.1 草案）

2026-09-18

## 0. 背景と方針

TogoCoord は、ゲノム・転写産物・タンパク質・立体構造など、生命科学の異なるレイヤーにまたがる配列座標を変換するサービスである。既存のコンセプト検証実装（`html/`、`sparqlist/`）は参考にとどめ、次の方針で作り直す。

| 課題（PoC） | 方針（本設計） |
|---|---|
| Ensembl REST・TogoWS・UniParc など、ヒトやマウスで充実した外部 API に依存している | 一次リポジトリ（GBFF、GFF3+FASTA）を基盤にして、全生物種に一般化する |
| delta 方式の対応表（負の座標や complement の特別扱い）が複雑で、バグの温床になっている | 「location = 写像」という考え方とブロック列の代数に置き換える |
| SPARQList 同士を HTTP で連鎖させていて、遅く、不具合も起きやすい | 純粋なコアライブラリと、それを使うサービスに分ける |
| 種をまたいだ変換が未実装 | compose による連鎖と、alignment edge の追加で対応する |
| ヒトやマウスの豊富なリソース | 一般化したコアの上に、Enrichment 層として後から追加できるようにする |

Location ID は、規則を厳密にしたうえで引き続き使う（§3）。

---

## 1. 用語

| 用語 | 定義 |
|---|---|
| 配列（Sequence） | 実在する残基の文字列。レジストリに登録し、名前空間・accession.version・分子種・長さ・ダイジェストを持つ |
| Location ID | 配列上の位置や区間の集合を表す文字列。INSDC の location 記法に基づく |
| 写像（Mapping） | ある配列の座標から別の配列の座標への対応。ブロック列で表す |
| ブロック（Block） | ギャップを含まない1対1の対応区間 |
| edge | レジストリ上の配列どうしを結ぶ写像に、種別と由来を付けたもの |
| アノテーション | ある配列上に Location ID で置かれた情報（ドメイン、PTM、CRE、バリアントなど） |
| アダプタ | データ源を読み込み、配列・edge・アノテーションを出力するモジュール |

---

## 2. 全体構成

```
┌─ Enrichment（任意・種限定・後から追加）──────────────────────┐
│ MANE/GENCODE, UniProt isoform, SIFTS, AlphaFold DB,           │
│ UCSC chain, Ensembl Compara, HPRC pangenome, FANTOM CAGE,     │
│ ChIP-Atlas, TogoVar/VEP, ドメイン・二次構造 …                 │
└───────────────┬───────────────────────────────────────────────┘
                │ 同じ形式で出力（配列 / edge / アノテーション）
┌─ Core（全生物種・必須）──┴────────────────────────────────────┐
│ アダプタ: GBFF, GFF3+FASTA（INSDC / RefSeq）                    │
│ 正規化ストア: Sequence registry / Mapping edges / Annotations   │
└───────────────┬───────────────────────────────────────────────┘
┌─ コアライブラリ（純粋・I/Oなし）─┴───────────────────────────┐
│ Location ID のパース・正規化 / ブロック演算 / 意味論            │
└───────────────┬───────────────────────────────────────────────┘
┌─ サービス ───────┴────────────────────────────────────────────┐
│ 経路探索・taxonごとのプロファイル / REST API / ワークスペース   │
└───────────────────────────────────────────────────────────────┘
  利用側: Web UI, TogoStanza, SPARQList（API を呼ぶだけ）, CLI
```

- **コアライブラリ**は TypeScript で書く（推奨）。ブラウザ、Node、SPARQList、TogoStanza のすべてで同じコードを動かせるため。CLI としても配布し、大規模データを持つ利用者が自分の環境で使えるようにする。
- **依存は一方向にする**。Core は Enrichment がなくても完全に動く。Enrichment の配列は、必ず Core の accession.version か refget ダイジェストに結び付ける。
- 変換を実行するときに外部 API を呼ばない。外部のデータは、取り込み時にスナップショットとして保存する（遅延取得してキャッシュする方式も可）。

---

## 3. Location ID 仕様

### 3.1 全体形式

```
<namespace>:<accession>[.<version>]:<location>
例: refseq:NM_014739.3:join(233..235,5958..6116)
    insdc:NC_000001.11:complement(join(201..300,401..500))
    uniprot:Q9BYF1-1:60
    refget:SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2:1..100
```

- `namespace` は必須とする。bioregistry などに登録された接頭辞から選ぶ（例: `insdc`、`refseq`、`ensembl`、`uniprot`、`pdb`、`refget`）。
- **version の扱いは名前空間ごとに決める**（§3.7）。version のある体系（INSDC、RefSeq、Ensembl）では、出力に必ず version を付ける。入力に version がない場合は最新版として解決する。UniProt と PDB はバージョンなしで扱う。
- accession に `:` は含めない。最初の `:` までを namespace、次の `:` までを accession とみなしてパースする。
- **assembly は ID に含めない**。INSDC や RefSeq の配列 accession.version はそれだけで一意に決まるため。assembly はレジストリにメタデータとして持つ。
- **座標は、実在する配列にだけ付ける**。`mrna#`、`-pre_mRNA`、`-cDNA` のようなレイヤー指定は使わない。レイヤーはレジストリの分子種から決まる。pre-mRNA はゲノム上の location として表す。

### 3.2 location の文法（INSDC 記法のサブセット）

```ebnf
location   = complement | join | order | span ;
complement = "complement(" location ")" ;
join       = "join(" item { "," item } ")" ;
order      = "order(" item { "," item } ")" ;
item       = [ remote ":" ] location ;          (* remote = accession.version、namespace は外側と同じ *)
span       = range | position | between | oneof ;
range      = [ "<" ] pos ".." [ ">" ] pos ;
position   = pos ;
between    = INT "^" INT ;                      (* 隣接する2残基の間 *)
oneof      = INT "." INT ;                      (* 範囲内のいずれか1残基 *)
pos        = INT [ codon ] ;
codon      = "c" ( "1" | "2" | "3" ) ;          (* 拡張: タンパク質参照のみ。§3.4 *)
INT        = 1以上の整数 ;
```

- 負の座標、0、空白は使えない。
- 環状配列で原点をまたぐ場合は `join(4000..4641652,1..100)` と書く。

### 3.3 正規形

同じ区間集合は1つの文字列にそろえる。ID として比較・保存するときは正規形を使う。

1. 空白を除く。`n..n` は `n` と書く。
2. 外側と同じ参照を指すリモート参照は省く。
3. join の要素がすべて complement で、参照が同じ場合は、`complement(join(…))` の形にまとめる（要素の順序は INSDC の意味論に従って反転する）。それ以外は要素ごとに complement を付ける。
4. 隣接する区間は**結合しない**（エキソン境界の情報を残すため）。等価かどうかの判定は、別の関数（区間集合としての比較）で行う。
5. コドン拡張が1〜3すべてを覆う場合（`12c1..12c3` など）は拡張を省く。

### 3.4 コドン内の位置の拡張（任意）

- nt 由来の位置を aa の ID として外に出すと、コドン内の位置が失われる。これを保持したい場合に使う。
- 構文は `<残基番号>c<1|2|3>` とする。タンパク質を参照する ID に限って使え、complement とは組み合わせられない。

```
uniprot:Q9BYF1-1:60c2          60番残基のコドンの2文字目
uniprot:Q9BYF1-1:60c2..63c1    範囲
uniprot:Q9BYF1-1:60c2..63      省略した場合、始点は c1、終点は c3 とみなす
uniprot:Q9BYF1-1:60c1^60c2     コドン内の塩基と塩基の間
```

- 正規形では、省略できる `c1`（始点）と `c3`（終点）は省く。
- nt 単位に換算すると、CDS 上の位置 `3(N−1)+k` に当たる。コドン内の位置は、そのタンパク質をどの CDS がコードしているかによらず一意に決まる。
- aa↔aa の写像（UniProt→PDB、アイソフォーム間、オーソログ間）では、コドン内の位置をそのまま引き継ぐ。そのため、ゲノム→UniProt→PDB→ゲノムと往復しても、元の1塩基に戻る。
- 採用しなかった候補: `60.2`（INSDC の「範囲内のいずれか」と衝突）、`60+2`（HGVS のイントロン内の位置と紛らわしい）、`60:2`（参照の区切りと衝突）、`60#2`（IRI のフラグメントと衝突）。
- **拡張部分を取り除くと、そのまま正しい aa の ID（残基を覆う範囲）になる**ことを保証する。拡張を知らない処理系でも、aa の解像度で正しく読める。
- 運用上は、nt 由来のデータは nt 側の ID（ゲノムや mRNA 上の location）を正本とし、aa の ID は派生した見え方とする（HGVS の `c.`/`g.` と `p.` の関係と同じ）。

### 3.5 IRI（identifiers.org 風）

```
https://<togocoordのドメイン>/refseq:NM_014739.3:join(233..235,5958..6116)
                             └─ identifiers.org でも解決できる部分 ─┘└─ location ─┘
```

- identifiers.org と同じ `<prefix>:<accession>` の形と prefix 名を使い、解決は TogoCoord 自身のドメインで行う。content negotiation で HTML、JSON、JSON-LD を返し分ける。
- **location を取り除くと、identifiers.org でそのまま解決できる短縮形の識別子になる**。
- 将来、identifiers.org / bioregistry に `togocoord` という prefix を登録して、そこから転送する形も検討する（ローカル ID に括弧やコロンが入る登録を受け付けてもらえるかは要確認）。
- エンコード: RFC 3986 では、パス中に `( ) , : .` をそのまま書ける。**パーセントエンコードするのは `<`（%3C）、`>`（%3E）、`^`（%5E）の3つだけ**とする。
- 正規形の ID から IRI を作るので、「IRI が一致すること」と「区間集合が一致すること」が同じ意味になる。
- Turtle では、括弧を含む IRI を prefix の短縮形で書けないので、`<…>` の完全な形で書く。

### 3.6 FALDO JSON-LD

FALDO の正式な語彙に合わせる（詳細は [spec-service.md](spec-service.md) §7）。鎖の向きは位置の型（`faldo:ForwardStrandPosition` / `faldo:ReverseStrandPosition`）で表し、逆鎖では begin を数値の大きい側にする。`join` は `faldo:ListOfRegions`、`order` は `faldo:BagOfRegions` で、要素は `rdf:_n` で並べる。`<` と `>` は `faldo:FuzzyPosition`、`^` は `faldo:InBetweenPosition`、`.` は `faldo:InRangePosition`。（当初この節に書いた `NegativeStrand` は FALDO の語彙ではなかったので訂正した。）

- Location ID と JSON-LD は 1:1 で対応させる。
- コドン拡張は、`faldo:Position` に `faldo:codonPosition`（値は 1..3）を追加する形で FALDO に提案する。FALDO に入らなかった場合は、独自の名前空間で `tgc:codonPosition` として定義する。このプロパティを知らない処理系は、aa の位置としてそのまま正しく読める。現在の実装は `tgc:codonPosition` を出力する。

```turtle
[] a faldo:ExactPosition ; faldo:position 60 ;
   faldo:reference <…/Q9BYF1-1> ; tgc:codonPosition 2 .
```

### 3.7 名前空間ごとの ID 規則

名前空間ごとに、次の4点を決める。

1. **接頭辞**: identifiers.org / bioregistry の prefix に合わせる（決定）。
2. **accession の構文**: 名前空間ごとの正規表現。区切り文字の意味もここで決める。
3. **配列の固定方法**: どの表記で配列を1本に固定するか。version がない場合の扱い。
4. **座標を付ける配列と番号体系**。

| 名前空間 | 入力として受け付ける形 | 出力の正規形 | 座標 |
|---|---|---|---|
| insdc / refseq / ensembl | `acc.ver`、`acc`（version なしは最新版として解決する） | `acc.ver` | その配列の残基番号 |
| uniprot | `P12883`、`P12883-2`、代表アイソフォームを番号付きで書いた `P12883-1` など | バージョンなし。**代表アイソフォームは番号なし**（`P12883`）、それ以外は `P12883-2` | 最新版の配列の残基番号 |
| pdb | `4HHB.A`（chain）。拡張 PDB ID（`pdb_00004hhb.A`）も受け付ける | バージョンなし。`4HHB.A` | その chain の **label_seq_id** |
| refget（GFA segment、ユーザのデータ） | `SQ.<digest>` | 同じ | その配列の残基番号 |

**UniProt**
- バージョンや取得日まで管理して使う利用者は少ないため、バージョンなしで扱い、常に最新版の配列として解決する。
- 代表アイソフォームは、UniProt が代表として指定したものとする（`-1` とは限らない）。入力に番号付き（`P12883-1` など）が来た場合も、代表アイソフォームなら番号なしと同じものとして扱う。内部でどちらの形で持つかは、実装しやすい方でよい。
- 配列が更新されると、同じ ID が指す位置がずれることがある。再現性が必要な利用者向けに、レスポンスのメタデータとして、解決に使った配列のダイジェストと UniParc の ID（UPI）、UniProt のリリースを返す。ID 自体には含めない。

**PDB**
- バージョンなしで扱う。
- PDB の名前空間に限り、**例外として `.` を chain の区切りとして扱う**（`4HHB.A`）。参照の単位は chain とする。
- chain ID は **auth_asym_id**（論文やビューアで表示される chain 名。例えば抗体の H/L 鎖）とする（決定）。RCSB の instance 表記（`4HHB.A`）は label_asym_id を使っており、両者が異なる場合がある。label_asym_id との対応はレジストリで持つ。
- 座標は label_seq_id（1から始まる連番）とする。著者番号（挿入コード付きの `52A` など）は座標に使わず、レスポンスで対応表として返す。
- 観測されている残基と未解像の残基は、chain 上のアノテーションとして持つ。

**共通**
- 「人が読める形（`ns:acc[.ver]:loc`）」と「ダイジェストの形（`refget:SQ…:loc`）」を、レジストリで相互に変換できるようにする。
- 入力には寛容に、出力は正規形で返す。

---

## 4. 座標モデルと写像

### 4.1 中心となる考え方

> **INSDC の location 文字列は、それ自体が「feature 配列（1..L）→ 参照配列」への写像である。**

GBFF の CDS、mRNA、exon の各 feature は、書かれている location がそのまま edge になる。変換表を作るための専用コードは要らない。GFA の path も同じ形で表せる（§7）。

### 4.2 内部表現

```ts
// 内部はすべて 0-based half-open。1-based closed との変換は入出力の境界でだけ行う
type Block   = { src: number; tgt: number; len: number; rev: boolean };
type Mapping = { from: SeqRef; to: SeqRef; blocks: Block[]; unit: Unit; provenance: Provenance };
// 不変条件: blocks は src の昇順。重なりは両側とも許す（ribosomal slippage や、アライメントでの重複に対応するため。
//           また invert で src と tgt が入れ替わるため）。重なりがある場合、map は該当する位置をすべて返す
```

写し方の規則は、`rev = false` のとき `y = tgt + (x - src)`、`rev = true` のとき `y = tgt + len - 1 - (x - src)` の2通りだけ。負の座標や complement 専用の delta は使わない。

例: `complement(join(201..300,401..500))`（長さ200）

| src | tgt | rev |
|---|---|---|
| [0,100) | [400,500) | ✓ |
| [100,200) | [200,300) | ✓ |

### 4.3 演算

| 演算 | 内容 |
|---|---|
| `fromLocation(loc)` / `toLocation(mapping)` | Location ID とブロック列を相互に変換する |
| `map(m, interval)` | 区間を写す。結果は写像できた部分と、できなかった部分（unmapped）に分かれる |
| `invert(m)` | src と tgt を入れ替える。逆方向の変換を別に実装する必要はない |
| `compose(a→b, b→c)` | 区間の交差で a→c を作る。すべての連鎖変換はこれで表す |
| `scale(protein ↔ CDS)` | aa の i 番目（0-based）を CDS の nt 区間 [3i+φ, 3i+φ+3) に写す。φ は `/codon_start` から決まる位相 |

- aa と nt の違いは `scale` の中だけで扱う。それ以外の写像はすべて nt↔nt（または aa↔aa）の単位比1で扱う。
- 連鎖変換では、途中の結果を ID として取り出さずに compose するので、nt の解像度が保たれる。

### 4.4 変換の意味論

| 状況 | 規則 |
|---|---|
| 区間の一部しか写せない | 写せた部分ごとに分割して返す。切り詰められた端には `<`（始点側）/`>`（終点側）を付ける。unmapped の部分は別に返す |
| nt → aa | 区間を覆うコドンに対応する aa の区間に広げる。端のコドン内の位置はメタデータとして返し、要求があればコドン拡張（§3.4）で ID に含める |
| aa → nt | コドンの3塩基すべてに写す。エキソン境界をまたぐコドンは join になる |
| `^`（残基の間） | 長さ0の区間 [k,k) として扱う。両側の残基が変換先で隣り合っている場合だけ写し、そうでなければ unmapped にする |
| `a.b`（範囲内のいずれか） | 範囲として写し、「不確か」のフラグを保つ |
| `order()` | 要素ごとに写し、順序を保つ |
| 終止コドン | CDS には含めるが、タンパク質には写さない |
| 部分的な CDS（`<`/`>`） | 位相を考慮して写す。端の残基は不確かとして扱う |
| 出力が複数の区間になる | 変換元の順序で join にし、§3.3 の正規形で出力する |

---

## 5. データモデル

### 5.1 Sequence registry

| 項目 | 内容 |
|---|---|
| `namespace`, `accession`, `version` | ID の構成要素 |
| `digest` | GA4GH refget の sha512t24u。DB 間で配列が同じかどうかの判定に使う |
| `moltype` | DNA / RNA / protein |
| `length`, `topology` | 長さと、線状か環状か |
| `taxon`, `assembly` | NCBI Taxonomy ID と、所属する assembly |
| `aliases` | UCSC の染色体名（`chr1`）など、別名の一覧 |

### 5.2 Mapping edge

| 項目 | 内容 |
|---|---|
| `from`, `to`, `blocks`, `unit` | §4.2 のとおり |
| `kind` | `annotation` / `identity` / `alignment` / `liftover` / `orthology` / `graph` など。後から種別を追加できる |
| `provenance` | データ源、ファイル、リリース、feature や qualifier、どの段階（T0〜T3）か、取り込み元の種別（配布 / 自前の計算 / ユーザ） |
| `validation` | 取り込み時の自己検証の結果（§6.3） |

### 5.3 Annotation

- 配列上に Location ID で置かれた情報（ドメイン、PTM 部位、CRE、バリアント、二次構造、PDB の観測残基など）を、写像とは別の層として持つ。
- **アノテーションの伝播は、アノテーションの Location ID を `map` するだけで済む**。

### 5.4 保存

assembly 単位で SQLite、DuckDB、Parquet などに保存する。ブロック列は配列型の列で持つ。RDF（FALDO）としての公開は別途行う。

---

## 6. データ源とアダプタ

### 6.1 アダプタのインターフェース

```ts
interface Adapter {
  meta: { id: string; layer: "core" | "enrichment"; taxa: number[] | "any";
          license: string; updateCycle: string; tier: "T0" | "T1" | "T2" | "T3" };
  ingest(input): AsyncIterable<SequenceRecord | Mapping | Annotation>;
}
```

ヒトやマウスのリソースを追加するときは、アダプタを1つ足すだけで済み、コアのコードには手を入れない。

### 6.2 段階とデータの入手元

段階は「どの粒度の対応ができるか」で決める。データをどこから得るかは別の軸として扱う。

| 段階（能力） | 配布データ | 自前で計算 | ユーザが持ち込む |
|---|---|---|---|
| **T0** 注釈と同一配列 | GBFF / GFF3 | ― | GFF3 + FASTA |
| **T1/T2** 配列どうしの対応 | NCBI `cDNA_match`, SIFTS | タンパク質のペアワイズ、スプライスアライメント | PAF など |
| **T3** ゲノム全体の対応 | UCSC chain, Ensembl Compara, GRC alignments, HAL | minimap2, wfmash | GFA, chain, PAF, MAF |

段階ごとにできることは次のとおり。

- **T0**（全生物種）
  - 同じ assembly の中での genome ↔ mRNA ↔ CDS ↔ protein ↔ exon の変換（注釈されたアイソフォームや alt. ORF を含む）
  - 配列が同一の DB 間（GenBank protein、RefSeq、UniProt、Ensembl、AlphaFold DB）の対応
  - 配列が同一の GCA/GCF 間の対応
- **T1**（主にモデル生物）
  - RefSeq 転写産物とゲノムの不一致の吸収
  - UniProt と PDB の残基単位の対応
  - assembly 間の liftover、ヒトとマウスの非コード領域を含む種間の変換
- **T2**（全生物種、コストは小さい）
  - 配列がわずかに違う DB 間の対応（MANE と UniProt の不一致など）
  - アイソフォーム間の対応
  - オーソログ対のタンパク質アライメントによる、**コード領域に限った種間の変換**
  - ゲノムに注釈されていない mRNA をゲノムに置くこと

  いずれもオンデマンドで計算し、キャッシュする。
- **T3**
  - **非コード領域を含む種間の変換**と、非モデル生物での assembly 間の変換
  - まず配布データとユーザの持ち込みで対応する。自前で計算するかは、需要を測ってから判断する

### 6.3 取り込み時の自己検証

- CDS の配列を取り出して翻訳し、`/translation` やタンパク質配列と照合する。不一致のものはフラグを付け、alignment edge に切り替える。
- 扱う必要がある特殊ケース:
  - `/codon_start`
  - `/transl_except`
  - ribosomal slippage
  - RNA editing
  - 部分的な CDS
  - `/exception`
- identity edge は、ダイジェストが一致した場合にだけ作る。

### 6.4 Core と Enrichment の対応表（ポスターの図の要素）

| 図の要素 | 層 | 表現方法 |
|---|---|---|
| Genome / mRNA / CDS / Protein / Exon | Core | feature の location をそのまま写像にする |
| Protein isoform, alt. Splicing | Core（UniProt isoform は Enrichment） | 転写産物・タンパク質ごとに別の配列ノードを持つ |
| alt. TSS | Enrichment（FANTOM CAGE など） | アノテーションとして持つ。新しい転写産物モデルを作る場合は配列ノードとして追加する |
| alt. ORF | Core / Enrichment | 1つの mRNA に複数の CDS 写像を持たせる |
| Structure, Unresolved residues | Enrichment（SIFTS, mmCIF, AlphaFold DB） | chain（`4HHB.A`）をノードにし、座標は label_seq_id とする。観測残基は chain 上のアノテーションとして持つ |
| Domain, α helix / β sheet | Enrichment | タンパク質上のアノテーションとして持つ |
| Genome of other organism | Enrichment（chain, Compara） | `liftover` / `orthology` の edge |
| Pangenome graph | Enrichment / ユーザの持ち込み | `graph` edge として持つ（§7） |
| CREs, Variant | Enrichment（ChIP-Atlas, TogoVar） | ゲノム上のアノテーションとして持ち、写像で伝播させる |

---

## 7. GFA とユーザによるデータの持ち込み

### 7.1 GFA の扱い

GFA の path は、INSDC location のリモート参照と complement でそのまま書ける。

```
P  hapA  s1+,s2-,s3+  *
→ hapA = join(s1:1..L1, complement(s2:1..L2), s3:1..L3)
```

- segment を、実在する配列（ダイジェスト付き）としてレジストリに登録する。
- path（W 行も同様）は、「path 配列 → segment 群」という写像 edge にする。
- 2本の path のあいだの対応は `compose(pathA→segs, invert(pathB→segs))` で求める。新しい演算は要らない。
- 片方の path にしか現れない segment は unmapped になる。変換先では `^` で位置を示す。
- L 行にオーバーラップがあるグラフは、取り込み時に切り詰めて blunt に変換する。
- path 名は PanSN 命名（`sample#hap#contig`）を読み取る。さらに、ダイジェストの一致で既知の配列に結び付け、公開データの注釈を伝播できるようにする。
- 大規模なグラフでは、全 path の組み合わせを事前に計算しない。問い合わせのあった組だけを、その都度 compose する。rGFA のタグ（SN/SO/SR）があれば、参照配列上の座標への近道として使う。

### 7.2 ワークスペース

- 持ち込まれたデータはワークスペース単位で分離する。既定で非公開とし、由来の種別は `user` として記録する。
- ユーザ独自の配列は `refget:` ダイジェストで参照する。同じ配列を別の場所で読み込んでも、同じ ID になる。
- 取り込み時に、ダイジェストの一致で検証する。
- 取り込めるファイルの大きさや、問い合わせあたりの compose の規模に上限を設ける。上限を超える規模のデータは、配布する CLI で扱ってもらう。

---

## 8. 経路探索

- 配列を節点、edge を辺とするグラフ上で、重み付きの最短経路を探す。重みは種別ごとに決め、概ね identity < annotation < alignment（配布） < alignment（自前の計算） < liftover/orthology の順とする。
- **taxon ごとのプロファイル**で、優先順位を設定として上書きできるようにする。例えばヒトでは MANE Select を優先する。設定がない種では、既定の方針（一次リポジトリの annotation を優先）で動く。
- 経路は1本だけでなく、候補を複数返せるようにする（例: 対応する PDB 構造がすべて欲しい場合）。
- **変換経路が見つからなかった問い合わせ**（対象の種、変換元と変換先の層）を記録し、T3 への投資判断の材料にする。

---

## 9. API（案）

| メソッド | パス | 内容 |
|---|---|---|
| GET/POST | `/v1/convert` | `loc`（Location ID）、`to`（名前空間・配列・レイヤー）、`via`、`profile`、`codon`。POST ではバッチ処理 |
| GET | `/v1/location/parse` | 正規形と、ブロック列の形の JSON を返す |
| GET | `/v1/location/faldo` | FALDO JSON-LD を返す |
| GET | `/v1/sequences/{ref}` | レジストリの情報 |
| GET | `/v1/sequences/{ref}/neighbors` | 直接つながる edge の一覧 |
| GET | `/v1/annotations` | 指定した区間にかかるアノテーションを、変換先の座標に伝播して返す |
| POST | `/v1/workspaces/{id}/datasets` | ユーザがデータを持ち込む |

レスポンスの例:

```json
{
  "input": "uniprot:Q9BYF1-1:60",
  "results": [{
    "location": "insdc:NC_000020.11:complement(join(…))",
    "pieces": [{ "source": "60", "target": "…", "codon": {"begin": 1, "end": 3} }],
    "unmapped": [],
    "path": [{ "from": "…", "to": "…", "kind": "identity", "provenance": {…} }],
    "faldo": {…}
  }]
}
```

---

## 10. テスト方針

1. **素朴な実装を正解（オラクル）にする**: すべての写像を残基ごとの配列に展開して計算する実装を別に用意し、ブロック演算の結果と突き合わせる（ランダム入力による property-based test）。
2. **演算の性質を検証する**:
   - `invert(invert(m)) = m`
   - `mapped(x) ⊆ map(invert(m), map(m, x))`。m が単射なら等号が成り立つ
   - compose が結合的であること
   - パースと出力を往復させると正規形が一致すること
3. **難しいケースを集めたテストデータ**:
   - 逆鎖
   - エキソン境界で分断されるコドン
   - codon_start が 2 または 3
   - ribosomal slippage
   - transl_except
   - 部分的な CDS
   - 環状ゲノムで原点をまたぐもの
   - RefSeq とゲノムの不一致（cDNA_match）
   - PDB の挿入コードと未解像残基
   - GFA の逆向き segment と bubble
4. **取り込み時の自己検証**（§6.3）を、データに対する回帰テストとしても使う。

---

## 11. ロードマップ

| フェーズ | 内容 |
|---|---|
| 0 | 本仕様を確定する。テストデータを整備する。**完了（2026-09-18）**: [spec-core.md](spec-core.md)、`core/test/corpus/` |
| 1 | コアライブラリ（パーサ、正規化、ブロック演算、意味論）とオラクルテスト。**v0.1 完了（2026-09-18）**: `core/`（テスト104件） |
| 2 | GBFF と GFF3 のアダプタ。汎用性を確かめるため、ヒト・マウスに加えて、性質の異なる種（環状ゲノムの細菌、シロイヌナズナなど）を最初から対象に含める。**v0.1 完了（2026-09-18）**: `ingest/`、[spec-ingest.md](spec-ingest.md)。ウイルス、ヒトのミトコンドリアゲノム、アデノウイルス、プラスミド、ヒト GRCh38 の cDNA_match の実データで検証。ヒトやマウスの全ゲノム規模の GFF3、シロイヌナズナ、Ensembl の seqid は未検証 |
| 3 | REST API、経路探索、Web UI。T1 の Enrichment（SIFTS、cDNA_match、UCSC chain、MANE）。**3a・3b 完了（2026-09-18）**: `service/`（経路探索、複数の保存先）、[spec-service.md](spec-service.md)。**3c 完了**: REST API、FALDO JSON-LD、同一配列（ダイジェスト）による経路。**3e（SIFTS）完了**: Ensembl・UniProt・SIFTS の取り込みと、構造との往復の検証（spec-service §8）。UI と、SIFTS 以外の拡張データ（UCSC chain、MANE）は未着手 |
| 4 | T2（オンデマンドのアライメントとキャッシュ）、アノテーションの伝播（ポスターのユースケースの再現） |
| 5 | T3（GFA と chain の持ち込み、ワークスペース）。需要に応じて、自前の計算も検討する |

---

## 12. 未決事項

- [ ] FALDO 開発者への `faldo:codonPosition` の提案（§3.6）
- [ ] TogoCoord のドメイン名。identifiers.org への prefix 登録の可否（§3.5）。IRI の形は identifiers.org 風で決定
- [ ] 運用する場所（DBCLS / DDBJ）と、計算資源の規模
- [ ] オーソログ情報の入手元（OrthoDB、eggNOG、Ensembl Compara）
- [ ] 経路が見つからなかった問い合わせのログの扱い（プライバシー）
