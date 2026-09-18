# TogoCoord コア仕様（v0.1.1）

2026-09-18。[design.md](design.md) の §3（Location ID）と §4（座標モデル）を、実装できる粒度まで厳密にしたもの。コアライブラリ（`core/`）はこの仕様に従い、テストデータ（`core/test/corpus/`）で検証する。

---

## 1. 範囲

- 扱う: Location ID のパース、意味値への変換、正規形での出力。写像（ブロック列）の構築・反転・合成。写像による location の変換。
- 扱わない: 配列の取得、データ源の読み込み（アダプタ）、経路探索、ID の解決（最新版の version の補完、UniProt の代表アイソフォームの判定など）。これらはサービス層で扱う。コアは、必要な情報（配列の単位など）を呼び出し側から受け取る。

---

## 2. 内部座標

- すべて **0-based half-open** の区間 `[start, end)` で扱う。1-based closed の表記との変換は、パースと出力のときにだけ行う。
- 座標の単位（unit）は配列の種類で決まる。
  - 塩基配列（`nt`）: 1塩基が1単位。
  - アミノ酸配列（`aa`）: **1残基を3単位**（コドン単位）とする。残基 n（1-based）のコドン内の k 文字目（1..3）は、単位 `3(n−1)+(k−1)` に当たる。
- この規約により、aa↔nt の変換はコドン単位と塩基の1対1の対応として、aa↔aa の変換は3単位ずつの対応として扱える。そのため、**コドン内の位置は aa↔aa の写像を通しても保たれる**。
- 配列の単位は、呼び出し側が関数 `unitOf(ref)` で与える。与えない場合は、名前空間の既定値を使う（§3.4）。

---

## 3. Location ID

### 3.1 全体

```
<namespace>:<accession>:<location>
```

- 最初の `:` までが namespace、次の `:` までが accession、残りが location。accession に `:` は含まれない。
- 配列の内部キーは `<namespace>:<accession>`（例: `refseq:NC_045512.2`、`pdb:4HHB.A`）とする。
- 空白は、パースの前にすべて取り除く。

### 3.2 location の文法

```ebnf
location   = [ remote ":" ] body ;
body       = "complement(" location ")"
           | "join(" location { "," location } ")"
           | "order(" location { "," location } ")"
           | span ;
span       = range | point | between | oneof ;
range      = [ "<" ] pos ".." [ ">" ] pos ;
point      = pos ;
between    = pos "^" pos ;
oneof      = INT "." INT ;
pos        = INT [ "c" ( "1" | "2" | "3" ) ] ;
INT        = 1以上の10進整数 ;
remote     = 英字を1文字以上含む [A-Za-z0-9_.-]+ ;    (* 外側と同じ名前空間の accession *)
```

**パースエラーとする条件**

| 条件 | 例 |
|---|---|
| 0 や負の数 | `0`、`-5` |
| 範囲の始点が終点より大きい | `12..5`（原点をまたぐ場合は join で書く） |
| `<` が終点側に、`>` が始点側にある | `1..<5`、`>1..5` |
| 1点の位置に `<` や `>` が付いている | `<5` |
| コドン拡張が 1..3 以外 | `5c4` |
| 残基の間（`^`）で、両側が隣り合っていない | `5^7` |
| `oneof` にコドン拡張や `<`、`>` が付いている | `5c1.7` |
| `order` が入れ子になっている | `join(order(1..2),3..4)` |

### 3.3 意味値（セグメント列）への変換

パースした構文木は、**向き付きのセグメントを並べた列**に変換する。

```ts
Segment = { ref, start, end, strand: +1 | -1, fuzzyLow?, fuzzyHigh?, uncertain? }
Location = { outer: ref, kind: "join" | "order", segments: Segment[] }
```

- **並び順（traversal order）**: セグメントは、feature を 5'→3'（N末→C末）にたどる順に並べる。
- `complement(X)`: X のセグメント列を逆順にし、各セグメントの strand を反転する。`complement(join(A,B))` は `[B⁻, A⁻]` になり、`join(complement(B),complement(A))` と同じ値になる。
- `join` と `order`: 子のセグメント列を連結する。入れ子の join は平坦にする。最上位に `order` がある場合だけ、kind を `order` とする。
- **fuzzy は数値の上での端として持つ**。`<` は数値の小さい側（`fuzzyLow`）、`>` は大きい側（`fuzzyHigh`）に付く。INSDC と同じく、生物学的な 5'/3' ではなく数値の大小で区別する。
- `remote:` は、その部分のセグメントの ref を `<外側の名前空間>:<remote>` にする。
- **単位の変換**（aa の場合）:
  - 位置 `n` → `[3(n−1), 3n)`
  - 位置 `nck` → `[3(n−1)+k−1, 3(n−1)+k)`
  - 範囲の始点 `n[ck]` → `3(n−1) + (k−1、省略時は 0)`
  - 範囲の終点 `n[ck]` → `3(n−1) + (k、省略時は 3)`（排他的な終端）
- **残基の間 `a^b`**: 長さ0のセグメント `[k, k)` にする。k は右側の単位の位置。
  - nt: `b = a+1` が必要で、k = a。
  - aa: 左の単位を `L = 3(a−1)+(x−1、省略時は 2)`、右の単位を `R = 3(b−1)+(y−1、省略時は 0)` とし、`R = L+1` が必要で、k = R。したがって `60^61`、`60c3^61c1`、`60c1^60c2` はすべて有効。
- `oneof a.b`: `[a の始点, b の終点)` とし、`uncertain` を立てる。
- **意味上のエラー**:
  - aa の配列に strand −1 がある（タンパク質に complement は使えない）
  - nt の配列にコドン拡張がある
  - 長さがわかっている配列で、区間が配列の範囲を超える

### 3.4 名前空間

| prefix | accession のパターン（概略） | 既定の単位 |
|---|---|---|
| `insdc` | `[A-Z]{1,6}\d{5,}(\.\d+)?` | 英字3文字＋数字（例: `AAF99721`）は aa、それ以外は nt |
| `refseq` | `[A-Z]{2}_[A-Z0-9]+(\.\d+)?` | `NP_`・`XP_`・`YP_`・`WP_`・`AP_` は aa、それ以外は nt |
| `ensembl` | `ENS[A-Z]*[EGTP]\d{11}(\.\d+)?` | `…P` は aa、それ以外は nt |
| `uniprot` | UniProt の accession ＋ `(-\d+)?` | aa |
| `uniparc` | `UPI[0-9A-F]{10}` | aa |
| `pdb` | `[0-9][A-Za-z0-9]{3}\.<chain>` または `pdb_\d{4}[0-9][A-Za-z0-9]{3}\.<chain>` | aa |
| `refget` | `SQ\.[A-Za-z0-9_-]{32}` | 既定値なし（呼び出し側が与える） |

- この表のパターンは、ID の形としての最低限の検査である。名前空間は、登録用の関数で追加できる（テスト用の `test` もこれで登録する）。
- 正規化として、prefix は小文字にする。PDB の ID は、4文字形式なら大文字にし、拡張形式（`pdb_0000xxxx`）で4文字形式に直せるものは4文字形式にする。chain は大文字と小文字を区別する（auth_asym_id）。
- **要確認**: prefix 名は bioregistry と照合する必要がある。例えば bioregistry には `uniprot.isoform` が別の prefix として存在する。現時点では `uniprot` の下でアイソフォームも受け付ける。

### 3.5 正規形での出力

意味値から文字列を作る規則。**意味値が同じなら、出力される文字列も必ず同じになる**。

1. セグメントが1つなら、そのセグメントを書く（strand が −1 なら `complement(…)` で囲む）。
2. セグメントが2つ以上で、**すべてが strand −1 かつ同じ ref** の場合は、`complement(join(…))` の形にする。中の要素は、セグメント列を逆順にし、strand を + にして書く。
3. それ以外は `join(…)`（kind が order なら `order(…)`）とし、要素ごとに必要なら `complement` を付ける。
4. ref が外側と異なるセグメントには `accession:` を前に付ける。
5. 隣接している区間は結合しない。
6. セグメントの書き方:
   - nt: 長さ1で fuzzy がなければ `n`、そうでなければ `[<]a..[>]b`。長さ0は `k^k+1`。uncertain なら `a.b`。
   - aa: 始点を `残基[c k]`（k=1 なら省く）、終点を `残基[c k]`（k=3 なら省く）とする。1残基をちょうど覆う場合は `n`、1単位だけなら `nck`。
   - aa の長さ0（位置 k）: k が3の倍数なら `k/3 ^ k/3+1`、そうでなければ `LcX^RcY`（コドン内の位置）。
7. **コドン拡張を省く出力モード**（`codon: "never"`）: `c` を取り除き、残基全体を覆う形にする。コドン内の `^` は、その残基の位置 `n` とする。

---

## 4. 写像

### 4.1 ブロック

```ts
Block = { srcRef, src, tgtRef, tgt, len, rev }
```

- src 側の単位 `x ∈ [src, src+len)` を、次の位置に写す。
  - `rev = false`: `tgt + (x − src)`
  - `rev = true`: `tgt + len − 1 − (x − src)`
- 不変条件: `len > 0`、`src ≥ 0`、`tgt ≥ 0`。**src 側でも tgt 側でも重なりを許す**（ribosomal slippage、アライメントの重複、および invert に対して閉じているため）。
- 写像（Mapping）はブロックの集合である。ある単位の写像先は0個以上ある。

### 4.2 演算

| 演算 | 定義 |
|---|---|
| `invert(m)` | 各ブロックの src と tgt を入れ替える（`rev` はそのまま） |
| `compose(ab, bc)` | ab の各ブロック a と bc の各ブロック b について、a の tgt 区間と b の src 区間の交差 `[lo, hi)` があれば、ブロックを1つ作る。src 側の始点は、`a.rev` なら `a.src + (a.tgt + a.len − hi)`、そうでなければ `a.src + (lo − a.tgt)`。tgt 側の始点は、`b.rev` なら `b.tgt + (b.src + b.len − hi)`、そうでなければ `b.tgt + (lo − b.src)`。長さは `hi − lo`、`rev` は `a.rev ≠ b.rev` |
| `fromLocation(F, loc)` | feature 配列 F から参照配列への写像。セグメントを並び順にたどり、累積の位置 `off` について `{F, off, seg.ref, seg.start, len, seg.strand = −1}` を作る |
| `cdsMapping(P, cds, codonStart, aaLength)` | `compose(scale, fromLocation(P#cds, cds))`。scale は `{P, 0, P#cds, codonStart−1, 3·aaLength, false}` |

- **先頭の不完全なコドン**（`leadingPartialCodon`、v0.1.2）: Ensembl は、5' 側が欠けた CDS のタンパク質の先頭に、欠けたコドンを表す `X` を1残基置く（INSDC の `/codon_start` にはこの残基がない）。この選択肢を指定すると、1番目の残基の最後の `codonStart − 1` 単位を、CDS の先頭の塩基に対応させる。scale は `{P, 3−(codonStart−1), P#cds, 0, 3·aaLength−(3−(codonStart−1)), false}`。1番目の残基をゲノムに変換すると、欠けた部分があるので始点に切り詰めの印が付く（例: `<930312..930313`）。
- `aaLength` は、呼び出し側が必ず与える（コアでは推定しない）。終止コドンの有無、不完全な終止コドン（`transl_except` の TERM）、3' 側が部分的な CDS などを、コアでは判断できないため。
- `3·aaLength + codonStart − 1` が CDS の長さを超える場合はエラーとする。

---

## 5. 変換（`mapLocation`）

入力の Location の各セグメント s について、並び順に次の処理を行う。

1. **断片（piece）**: s と交差するすべてのブロックについて、交差する部分を写す。写像先の strand は `s.strand × (rev ? −1 : +1)`。
2. **断片の順序**: s の並び順に従う（+ なら src の始点の昇順、− なら src の終点の降順）。同じ位置の場合は、ブロックの登録順とする。
3. **写像できなかった部分（unmapped）**: s のうち、どの断片にも覆われない部分を、s と同じ strand で、並び順に返す。
4. **結合**: 並び順で連続する断片 p、q（入力のセグメントをまたいでもよい）について、次の条件をすべて満たせば1つにまとめる。
   - tgt の ref と strand が同じ
   - tgt が並び順で連続している（+ なら `q.start = p.end`、− なら `q.end = p.start`）
   - 同じセグメント内なら、src が並び順で後ろに進んでいる（重なっていない）
   - 別のセグメントなら、入力の kind が `join` で、どちらのセグメントも uncertain（`a.b`）でない

   イントロンをまたぐゲノム区間を写すと、タンパク質側では連続した1区間になる。エキソンごとの join で書かれた CDS の location や、エキソン境界で分断されたコドンも、タンパク質側では1区間になる。一方、slippage で同じ塩基が2つの位置に写る場合は、結合しない（v0.1.1 で、セグメントをまたぐ結合を追加した）。
5. **fuzzy**: まとめた断片の並び順での始端（最初の断片の、そのセグメントでの始端）と終端（最後の断片の、そのセグメントでの終端）ごとに、「切り詰めがあったか」を判定する。
   - そのセグメントの端と一致する場合は、入力の fuzzy を引き継ぐ（+ なら始端が `fuzzyLow`、− なら始端が `fuzzyHigh`）。
   - セグメントの内側にある場合は、並び順ですぐ外側の単位が、**同じ変換先の配列への**断片に覆われていなければ、切り詰めありとする。ほかの配列に写っているかどうかは関係しない。そのため、ある配列への結果が、無関係な edge の有無に左右されない（v0.1.1 で変更）。
   - 切り詰めを tgt の数値の端に変換する（tgt の strand が + なら始端を Low に、− なら始端を High に）。
6. **残基の間（長さ0のセグメント `[k,k)`）**: 単位 k−1 と k をそれぞれ写す。**変換先の配列ごとに**、次の条件をすべて満たせば、その配列上の `[max, max)` に写す。
   - その配列の中で、どちらも写像先がちょうど1つである
   - 写像先が隣り合っている（位置の差の絶対値が1）

   どの配列でも条件を満たさない場合は unmapped とする。配列ごとに判定するのは、読み枠の異なる重なり合った遺伝子（例: ヒトのミトコンドリアゲノムの ATP8 と ATP6）があっても、それぞれの配列では位置が一意に決まるためである（v0.1.1 で変更。v0.1 では全体でちょうど1つを条件にしていた）。
7. **uncertain**: s が uncertain で、まとめた結果が1区間になった場合は、その区間に uncertain を引き継ぐ。複数の区間になった場合は、結果全体に `uncertain: true` を立てる。

**出力**
- まとめた断片を、tgt の ref ごとにまとめる（最初に現れた順）。ref ごとに1つの Location（kind は入力と同じ）を作る。
- **tgt が aa の場合**: strand がすべて −1 なら、並びを逆にして strand を + にし、`orientation: "reverse"` とする。strand が混在する場合は、すべて + にして `"mixed"` とする。数値の上での fuzzy は、そのまま保つ。
- 結果として、断片の一覧（表示用、まとめる前のもの）、写像先の Location の一覧、unmapped の Location を返す。

---

## 6. v0.1 で扱わないもの

- 環状配列での `n^1`（末尾と先頭の間）。パースは通るが、変換はエラーとする。
- 配列の長さの検査。長さが与えられた場合だけ行う。
- ID の解決（version の補完、アイソフォームの正規化、PDB の chain の対応付け）。これはサービス層で扱う。

---

## 7. テストデータ（`core/test/corpus/`）

| ファイル | 由来 | 検証すること |
|---|---|---|
| `sars2_orf1ab_slippage.json` | NCBI の NC_045512.2（2026-09-18 に取得）。CDS `join(266..13468,13468..21555)`、YP_009724389.1（7096 aa） | slippage（1塩基が2つのコドンに属する）、終止コドンが unmapped になること、切り詰めの fuzzy |
| `human_mt_nd6_minus.json` | NC_012920.1 の ND6 `complement(14149..14673)`、YP_003024037.1（174 aa） | マイナス鎖、aa を逆向きに出力すること、fuzzy の向きの変換 |
| `human_mt_nd1_partial_stop.json` | NC_012920.1 の ND1 `3307..4262`、`transl_except=(pos:4261..4262,aa:TERM)`、318 aa | 不完全な終止コドン（aaLength を明示すること） |
| `human_mt_dloop_circular.json` | NC_012920.1 の D-loop `complement(join(16024..16569,1..576))` | 原点をまたぐこと、complement(join) の正規形、原点をまたぐ入力が feature 上で1区間にまとまること |
| `synthetic_plus_split_codon.json` | 合成データ | + 鎖で、コドンがエキソン境界で分断される場合、イントロンをまたぐ区間をまとめること、`^` |
| `synthetic_minus_split_codon.json` | 合成データ | − 鎖で、コドンが分断される場合 |
| `synthetic_aa_to_aa.json` | 合成データ（UniProt→PDB を想定） | aa↔aa でコドン内の位置が保たれること、未解像の領域、aa の `^` |

**フェーズ2（アダプタ）で実データを確認する候補**
- セレノプロテイン（`transl_except` で Sec を指定しているもの）
- RefSeq の転写産物とゲノムの不一致（NCBI GFF の `cDNA_match` と Gap）
- `codon_start` が 2 または 3 の部分的な CDS（INSDC）
- PDB の挿入コード（Kabat 番号の抗体の構造）と未解像の残基
- 細菌の環状ゲノムで、原点をまたぐ CDS
- GFA（逆向きの segment、bubble）
