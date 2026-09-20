# TogoCoord アダプタ仕様（v0.3：GBFF / GFF3 / FASTA / SIFTS）

[English](spec-ingest.md) | 日本語

2026-09-18。[design.md](design.ja.md) の §6 のうち、フェーズ2で実装した部分の規則。実装は `ingest/`。

---

## 1. 出力

アダプタは次の3種類を出力する（`ingest/src/model.ts`）。

| 種類 | 内容 |
|---|---|
| `SequenceRecord` | 参照キー、分子種、単位、長さ、トポロジー、taxon、refget ダイジェスト（配列が手元にある場合） |
| `Edge` | `from → to` のブロック列（コアの単位）、表示用の Location ID、属性、由来、**自己検証の結果** |
| `Annotation` | feature を Location ID で置いたもの（種類と、qualifier や attribute） |

CLI（`togocoord-ingest`）は、これらを JSON Lines（`{"record": "sequence" | "edge" | "annotation" | "warning", …}`）で標準出力に書き、検証結果の集計を標準エラーに書く。

---

## 2. 参照キー

- Ensembl の安定 ID（`ENST…`、`ENSP…`、他種の `ENSMUSP…` など）は `ensembl:`、accession に `XX_` の形の接頭辞があれば `refseq:`、それ以外は `insdc:` とする（例: `refseq:NC_012920.1`、`insdc:AAF99721.1`、`ensembl:ENSP00000407375.1`）。
- GFF3 の seqid が accession でない場合（Ensembl の `1` など）は、`seqidToRef` で対応を与える。CLI では `--seqid-map` に NCBI の assembly report を渡すと、配列名（`1`、`MT`）、GenBank の accession（`CM000663.2`、`KI270728.1`）、UCSC の名前（`chr1`）を、すべて RefSeq の accession に対応させる。対応がない場合は警告を出して、その seqid の feature を読み飛ばす。
- Ensembl の GFF3 は version を別の属性（`transcript_id=ENST00000419783;version=3`、CDS 行の `protein_id=…;version=1`）に書くので、結合して `ensembl:ENST00000419783.3` とする。version がない ID は、`--fasta` の配列の ID から補う（`VersionResolver`）。

---

## 3. GenBank / GenPept（`ingestGenBank`）

| feature | 処理 |
|---|---|
| `source` | taxon と organism を SequenceRecord に入れる |
| `CDS`（核酸レコード、`/pseudo` なし） | `/protein_id` のタンパク質 → レコードへの edge（`cdsMapping`）。`/codon_start` と `/transl_table` を使う。aaLength は `/translation` の長さ |
| `CDS`（GenPept） | `/coded_by` の location を CDS として、レコード自身（タンパク質）→ 核酸配列への edge。aaLength はレコードの長さ |
| RNA 系（mRNA、ncRNA など）で `/transcript_id` があるもの | 転写産物 → レコードへの edge（`mappingFromLocation`） |
| すべて（`source` を除く） | Annotation（`/translation` は除く） |

- location をパースできない feature は、警告を出して読み飛ばす。例: `bond()` や、v0.1 で扱わない `n^1`。

---

## 4. GFF3（`ingestGff3`）

- **行のまとめ方**: 同じ `ID` を持つ行を1つの feature とする。
- **区間の並べ方**: 行を開始位置で並べる。`-` 鎖の場合は、それを逆順にする。
- **環状配列**（`region` の `Is_circular=true`）: 配列長を超える行（例: D-loop の `16024..17145`、配列長は 16569 nt）は原点で折り返し、2区間に分ける。線状の配列で配列長を超える行は、警告を出して読み飛ばす。
- **部分的な feature**: `start_range=.,N` は `<`（数値の小さい側）、`end_range=N,.` は `>`（大きい側）とする。
- **CDS**:
  - `protein_id` のタンパク質 → seqid の配列への edge にする。`protein_id` がない場合、偽遺伝子（`pseudo=true`）は警告なしで読み飛ばす。それ以外（免疫グロブリンの遺伝子断片など）は警告を出す。
  - codon_start は、5' 側の端の行の `phase + 1` とする。
  - GFF3 には翻訳配列がないので、aaLength は §5 の規則で推定する。
- **転写産物**: RNA 系の型、または `transcript_id` を持つ feature（Ensembl の gene segment など）。ただし、exon・CDS・UTR など転写産物の部分を表す型は、NCBI では `transcript_id` を持っていても転写産物として扱わない。
- **RNA 系**（mRNA、ncRNA、lnc_RNA など）: `Parent` で紐づく **exon の行から location を組み立てる**。NCBI の GFF3 では、mRNA 自体は遺伝子の範囲全体を表す1行にすぎないため。エキソンは親より後に来るので、転写産物は保留しておく。一定数（1000件）の feature が流れても更新されないとき、別の配列に移ったとき、または最後に確定させる。エキソンがない場合は、転写産物自身の行を使う。`transcript_id` があれば、転写産物 → 配列への edge にする。
- **アライメント**（`cDNA_match`、`match` など、`Target` を持つ行）: `Target` の配列 → seqid の配列への alignment edge にする。
  - 自己検証: ブロックごとに同一塩基数を数え、**自分で数えた一致率が5割以上なら ok** とする（Gap の読み違いで座標がずれると25%前後に落ちる）。NCBI の `num_mismatch` と数が違う場合は、detail に `(counted N)` と記録する（例: NM_001291281.3 は、数えると3、`num_mismatch=4`）。v0.2 では報告値との一致で判定していたため、座標は正しいのに不一致とした例がヒトで23件あった。
  - Gap（CIGAR）の各操作の意味:

    | 操作 | 意味 |
    |---|---|
    | `M` | 対応する（ブロックを作る） |
    | `I` | Target 側だけが進む |
    | `D` | ゲノム側だけが進む |

  - **Gap の操作は、行の鎖の向きにゲノムをたどる順に書かれている**（+ 鎖の行ではゲノムの昇順、− 鎖の行では降順）。Target がどちらの鎖でも同じ規則で読む。GFF3 の仕様からは読み取れないため、NCBI GRCh38 の実データで確認した。
    | 行の鎖 / Target の鎖 | 確認に使ったデータ | この規則での一致 | 逆の順序での一致 |
    |---|---|---|---|
    | − / + | NM_012234.7（`Gap=M536 I1 M1191 I3 M2314`） | 4036/4041（不一致5 = `num_mismatch=5`） | 1977/4041 |
    | + / − | NG_162589.1（`M217 D5 M53`）、NG_043318.1（`M216 I3 M85`、2か所） | 270/270、301/301 | 150/270、204/301 |
    | + / −（Gap なし） | NG_044083.1 | 304/304 | — |
  - v0.2 で扱わないもの: `F`/`R`（フレームシフト）。警告を出す。
- **両方の鎖にまたがる feature**（トランススプライシング。例: 葉緑体やミトコンドリアの rps12）: 位置では並べられないので、ファイルの行の順をそのまま使い、行ごとの鎖を保つ。NCBI の GFF3 は転写産物の順に書いていることを、シロイヌナズナの NP_051038.1（`join(complement(69611..69724),139856..140087,140625..140650)`）の翻訳が公開配列と一致することで確認した。

---

## 5. タンパク質長の推定（翻訳配列がない場合）

`coding = CDS の長さ − (codon_start − 1)` とする。

| 条件 | aaLength |
|---|---|
| 3' 側が部分的（`>`） | `floor(coding / 3)` |
| 3' 側が完全で、`coding` が3の倍数 | `coding / 3 − 1`（最後のコドンは終止コドン） |
| 3' 側が完全で、端数がある | `floor(coding / 3)`（端数は、ポリA 付加で完成する不完全な終止コドンとみなす） |

この規則による推定値は、テストに使った実データ（SARS-CoV-2、ヒトのミトコンドリアゲノム、アデノウイルス）の63件の CDS すべてで `/translation` の長さと一致した。

---

## 6. 取り込み時の自己検証

| 種類 | 方法 | 結果 |
|---|---|---|
| CDS | CDS の配列を取り出して翻訳し、比較する | 比較する相手: `/translation`（GBFF）、タンパク質レコードの配列（GenPept）、`--fasta` で与えた公開タンパク質配列（GFF3）。どれもない場合は、途中に終止コドンがないことを確認する。NCBI の `/exception`（例: "annotated by transcript or proteomic data"）が付いた CDS の不一致は、detail に `expected: /exception=…` と記録する |
| 転写産物のモデル | ゲノムから組み立てた配列と、転写産物自身の配列（`--fasta rna.fna`）を比べる | 長さが同じで、塩基置換が5%以下なら ok（座標は保たれる）。転写産物が長く、はみ出した部分が9割以上 A なら、ポリA鎖とみなしてその手前で比べる。それ以外の長さの違いは、挿入・欠失ありとして mismatch |
| alignment | ブロックごとに同一塩基数を数える | `num_mismatch` があれば、それと一致するかを判定する |

CDS の検証では、次の特殊ケースを考慮する（Ensembl GRCh38 release 116 の全 CDS、約37万件がすべて ok になることを確認した）。

- **CDS のタンパク質長**: 公開タンパク質配列が手元にあれば、その長さを使う（推定では、終止コドンのない不完全な CDS の最後の残基が欠けるため）。
- **先頭の `X`**: phase が0以外で、公開配列の先頭が `X` なら、Ensembl の慣習（欠けた先頭コドンを1残基とする）とみなし、`leadingPartialCodon` で写像を作る（spec-core §4.2）。これがないと、ヒトの Ensembl で7,700件以上が1残基ずれていた。
- **翻訳開始コドン**: 5' 側が完全な CDS で、先頭のコドンが翻訳表の開始コドン、または公開配列の1番目が M なら、M として読む（GTG、ACG などの AUG 以外の開始。Ensembl で82件）。
- **読み替えられた終止コドン**: 公開配列が U（セレノシステイン）か O（ピロリシン）で、ゲノムの翻訳が終止コドンになっている位置は、一致とみなして記録する（Ensembl の GFF3 には `transl_except` がない）。
- **翻訳表の推定**: 翻訳表が明示されていない場合（Ensembl の GFF3 のミトコンドリアなど）で不一致になったときは、公開配列と完全に一致する翻訳表を探して採用し、`translTable` 属性と detail に記録する。GBFF では、INSDC の規約どおり `/transl_table` がなければ翻訳表1とし、推定はしない。

- 開始コドン: 5' 側が完全で、NCBI の翻訳表で開始コドンになっているものは M とする（例: ND2 の ATT）。
- `/transl_except`: 位置を**コアの逆写像で残基番号に変換**して、アミノ酸を置き換える（例: GPX1 の Sec が49番残基）。TERM はタンパク質の外側なので無視する。アミノ酸名は大文字と小文字を区別しない（GBFF は `OTHER`、NCBI の GFF3 は `Other`）。
- 翻訳表: NCBI の `gc.prt` から生成した全27表を使う（`ingest/src/genetic-codes.ts`）。
- 配列が手元にない場合は `skipped` とする。
- **検証の根拠（`basis`）**: 配列と全体で照合したものは `full`、終止コドンがないことだけを確かめたものは `partial`。経路探索では、`full` で ok の edge だけを、exception の目印があっても照合済みとして扱う（spec-service §3）。
- `/exception` の付いた feature の不一致は、CDS でも転写産物でも、detail に `expected: /exception=…` と記録する。

---

## 7. テストに使った実データ（`ingest/test/fixtures/`、2026-09-18 に NCBI から取得）

| データ | 確かめたこと |
|---|---|
| NC_045512.2（SARS-CoV-2、GBFF と GFF3） | slippage、CDS 12件 |
| NC_012920.1（ヒトのミトコンドリアゲノム、GBFF と GFF3） | 翻訳表2、不完全な終止コドン、代替開始コドン、原点をまたぐ D-loop、CDS 13件 |
| NC_001405.1（ヒトアデノウイルス C、GBFF と GFF3） | マイナス鎖でスプライスされた CDS、CDS 38件 |
| NC_002127.1（環状プラスミド） | 翻訳表11 |
| NM_000581.4（GPX1） | セレノシステイン |
| NM_014739.3 と NP_055554.1（GenPept） | `/coded_by` が mRNA レコードの CDS と一致すること |
| NC_000003.12 の `cDNA_match` 12行（GRCh38）と関係する配列 | マイナス鎖の I（挿入）とプラス鎖の D（欠失）を含む alignment edge |

**GFF3 と GBFF の一致**: 3種の生物の CDS 63件で、Location ID、ブロック列、codon_start、aaLength がすべて一致した。

**端から端までの検証**: 全フィクスチャのタンパク質の全残基（2万残基以上）について、「コアでゲノムの位置に変換 → 配列を取り出す → 翻訳」した結果が `/translation` と一致した。

---

## 8. 大規模データの取り込み（v0.2）

設計は [scaling.md](scaling.ja.md)。

- **ストリーミング**: `ingestGff3File` と `ingestGenBankFile` は、行単位で読む（`.gz` はそのまま展開しながら読む）。GFF3 は、同じ ID の行を最大1000 feature の範囲でまとめる（NCBI の GFF3 では常に隣り合っている）。確定したあとに同じ ID が再び現れた場合は警告を出す。ストリーミング時は、`##FASTA` 節を読まない（`--fasta` で渡す）。
- **出力先（Sink）**: `MemorySink`（テストと小規模データ用）、`JsonlSink`、`SqliteSink`。アダプタのロジックは出力先によらず同じ。
- **FASTA**: 64MB を超えるファイルや、`.fai` があるファイルは、`.fai`（samtools と同じ形式。なければ自動で作る）を使ってランダムアクセスする。
- **保存先**: SQLite（`node:sqlite`）。スキーマは scaling.md §4（v0.2 で edge に `basis` 列を追加し、スキーマのバージョンを2にした）。R*Tree は、全行を書き込んだあとに位置順に構築する。
- **exon**: 既定では annotation として保存しない（`--all-annotations` で保存する）。

---

## 9. FASTA アダプタ（配列の同一性）

- `ingestFastaFile`: 各配列を、長さ・refget ダイジェスト・MD5 を持つ SequenceRecord として出力する（配列そのものは保存しない）。
- ヘッダの解釈: `sp|P07203|…`・`tr|…` → `uniprot:`（アイソフォームの接尾辞は残す）、`ENSP…` → `ensembl:`、`101m_A mol:protein` → `pdb:101M.A`（wwPDB の `pdb_seqres.txt`。`mol:na` は読み飛ばす）、それ以外は accession として解釈する。
- 同じ配列が先に（ダイジェストなしで）登録されていた場合、SQLite の保存先は、空の項目だけを後のレコードで埋める。
- CLI: `.fa`、`.faa`、`.fasta`（`.gz` も可）を入力に渡すと、このアダプタで取り込む。

## 10. SIFTS アダプタ

- 入力: `uniprot_segments_observed.tsv(.gz)`（EBI）。1行 = UniProt の区間と PDB 鎖の SEQRES の区間。
- edge: `uniprot:<SP_PRIMARY>` → `pdb:<PDB>.<CHAIN>` の alignment。座標は SEQRES の番号（`RES_BEG..RES_END`、label_seq_id に相当）で、著者番号（`PDB_BEG-PDB_END`）は `authorNumbering` 属性に残す。CHAIN は著者の chain ID（auth_asym_id）。
- 同じ（UniProt, 鎖）の連続する行を1つの edge にまとめる。UniProt 側と SEQRES 側で長さが違う行は、1対1にできないので読み飛ばす（ヒトで867行、0.06%）。
- 自己検証: UniProt と `pdb_seqres.txt` の配列で、一致する残基の数を数える。5割未満なら mismatch（人工的な変異は多くても数残基なので、座標のずれだけを検出する）。ヒトでは24万 edge のうち、mismatch は137件（いずれも短い区間）。
- CLI: ファイル名に `sifts` か `uniprot_segments` を含む `.tsv(.gz)` を、このアダプタで取り込む。`--sifts-known-only` を付けると、`--fasta` で与えた UniProt 配列にある accession の行だけを取り込む。

## 11. 配列のタグと MANE（v0.4）

- SequenceRecord に `tags`（文字列の一覧。例: `MANE Select`）と `gene`（遺伝子名）を追加した。保存先のスキーマはバージョン4。同じ保存先の中では、空の項目だけを後のレコードで埋める。保存先をまたぐとタグは和集合にする（spec-service §1）。
- MANE アダプタ（`ingestManeSummary`）: NCBI の `MANE.GRCh38.*.summary.txt(.gz)` を読み、各遺伝子の RefSeq の NM・NP と、Ensembl の ENST・ENSP に、`MANE Select` または `MANE Plus Clinical` のタグと遺伝子名を付けた SequenceRecord を出力する。長さは0（不明）とし、ほかの保存先の値を使う。CLI は、ファイル名が `MANE…summary.txt` のものをこのアダプタで取り込む。
- MANE の RNA 配列（`refseq_rna.fna`、`ensembl_rna.fna`）を FASTA アダプタで取り込むと、ダイジェストの一致によって、MANE の NM と ENST を同一配列として行き来できるようになる（例: NM_000581.4 と ENST00000419783.3 は899塩基が完全に一致する）。FASTA アダプタは、転写産物の accession（NM・NR・XM・XR、ENST）の分子種を RNA とする。

## 12. 保存先のメタデータと内容の集計（v0.4）

- CLI の指定: `--label`（表示名）、`--taxon`・`--organism`・`--assembly`、`--assembly-report FILE`（NCBI の assembly report の見出しから、学名、taxon、アセンブリ名、アセンブリの accession を読む。`--seqid-map` も同じ情報を読む）。値は保存先の `meta` 表に記録する。
- 生物種は、データからも拾う。GBFF の `source` と NCBI の GFF3 の `region`（`Dbxref=taxon:N`）、UniProt の FASTA のヘッダ（`OS=`、`OX=`）から読む。レコードから作ったタンパク質は、その生物種を引き継ぐ。
- 保存先を閉じるときに、内容の集計（`summary`）を `meta` 表に記録する。内容は、分子種ごとの配列の数、種類ごとの edge の数、ブロックの数、アノテーションの数、生物種ごとの配列の数、試せる例（検証済みの CDS のタンパク質の残基と、そのタンパク質全体、アライメントの位置）。集計を持たない古い保存先では、サービスがその場で計算する。

## 13. 生物種を増やしたときの検証（2026-09-18）

| 生物種 | データ | CDS の自己検証 |
|---|---|---|
| マウス | RefSeq GRCm39（GFF3、ゲノム・RNA・タンパク質の配列）、RefSeq RNA、UniProt UP000000589 | 説明のつかない不一致0件（exception 付き37件） |
| シロイヌナズナ | RefSeq TAIR10.1、RefSeq RNA、UniProt UP000006548 | 説明のつかない不一致0件（exception 付き64件。葉緑体・ミトコンドリアの RNA 編集など） |
| ゼニゴケ | **INSDC のみ**（GCA_003032435.1 の GenBank ファイル）、UniProt UP000244005 | 24,674件すべて ok |

- ゼニゴケでは、UniProt → 同一配列の INSDC タンパク質 → scaffold と、完全一致の経路でゲノムに到達できる（例: `uniprot:A0A2R6VYA2:50` → `insdc:KZ773422.1:208..210`）。
- シロイヌナズナの RefSeq RNA（GenBank）には、`/translation` が CDS と合わないレコードが4件あった（例: NM_117687.2。CDS は619塩基だが、翻訳配列は844残基で、`J` を含む）。同じ翻訳配列が4件に入っており、元データの誤りとみられる。アダプタは警告を出して取り込まない。
- 大きな gzip の FASTA（例: マウスの RNA、展開すると約5.4億文字を超える）は、CLI が隣に展開したファイルを作り、`.fai` でランダムアクセスする。

## 14. UCSC chain アダプタと、chunk による保存（v0.5）

- 入力: UCSC の chain 形式（`*.over.chain(.gz)`）。ファイル名が `.chain` または `.chain.gz` で終わるものを、このアダプタで取り込む。chain の配列名（`chr1` など）は UCSC の名前なので、CLI では `--from-report`（変換元のアセンブリ）と `--to-report`（変換先）に NCBI の assembly report を指定して、RefSeq の accession に読み替える（`UCSC-style-name` の列を使う）。対応のない配列（alt、未配置の scaffold など）の chain は、警告を出して飛ばす。
- 1本の chain を、`kind: "liftover"`、`directional: true` の edge 1本にする。属性は chain の id、スコア、変換先の向き（`qStrand`）。変換先が `-` の chain は、座標を順方向に直し（`qSize - q - size`）、ブロックを `rev` にする。
- **自己検証**: `--fasta` で両方のゲノム配列を与えると、chain ごとに等間隔に最大20ブロック（各200塩基まで）を取り出して比べる。一致率が0.5以上なら ok（`basis: "partial"`）。整列した相同配列では0.6〜0.9、位置がずれていれば0.25前後になる。hg38 ↔ mm39 では全体で約70%だった。
- **保存（スキーマ5）**: 向きを持つ edge のブロックは、`block` 表ではなく、256ブロックごとの `chunk` 表に入れる。chunk には、両方の配列、元の側の範囲（`lo`、`hi`）、最初のブロックの位置、ブロックの数と、各ブロックを「前のブロックからの差分（元の側、先の側）と長さ」の zigzag LEB128 で並べたバイト列を持つ。元の側の範囲を R*Tree（`chunk_src`）で引き、該当する chunk だけを展開する。ブロック1個あたり約5バイトで、hg38 → mm39 の3,177万ブロックは158MB になる（`block` 表では約7GB を見込んでいた）。読み手はスキーマ4と5を受け付ける。
- 集計（§12）の例には、chain の範囲から30塩基を選ぶ。

## 15. アセンブリの配列名と、アノテーションのないアセンブリ（v0.5）

- `--assembly-report FILE` を指定すると、保存先の `meta` に、アセンブリの配列名（Sequence-Name、UCSC 名、GenBank の accession）から RefSeq の accession への対応（`aliases`、JSON）と、UCSC のデータベース名（`ucsc`。GRC のアセンブリは assembly report に載っていないので、組み込みの表 `UCSC_DATABASES` から決める。GRCh38 → hg38、GRCh37 → hg19、GRCm39 → mm39、GRCm38 → mm10）も記録する。サービスはこれを使い、`hg19:chr7:140453136` のような入力を読み替える（spec-service §2.2）。
- NCBI の assembly report（`*_assembly_report.txt`）そのものも入力にできる。各配列を、RefSeq の accession、長さ、生物種、DNA（ミトコンドリアと葉緑体は環状）の SequenceRecord にする。アノテーションのないアセンブリ（例: GRCh37）を、ゲノムの FASTA（ダイジェスト用）と一緒に1つの保存先にする。
- chain の保存先は、`--from-report` と `--to-report` の両方のアセンブリの配列の記録（種と長さ）も持つ。
- RefSeq の番号のない（INSDC だけの）アセンブリでは、配列を GenBank の accession（`insdc:AP031342.1`）で表す。配列名の対応の値は、名前空間付きの配列の鍵にした（`refseq:NC_000007.13`、`insdc:AP031342.1`。以前の保存先の名前空間なしの値は、サービスが `refseq:` とみなす）。
- assembly report の公開日（`# Date`）を `meta.released` に記録する。既定のアセンブリの選択に使う（spec-service §2.2）。
- `--species-taxon N`: 保存先の taxon が属する種を明示する（亜種や株の taxon を種にまとめる規則で決まらないとき）。

**ヒト GRCh37 ↔ GRCh38（2026-09-18）**: `grch37.sqlite`（GRCh37.p13 の assembly report と genomic.fna）、`chain_hg19ToHg38.sqlite`、`chain_hg38ToHg19.sqlite`（UCSC の `hg19ToHg38.over.chain.gz`、`hg38ToHg19.over.chain.gz`）。ほかに、以前に作った `human.sqlite` には配列名の記録がないので、GRCh38 の配列名だけの保存先 `grch38_names.sqlite` を assembly report から作った（705配列）。

**マウス GRCm38 ↔ GRCm39（2026-09-18）**: `grcm38.sqlite`（GRCm38.p6 の assembly report と genomic.fna）、`grcm39_names.sqlite`（GRCm39 の配列名）、`chain_mm10ToMm39.sqlite`、`chain_mm39ToMm10.sqlite`（UCSC の `mm10ToMm39.over.chain.gz`、`mm39ToMm10.over.chain.gz`）。

## 16. PAF アダプタ（ゲノム全体のアライメント、T3、v0.5）

- 入力: PAF（`*.paf(.gz)`）。minimap2 などのゲノム間アライメントの出力で、CIGAR（`cg:Z`。minimap2 の `-c`）が必須。query を変換元、target を変換先とする（`minimap2 -c 変換先.fa 変換元.fa`）。配列名の読み替えには chain と同じく `--from-report`（query のアセンブリ）と `--to-report`（target）を使い、両方のアセンブリの配列の記録も持つ。
- CIGAR から、隣り合う一致をまとめたブロック列を作る（`=`、`X`、`M` は両方を、`I` は query を、`D`、`N` は target を進める）。`-` 鎖では、target を順方向に、query を逆向きにたどる。
- **1対1に絞る**: PAF は、反復配列やパラログに対して複数の対応を出す。UCSC の liftOver chain と同じく、変換元（query）の各塩基には1つの対応だけを残す。secondary（`tp:A:S`）と、query 上で 1kb 未満のアライメントを除き、スコア（`AS:i`、なければ一致塩基数）の高い順に採る。既に覆われた query の範囲は、後のアライメントから切り取る。target 側の重なりは許す。1本のアライメントを、向きのある `liftover` edge 1本にする（chain と同じ圧縮した保存）。
- 自己検証は chain と同じ（標本のブロックの一致率が0.5以上なら ok）。
- **取り込まないアライメント**（chain も同じ。2026-09-18）: (1) 変換元の配列が両方のアセンブリに属するもの（同じ accession。TAIR10 と TAIR10.1 の核の染色体、GRCh37 と GRCh38 の chrM や多くの未配置 scaffold）。対応は同一性で、アライメントは反復配列などの別の場所を指すだけ。(2) 核の配列とオルガネラのゲノム（ミトコンドリア、葉緑体）の間のもの。核に入り込んだオルガネラ由来の配列（NUMT、NUPT）との対応で、同じ位置ではない。分子の種類は assembly report の Assigned-Molecule-Location/Type で決める。例: TAIR10 のミトコンドリアを minimap2 で TAIR10.1 に並べると、Col-0 の2番染色体にあるミトコンドリア由来の大きな挿入によく一致し、1対1に絞るとそちらが選ばれていた。ヒト → マウスの chain では392本がこれに当たる。
- **再現のための記録**（すべての保存先）: 入力ファイルの MD5（`meta.inputs_md5`）、`--fasta` の配列ファイルの MD5（`meta.sequences_md5`。アライメントの元になったゲノム）、TogoCoord の commit（`meta.togocoord`。未 commit の変更があれば `+local changes`）、`--method TEXT`（入力の作り方。アライナーの版と引数、絞り方）。Loaded data に表示する。

**ゼニゴケ v3.1 ↔ v7.1（2026-09-18）**

| | v3.1 → v7.1 | v7.1 → v3.1 |
|---|---|---|
| minimap2 2.31-r1302 `-c --eqx -x asm5 -t 8` | 27秒、4.3GB | 56秒、5.3GB |
| PAF の行 | 6,372 | 7,150 |
| 残したアライメント（1対1） | 4,343（15,175ブロック） | 5,858（28,766ブロック） |
| よりよいアライメントに覆われて除いた query の塩基 | 55,948 | 1,781,960 |
| 標本の一致率 | 99.53% | 98.51% |
| 保存先 | `mp_v31_to_v71.sqlite` 2.7MB | `mp_v71_to_v31.sqlite` 3.4MB |

## 17. BED アダプタと、ID で引ける注釈（fanta.bio の CRE、v0.5）

- 入力: BED（`*.bed(.gz)`）。3〜12列の標準の列と、その後の独自の列。ゲノム上の領域を注釈（annotation）として保存する。BED12 のブロックは join に、`-` 鎖は complement にする。染色体名（`chr1`、`1` など）は `--assembly-report` の配列名で読み替える。assembly report に載らないパッチや未配置 scaffold の UCSC 名（`chr11_GL456060_alt`、`chr1_KI270706v1_random`）は、名前の中の GenBank の accession で引く（`lookupSeqid`。chain と PAF の配列名にも使う）。
- `--bed-type TYPE`（注釈の種類。既定 `region`）、`--bed-columns NAME,...`（標準の列の後の列の名前。`attributes` は `key:value|key:value` を属性に分ける）。BED の name 列は属性 `ID` にする。
- `--id-namespace NS`: 注釈を `NS:ID` で入力できるようにする（spec-service §6）。保存先に `annotation_id` 表（ID → 注釈）を作る。`--link URL`: `{id}` を ID に置き換えた URL を、注釈のリンクにする。

**fanta.bio CRE v1.2.1（2026-09-18）**: `https://data.fanta.bio/cre/v1.2.1/` の BED9+2（10列目が CRE 名、11列目が `directionality:…|class:PLA/ELA`）。prefix は bioregistry の `fanta`（パターン `^FC(HS|MM)_\d+$`、`https://fanta.bio/cre/$1`）。

| 保存先 | 入力 | 領域 | 取り込めなかった領域 | 大きさ | 時間 |
|---|---|---|---|---|---|
| `fanta_human_hg38.sqlite` | `human-CREv1.2.1.hg38.cre-peaks.bed.gz`、`--assembly-report` GRCh38.p14 | 513,895 | 0 | 196MB | 4秒 |
| `fanta_mouse_mm10.sqlite` | `mouse-CREv1.2.1.mm10.cre-peaks.bed.gz`、`--assembly-report` GRCm38.p6 | 307,621 | 206（mm10 のパッチのうち、GRCm38.p6 で版の上がったもの。座標が同じとは限らないので取り込まない） | 115MB | 3秒 |

```
togocoord-ingest --db fanta_mouse_mm10.sqlite --assembly-report GCF_000001635.26_GRCm38.p6_assembly_report.txt \
  --bed-type CRE --bed-columns Name,attributes --id-namespace fanta --link 'https://fanta.bio/cre/{id}' \
  mouse-CREv1.2.1.mm10.cre-peaks.bed.gz
```

活性の表（TPM）と JSONL の注釈（関連遺伝子、TF の結合など、446MB）は取り込まない。座標の対応には要らず、fanta.bio へのリンクで参照できる。

## 18. 同一配列のないタンパク質のアライメント（T2、v0.5）

UniProt のエントリは、同一配列（refget ダイジェスト）の RefSeq・Ensembl・INSDC のタンパク質を介してゲノムに届く。1残基でも違えば同一配列がなく、ID の関係（UniProt の相互参照）はあっても、どの残基がどの残基に当たるかという座標の関係がないので、届かない。例: マウス Nras の Swiss-Prot `P08556` は、GRCm39 の翻訳（`NP_035067.2`）と2残基違う（168 L/M、184 S/L）。このアダプタは、ID の関係で候補を選び、配列を並べて座標の関係（edge）を作る。

- 入力: UniProt の `*_idmapping_selected.tab(.gz)`（by_organism）と、`--fasta` の UniProt と候補のタンパク質の配列。
- **対象**: 候補の配列のどれとも同一でない UniProt のエントリ（アイソフォームを含む）。同一のものは同一配列の経路で届くので扱わない。
- **候補**（順に、見つかった段階で止める）: (1) エントリの行の RefSeq（4列目）、EMBL-CDS（INSDC のタンパク質、18列目）、Ensembl_PRO（21列目）。(2) 同じ遺伝子（GeneID の3列目、Ensembl gene の19列目）の他のエントリが挙げるタンパク質。(3) 同じ UniRef90 クラスタ（9列目）の他のエントリが挙げるタンパク質。UniProt は RefSeq を配列の一致するエントリに付けるので、`P08556` の行は RefSeq も遺伝子も空で、`NP_035067.2` は同じ UniRef90_P08556 の TrEMBL `A0A0G2JDN6` の行にある。
- **アライメント**: 両方に1回ずつ現れる6残基を錨にし、両方で順に並ぶ最長の錨の列（LIS）をとる。錨の間を Needleman-Wunsch（一致 2、不一致 -1、ギャップ -2）で埋める。両端は端のギャップを無料にする。埋めた区間は、錨の間なら長さが同じ（置換だけ）か半分以上一致、端なら半分以上一致かつ3残基以上一致のときだけ残す（違う最初・最後のエキソンを偶然に並べない）。ギャップのない列の連続をブロックにし、半分未満しか一致しない連続は捨てる。4百万セルを超える区間は並べない。
- 候補のうち一致した残基の最も多いものを採り、エントリ → 候補の `alignment` の edge を1本作る（provenance の adapter は `protein-alignment`）。属性: 一致率（`identity`）、エントリのうち並んだ割合（`coverage`）、置換（`substitutions`。`エントリの位置/候補の位置:残基>残基` を最大200）、候補の選び方（`candidates`）。一致率0.9以上かつ被覆0.5以上なら ok、そうでなければ mismatch（経路では `approximate`）。

| 種 | UniProt のエントリ | 同一配列あり | 候補なし | 並べた（ok） | うち同じ遺伝子 | うち同じ UniRef90 | 並ばず |
|---|---|---|---|---|---|---|---|
| ヒト（RefSeq と Ensembl） | 169,651 | 162,125 | 947 | 6,556（6,470） | 362 | 117 | 23 |
| マウス（RefSeq） | 63,328 | 29,794 | 1,819 | 31,257（30,009） | 24,983 | 926 | 458 |
| シロイヌナズナ（RefSeq） | 41,596 | 40,369 | 259 | 965（949） | 17 | 113 | 3 |

マウスで同一配列が少ないのは、Ensembl を読み込んでおらず、Ensembl 由来の TrEMBL エントリが RefSeq と一致しないため。それらは同じ遺伝子の RefSeq に並べてつなぐ。ゼニゴケは by_organism に ID mapping がないので対象外（UniProt の99.2%が同一配列でつながる）。

## 19. 違う塩基の記録（スキーマ6、2026-09-19）

保存先は配列そのものを持たない（容量、配布元との同期、再配布の条件、座標のサービスの役割を越えるため）。その代わり、**違いのある位置だけ**を持つ。

- 同じ種のアセンブリ間の chain と PAF（両方の assembly report の taxon が同じか、一方が他方の二名法の名前に種より下の階級を付けた名前）では、取り込み時に全ブロックの塩基を実配列で比べ、違う塩基の位置（変換元の0始まりの位置）、変換元の塩基、変換先の塩基（変換元の向きに直したもの）を `mismatch` 表に記録する（`(seq, pos)` の索引）。edge の属性 `mismatches` に件数。`N` は違いとしない。
- 種の違う chain（ヒト ↔ マウス）では、塩基の多くが違うので記録しない（spec-service §14 の注意で示す）。
- T2 のタンパク質のアライメントは、置換を edge の属性 `substitutions` に持つ（§18）。

| 保存先 | 違う塩基 | 大きさ |
|---|---|---|
| `chain_hg19ToHg38` / `chain_hg38ToHg19` | 66,673 / 185,228 | 3.5MB / 16MB |
| `chain_mm10ToMm39` / `chain_mm39ToMm10` | 26,056 / 39,541 | |
| `chain_tair10.1ToTair10` / `tair10_to_tair10.1` | 233 / 165 | |
| `mp_v31_to_v71` / `mp_v71_to_v31` | 42,516 / 336,134 | / 15MB |

スキーマ6の読み手は4〜6を受け付ける（`mismatch` 表のない保存先は、違いを記録していないだけ）。

## 20. 注釈のあるアセンブリを中心にしたアライメント（スター型、2026-09-20）

同じ種にアセンブリが増えると、総当たりのアライメントは N(N−1) 本になる。代わりに、**その種の注釈の最も多いアセンブリ（既定のアセンブリ）を中心にして、他のアセンブリとの間だけを両方向に並べる**（2(N−1) 本）。中心にない組み合わせは、中心を経由して変換する（spec-service §2.2 の「同じ種ではアセンブリを2回までまたぐ」）。

- 中心の選び方は既定のアセンブリと同じ（注釈があるもののうち、公開の新しいもの、次に注釈の多いもの）。ゼニゴケでは MpTak_v7.1（Tak-1 の常染色体に Tak-2 由来の chrU と Tak-1 の chrV、オルガネラを含むので、雄株・雌株のどちらのアセンブリも受けられる）。
- 必要な組には、後から直接のアライメント（弦）を足せる。中心を増やすより本数も経路も増えない。
- minimap2 の設定は、配列のずれの大きさで選ぶ（minimap2 の目安）: `asm5`（約0.1%。同じ株の版違い）、`asm10`（約1%。同じ種の別の株）、`asm20`（約5%。同じ種の離れた系統や亜種）。5%を大きく超える相手（別の種）には向かず、配布されている chain（UCSC は lastz で作る）を使う。

**ゼニゴケ（2026-09-20）**: 中心 MpTak_v7.1、spoke は v3.1、MpTak2_v7.1（Tak-2）、ASM993635v2（v5.1）、cmMarPoly1.2（注釈なし）、Mp_v4。

| spoke | 設定 | 残したアライメント（→中心 / 中心→） | 標本の一致率 |
|---|---|---|---|
| v3.1（GCA_003032435.1） | asm5 | 4,342 / 5,843 | 99.5% / 98.5% |
| MpTak2_v7.1（GCA_037833965.1） | asm5 | 1,461 / 1,538 | 98.7% / 98.6% |
| ASM993635v2（GCA_009936355.2、v5.1） | asm5 | 1,344 / 1,283 | 99.6% / 99.1% |
| cmMarPoly1.2（GCA_965642975.2） | asm20 | 13,501 / 11,314 | 93.4% / 93.7% |
| Mp_v4（GCA_001641455.1） | asm5 | 7,346 / 11,640 | 98.7% / 98.2% |

**中心に移る区間の割合**（`verify-genome-pair.ts`、spoke のゲノムから無作為に1,500区間）: v5.1 97.9%、MpTak2_v7.1 91.6%、Mp_v4 82.9%、cmMarPoly1.2 63.9%。cmMarPoly1.2 は別系統で、塩基が2〜3%ずれ、ゲノムも 265Mb と中心（248Mb）より大きい。設定を上げると移る割合も上がる（asm5 39.6% → asm10 54.7% → asm20 63.9%）。アライメント自体が覆うのは asm10 で 60.3%、asm20 で 71.5% なので、1対1に絞る処理ではなくアライメントの段階での限界。

**Mp_v4 を入れた理由（2026-09-20）**: 2016年の scaffold レベルのアセンブリ（Oxford、Tak-1 と Tak-2 のプール、205.7Mb）で、ゲノムとしては世代交代している。しかし UniProt の UP000077202（17,951件）がこのアセンブリに基づいており、到達できないゼニゴケのタンパク質として最大の塊だった。収載済みのタンパク質と同一配列なのは 20.1% だけで、残り8割は孤立した節点になる。Mp_v4 を入れると、17,950/17,951 がその CDS のタンパク質と同一になり、96% が中心のゲノム座標まで到達する（`uniprot:A0A176VNS3:10` → `insdc:OAE22043.1:10` → Mp_v4 の scaffold → `insdc:AP031346.1:complement(22320360..22320362)`）。遺伝子モデルは古いので、Mp_v4 のタンパク質のうち中心のものと同一なのは 18.8%（3,373 / 17,956）にとどまり、残りはゲノムのアライメント経由で中心に移る。`asm10` にしてもアライメントは 0.6% しか増えない（問い合わせ側の塩基で 179.2Mb 対 178.2Mb）ので `asm5` のままとした。

**違いの記録の上限**: 違う塩基が、揃った塩基の1%（`MAX_MISMATCH_RATE`）を超えるアライメントでは、位置を記録せず、割合（属性 `mismatchRate`）だけを残す。ほぼ同じ配列で「どこが違うか」を示すための記録であり、離れた相手では例外ではなくなるため。cmMarPoly1.2 では、これで保存先が 291MB から 11.6MB になった（記録が残るのは、違いの少ない一部の領域だけ）。

**ゼニゴケで入れていないアセンブリ**（2026-09-20）

| アセンブリ | 理由 |
|---|---|
| MpTak1_v7.1（GCA_037833805.1） | 配列が標準ゲノム（中心）と同じ accession（AP031342〜AP031350）で、CDS のタンパク質 20,184件（BFI03036.1 など）も中心の 20,412件にすべて含まれる。配列もIDも変換も増えない。プロテオーム UP001452901 も UniProtKB に配列を持たない |
| Col-CEN v1.2（アラビドプシス） | INSDC の accession がなく、配列の鍵を決められない |
| cmMarPoly1.2 以外の注釈のないアセンブリ（Marpolrud_CA_v1、亜種2つ、ASM1997375v1） | 注釈がなく、中心との対応も自前の計算が必要。cmMarPoly1.2 を代表の実例とする |

MpTak2_v7.1（Tak-2）は、標準ゲノムと同一配列のタンパク質が 81.8%（16,641 / 20,354）で、残りは Tak-2 にしかない。chrU は標準ゲノムが Tak-2 のものを取り込んでいるので共通。Tak-2 の座標で書かれたデータを受けるために入れている。

**ヒト T2T-CHM13v2.0（2026-09-20）**: `chm13.sqlite`（`GCF_009914755.1` の RefSeq 注釈 2025-08。GFF3、ゲノム、RNA、タンパク質。CDS の自己検証は説明のつかない不一致0件、`/exception` 付き36,846件）、`chain_hs1ToHg38.sqlite`（13,185 chain、856,771ブロック、違う塩基3,076,803件、119MB）、`chain_hg38ToHs1.sqlite`（7,679 chain、823,316ブロック、3,365,500件、128MB）。UCSC 名は `hs1`（`UCSC_DATABASES` に追加）。ヒトの中心は GRCh38 のままで、CHM13 と GRCh37 は spoke（`hs1:chr7:142067515` → GRCh38 の `NC_000007.14:140753336` → GRCh37 の `NC_000007.13:140453136`）。CHM13 は自身の注釈を持つので、タンパク質や UniProt には中心を経由せず1段で届く。
