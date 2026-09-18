# TogoCoord アダプタ仕様（v0.3：GBFF / GFF3 / FASTA / SIFTS）

2026-09-18。[design.md](design.md) の §6 のうち、フェーズ2で実装した部分の規則。実装は `ingest/`。

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

設計は [scaling.md](scaling.md)。

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

## 16. PAF アダプタ（ゲノム全体のアライメント、T3、v0.5）

- 入力: PAF（`*.paf(.gz)`）。minimap2 などのゲノム間アライメントの出力で、CIGAR（`cg:Z`。minimap2 の `-c`）が必須。query を変換元、target を変換先とする（`minimap2 -c 変換先.fa 変換元.fa`）。配列名の読み替えには chain と同じく `--from-report`（query のアセンブリ）と `--to-report`（target）を使い、両方のアセンブリの配列の記録も持つ。
- CIGAR から、隣り合う一致をまとめたブロック列を作る（`=`、`X`、`M` は両方を、`I` は query を、`D`、`N` は target を進める）。`-` 鎖では、target を順方向に、query を逆向きにたどる。
- **1対1に絞る**: PAF は、反復配列やパラログに対して複数の対応を出す。UCSC の liftOver chain と同じく、変換元（query）の各塩基には1つの対応だけを残す。secondary（`tp:A:S`）と、query 上で 1kb 未満のアライメントを除き、スコア（`AS:i`、なければ一致塩基数）の高い順に採る。既に覆われた query の範囲は、後のアライメントから切り取る。target 側の重なりは許す。1本のアライメントを、向きのある `liftover` edge 1本にする（chain と同じ圧縮した保存）。
- 自己検証は chain と同じ（標本のブロックの一致率が0.5以上なら ok）。
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
