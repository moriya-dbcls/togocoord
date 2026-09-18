# TogoCoord 座標・対応表仕様

基本的にはデータ取ってきて四則演算するだけ  
前提

* source: 変換元  
* target: 変換先

座標仕様

* source, target ともにそれぞれ 1-base begin の連続配列

対応表仕様

* 対応表は target の begin, end の range と、差分 delta で記述  
  * 複数 range の join、相補鎖の complement にも対応  
* delta  は range 毎の source 座標と taeget 座標との差分  
  * ただし、アミノ酸 \<-\> 塩基変換の場合は、塩基座標換算  
    * coord 変換内でもと座標を計算  
    * 問題点：塩基 \-\> アミノ酸の対応表作成計算がややこしくなりがち

delta 例

* 通常 (not comlement)  
  * source 全長が 1..100 で、target の 201..300 に対応している場合  
    * range:{begin: 201, end: 300, delta: \-200}  
      * begin \+ delta \= 1; (source begin)  
      * end \+ delta \= 100; (source end)  
  * source 全長が 1..100 で、target の join(201..300, 401..500) に対応している場合  
    * range:{begin: 201, end: 300, delta: \-200}  (同上)  
      * begin \+ delta \= 1; (range 1: source begin)  
      * end \+ delta \= 100; (range 1: source end)  
    * range:{begin: 401, end: 500, delta: \-300}  
      * delta \= 前のdelta \- gap(300,401) で計算  
      * begin \+ delta \= 101; (range 2: source begin)  
      * end \+ delta \= 200; (range 2: source end)  
* complement  
  * delta は source が \-end .. \-1 となるように計算して格納  
  * 計算の簡易化のため  
    * target 上での数値の大小関係と、source 側の begin, end を理解しやすくするためのの苦肉の策  
  * source 全長が 1..100 で、target の complement(201..300) に対応している場合  
    * range: {begin: 201, end: 300, delta: \-301}  
      * end \+ delta \= \-1; (source begin)  
      * begin \+ delta \= \-100;  (source end)  
  * source 全長が 1..100 で、target の complement(join(201..300,401..500)) に対応している場合  
    * 後ろの range から順に計算  
    * range: {begin: 201, end: 300, delta: \-401}  
      * delta \= 次のdelta \+ gap(301,400)  
      * begin \+ delta \= \-200;  
      * end \+ delta \= \-101;  
    * range:{begin: 401, end: 500, delta: \-501}  
      * begin \+ delta \= \-100;  
      * end \+ delta \= \-1;  \# ここが基準  
* あとは[座標変換 SPARQList](https://sparql-support.dbcls.jp/sparqlist/togocoord_coord_converter) でうまくやる

最近、様々なDICPのデータベースを組み合わせて、データやアノテーションを比較したり俯瞰したりするアプリケーションを作成している  
しかし、生命科学データはゲノム配列やアミノ酸配列など、様々なレイヤーにおいて、それぞれの座標で管理されていることが多い  
そのため、これらのデータ間の対応づけのためには、座標のすり合わせが必要になってくる  
ポスターでは座標変換についての取り組みや、複数のデータベースを組み合わせた使用例を紹介したい