# TogoCoord - convert ensembl protein to mRNA, transcript, reference genome

## Description

## Parameters

* `sources`
  * default: ["ENSP00000252519"]
  * example: ["ENSMUSP00000073626"] ["ENSP00000252519","ENSP00000389326","ENSP00000504103","ENSP00000505992"]
* `mode`
  * default: premrna
  * example: exon, mrna (cdna), premrna, genome
* `target` (Opt.)
  * example: NM_014739

## `ensembl_proteins`
```javascript
({sources}) => {
  return JSON.parse(sources);
}
```

## Endpoint

https://rdfportal.org/ebi/sparql

## `ensts`
```sparql
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX so: <http://purl.obolibrary.org/obo/so#>

SELECT DISTINCT ?protein ?transcript
FROM <http://rdfportal.org/dataset/ensembl>
WHERE {
  VALUES ?protein { {{#each ensembl_proteins}}"{{this}}" {{/each}} }
  [] so:translates_to / dct:identifier ?protein ;
     dct:identifier ?transcript .
}
```

## `mappings`
```javascript
async ({target, mode, ensts}) => {
  const transcripts = ensts.results.bindings.map(d =>  d.transcript.value);
  //console.log(transcripts);
  
  let mappings = [];
  for (let accession of transcripts) {
    const api = "https://rest.ensembl.org/lookup/id/" + accession + "?expand=1&content-type=application/json";
    //console.log(api);
    const tr = await fetch(api).then(r => r.json());

  // --- 基本情報 ---
  const proteinId = tr.Translation.id;
  const txGenStart  = tr.Translation.start;
  const txGenEnd    = tr.Translation.end;
  const chr         = tr.seq_region_name;
  const strand      = tr.strand; // +1 または -1
  let exons         = tr.Exon.slice();
    // strand == -1 のとき、C末側から処理
  if (strand == -1) exons.reverse();

  // --- cDNA ⇔ genome セグメント対応 ---
  let cdnaCursor = 1;
  const segments = exons.map(exon => {
    const len = exon.end - exon.start + 1;
    const seg = {
      cdnaStart:    cdnaCursor,
      cdnaEnd:      cdnaCursor + len - 1,
      genomicStart: exon.start,
      genomicEnd:   exon.end
    };
    cdnaCursor += len;
    return seg;
  });
    
  // genomic Translation.start/end を cDNA 上に射影
  function genomicToCdna(pos) {
    const seg = segments.find(s => pos >= s.genomicStart && pos <= s.genomicEnd);
    if (!seg) throw new Error(`Translation座標 ${pos} がいずれの exon にも含まれません`);
    const offset = pos - seg.genomicStart;
    return seg.cdnaStart + offset;
  }
  const cdsStart = genomicToCdna(txGenStart);
  const cdsEnd   = genomicToCdna(txGenEnd);
    
  // --- ranges を構築 ---
  let ranges = segments.map(seg => {
    // CDS 外の exon はスキップ
    if (seg.cdnaEnd < cdsStart || seg.cdnaStart > cdsEnd) {
      return null;
    }

    // exon 内で CDS 部分のみを抽出
    const cdnaBegin = Math.max(seg.cdnaStart, Math.min(cdsStart, cdsEnd));
    const cdnaEnd   = Math.min(seg.cdnaEnd,   Math.max(cdsStart, cdsEnd));
    
    // source (cds) の begin/end
    let aaCdsBegin = cdnaBegin - cdsStart;
    let aaCdsEnd   = cdnaEnd   - cdsStart;

    // source (aa) の begin/end
    const aaBegin = Math.floor((aaCdsBegin) / 3) + 1;
    const aaEnd   = Math.floor((aaCdsEnd) / 3) + 1;

    // codon index (1～3) の begin/end
    const idxBeg = ((aaCdsBegin) % 3 + 3) % 3 + 1;
    const idxEnd = ((aaCdsEnd) % 3 + 3) % 3 + 1;

    // genome (nt) の begin/end
    const offsetB = cdnaBegin - seg.cdnaStart;
    const offsetE = cdnaEnd   - seg.cdnaStart;  
    const ntBeg = seg.genomicStart + offsetB;
    const ntEnd = seg.genomicStart + offsetE;

    return {
      source: {
        begin: aaBegin,
        end:   aaEnd
      },
      sourceCodonIndex: {
        begin: idxBeg,
        end:   idxEnd
      },
      target: {
        begin: ntBeg,
        end:   ntEnd
      }
    };
  }).filter(x => x !== null);

  // strand -1 のときに、ranges を reverseし、アミノ酸ポジションをN末C末を入れ替える
    if (strand == -1) {
      const aaLen = ranges[ranges.length - 1].source.end;
      ranges = ranges.slice().reverse().map(d => {
        return {
          ...d,
          source: {
            begin: aaLen - d.source.end + 1,
            end: aaLen - d.source.begin + 1
          },
          sourceCodonIndex: {
            begin: 4 - d.sourceCodonIndex.end,
            end: 4 - d.sourceCodonIndex.begin
          }
        };
      });
    }
    // stop codon 分削る
   /*  if (strand == 1) {
      ranges[ranges.length - 1].source.end -= 1;
      ranges[ranges.length - 1].target.end -= 3;
    } else {
      ranges[ranges.length - 1].source.begin -= 1;
      ranges[ranges.length - 1].target.begin += 3;      
    } */

  	let mapping = {
      source: {
        reference: proteinId,
        type: "aa"
      },
      target: {
        reference: "GRCh38-" + chr,
        type: "nt",
        strand: strand === 1 ? "+" : "-"
      }
    };
    
    // genome -> pre-mRNA, cdna (mrna), exon
    if (mode != "genome") {
      mapping.target.reference = tr.id + "-pre_mRNA";
      mapping.target.strand = "+";
   	  if (strand == 1) {
        ranges = ranges.map(d => {
          return {
            ...d,
            target: {
              begin: d.target.begin - tr.start + 1,
              end: d.target.end - tr.start + 1
            }
          }
        });
      } else {
        ranges = ranges.map(d => {
          return {
            ...d,
            target: {
              begin: tr.end - d.target.end + 1,
              end: tr.end - d.target.begin + 1
            }
          }
        });        
      }
    }
    
    if (mode == "cdna" || mode == "mrna") {
      mapping.target.reference = tr.id + "-cDNA";
      ranges = [
        {
          source: {
            begin: 1,
            end: ranges[ranges.length - 1].source.end
          },
          sourceCodonIndex: {
            begin: 1,
            end: 3
          },
          target: {
            begin: ranges[0].target.begin,
            end: ranges[0].target.begin + ranges[ranges.length - 1].source.end * 3 - 1
    	  }
        }
      ];   
    }
    
    if (mode == "exon") {
      for (let [i, d] of ranges.entries()) {
        if (i != 0) {
          const pos = d.target.begin;
          d.target.begin = 1;
          d.target.end = d.target.end - pos + 1;
        }
        mapping.target.reference = tr.Exon[i].id;
        mappings.push({
          ...mapping,
          ranges: [{...d}]
        });
      } 
    } else {
      mappings.push({
        ...mapping,
        ranges
      });
    }
  }
  return mappings;
}
```

