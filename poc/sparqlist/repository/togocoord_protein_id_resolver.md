# TogoCoord - resolve protein IDs via UniParc xref

## Description

* [UniParc API](https://www.uniprot.org/api-documentation/uniparc)

## Parameters

* `protein_accession`
  * default: EAW47950
  * example: NP_055554, Q9NYF8-3, AT2G33550

## `uniparc_api`
* 
```javascript
async ({protein_accession}) => {
  let accession = protein_accession.replace(/\s/g, "").replace(/\.\d+$/, "");
  const api = "https://rest.uniprot.org/uniparc/";
  let json = await fetch(api + "search?size=10&query=" + accession).then(r => r.json());
  let upis = [];
  if (json.results.length == 1) {
    upis = [json.results[0].uniParcId] || [false];
  } else { // 複数の UniParc が Hit した場合
    for (let d of json.results) {
      // UniProt ID が完全一致を含む
      if (d.uniProtKBAccessions?.includes(accession)) {
        upis = [d.uniParcId];
      }
      break;
    }
    if (!upis[0]) {
      let count = 0;
      for (let d of json.results) {
        // UniProt entry が多い順
        if (d.uniProtKBAccessions?.length > count) {
          upis.unshift(d.uniParcId);
          count = d.uniProtKBAccessions.length;
        } else {
          upis.push(d.uniParcId);
        }
      }      
    }
  }
  if (!upis[0]) return {};

  let tax = false;
  for (let upi of upis) {
    json = await fetch(api + upi + "/?fields=upi%2Caccession%2Corganism_id").then(r => r.json());
    for (let d of json.uniParcCrossReferences) {
      if (d.id.replace(/\.\d+$/, "") == accession) {
        tax = d.organism?.taxonId || false;
        break;
      }
    }
    if (tax) break;
  }
  console.log(json);
  
  let res = {};
  for (let d of json.uniParcCrossReferences) {
    if (d.organism?.taxonId != tax || !d.active) continue;
    if (d.database == "UniProtKB/Swiss-Prot") (res.uniprot ??= []).push(d.id);
    else if (d.database == "UniProtKB/Swiss-Prot protein isoforms") (res.isoform ??= []).push(d.id);
    else if (d.database == "Ensembl") (res.ensembl ??= []).push(d.id);
    else if (d.database == "RefSeq") (res.refseq ??= []).push(d.id);
    else if (d.database == "EMBL") (res.insdc ??= []).push(d.id);
    else if (d.database == "EMBL_CON") (res.insdc_contig ??= []).push(d.id);
  }
  
  return res;
}
```
