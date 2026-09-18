# TogoCoord - convert protein to mRNA

## Description

* req.
  * [./togocoord_convert_refseq_protein_refseq_rna](./togocoord_convert_refseq_protein_refseq_rna)
  * [./togocoord_convert_ensembl_protein_ensembl_genome](./togocoord_convert_ensembl_protein_ensembl_genome)
  * [./togocoord_protein_id_resolver](./togocoord_protein_id_resolver)

## Parameters

* `source`
  * default: ENSP00000435210
  * example: Q9NYF8, Q9NYF8-2, NP_055554, AT2G33550
* `location`
  * default: 1..100
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example: NM_014739

## `target_align_ranges`
* 
```javascript
async ({source, location, target}) => {
  let refseq_ids = [];
  let insdc_ids = [];
  let ensembl_ids = [];
 // if (source.match(/^[ANXYWZ]P_\d+(?:\.\d+)?$/)) {
 //   refseq_ids.push(source.replace(/\.\d+$/,""));
 // } else {
    let json = await fetch("./togocoord_protein_id_resolver?protein_accession=" + source).then(r => r.json());
    console.log(json);
    Object.keys(json).forEach(key => {
      if (key == "refseq" || key == "insdc") {
        refseq_ids.push(...json[key]);
      } else if (key == "insdc") {
        insdc_ids.push(...json[key]);
      } else if (key == "ensembl") {
        ensembl_ids.push(...json[key]);
      }
    });
 // }

  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  let body = "location=" + encodeURIComponent(location);
  let target_ensembl_flag = false;
  if (target) {
    if (target.match(/^(?:(?:(?<id1>ENS[A-Z]*T\d{11})(?:\.\d+)?)|(?<id2>FBtr\d{7})|(?:(?<id3>MGP_[A-Za-z0-9]+_T\d{7})(?:\.\d+)?)|(?<id4>LRG_\d+t\d+(?:-\d)?)|(?<id5>(?:[A-Za-z0-9][A-Za-z0-9_]+)\.t?\d+[a-z]?(?:\.\d+)?)|(?<id6>Y[A-Z]{2}\d{3}[A-Z](?:-[A-Z])?(?:_[a-z]{1,3}RNA)?)|(?<id7>Q\d{4}_[a-z]{1,3}RNA)|(?<id8>[A-Z]{3}\d{1,3}(?:-\d)?_[a-z]{1,3}RNA)|(?<id9>snR\d+-?[A-Za-z]?_sno?RNA)|(?<id10>t[A-Z]\([ACGUX]{3}\)[A-Z]\d?_tRNA))$/)) {
      target_ensembl_flag = true;
    }
    body += "&target=" + target;
  }
  
  let res = [];
  if ((refseq_ids[0] || insdc_ids[0]) && !target_ensembl_flag) {
    options.body = body + "&mode=mrna&sources=" + encodeURIComponent(JSON.stringify([...refseq_ids, ...insdc_ids]));
    const json = await fetch("./togocoord_convert_refseq_protein_refseq_rna", options).then(r => r.json());
    res.push(...json);
  }
  if (ensembl_ids[0]) {
    options.body = body + "&mode=mrna&sources=" + encodeURIComponent(JSON.stringify(ensembl_ids));
    const json = await fetch("./togocoord_convert_ensembl_protein_ensembl_genome", options).then(r => r.json());
    res.push(...json);
  }
  return res;
}
```
