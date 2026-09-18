# TogoCoord - convert mRNA to protein

## Description

* req.
  * [./togocoord_convert_ensembl_transcript_ensembl_protein](./togocoord_convert_ensembl_transcript_ensembl_protein)
  * [./togocoord_convert_refseq_rna_refseq_protein](./togocoord_convert_refseq_rna_refseq_protein)
  * [./togocoord_protein_id_resolver](./togocoord_protein_id_resolver)

## Parameters

* `source`
  * default: NM_014739
  * example: ENST00000531224
* `location`
  * default: 1..900
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example:

## `uniprot`
```javascript
async ({source, location, target}) => {
  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  options.body = "sources=" + encodeURIComponent('["' + source + '"]') + "&mode=mrna&location=" + encodeURIComponent(location);
  if (target) options.body += "&target=" + target;
  let json = [];
  if (source.match(/^(?:(?:(?<id1>ENS[A-Z]*T\d{11})(?:\.\d+)?)|(?<id2>FBtr\d{7})|(?:(?<id3>MGP_[A-Za-z0-9]+_T\d{7})(?:\.\d+)?)|(?<id4>LRG_\d+t\d+(?:-\d)?)|(?<id5>(?:[A-Za-z0-9][A-Za-z0-9_]+)\.t?\d+[a-z]?(?:\.\d+)?)|(?<id6>Y[A-Z]{2}\d{3}[A-Z](?:-[A-Z])?(?:_[a-z]{1,3}RNA)?)|(?<id7>Q\d{4}_[a-z]{1,3}RNA)|(?<id8>[A-Z]{3}\d{1,3}(?:-\d)?_[a-z]{1,3}RNA)|(?<id9>snR\d+-?[A-Za-z]?_sno?RNA)|(?<id10>t[A-Z]\([ACGUX]{3}\)[A-Z]\d?_tRNA))$/)) {
    json = await fetch("./togocoord_convert_ensembl_transcript_ensembl_protein", options).then(r => r.json());
  }else {
    json = await fetch("./togocoord_convert_refseq_rna_refseq_protein", options).then(r => r.json());
  }
  //console.log(json);
  let protein_id = json[0].jsonld.id.match(/\/([^\/:]+):/)[1];
  let proteins = await fetch("./togocoord_protein_id_resolver?protein_accession=" + protein_id).then(r => r.json());
  let res = [];
  Object.keys(proteins).forEach(key => {
    for (let id of proteins[key]) {
      if (target && target != id) continue;
      let prefix = key;
      if (key == "isoform") prefix = "uniprot";
      if (key == "insdc_contig") prefix = "insdc";
      let data = JSON.parse(JSON.stringify(json[0]));
      data.jsonld.id = data.jsonld.id.replace(protein_id, id);
      if (data.jsonld.location.begin) data.jsonld.location.begin.reference = prefix + ":" + id;
      if (data.jsonld.location.end) data.jsonld.location.end.reference = prefix + ":" + id;
      if (data.jsonld.location.reference) data.jsonld.location.reference = prefix + ":" + id;
      data.jsonld.location.member?.forEach((m, i) => {
        if (m.begin) data.jsonld.location.member[i].begin.reference = prefix + ":" + id;
        if (m.end) data.jsonld.location.member[i].end.reference = prefix + ":" + id;
        if (m.reference) data.jsonld.location.member[i].reference = prefix + ":" + id;
      });
      data.correspondence.reference = data.correspondence.reference.replace(protein_id, id);
      res.push(data);
    }
  });
  return res;
} 
```