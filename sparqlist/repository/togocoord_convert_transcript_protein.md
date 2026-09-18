# TogoCoord - convert transcript to protein

## Description

* req.
  * [./togocoord_convert_ensembl_transcript_ensembl_protein](./togocoord_convert_ensembl_transcript_ensembl_protein)
  * [./togocoord_protein_id_resolver](./togocoord_protein_id_resolver)

## Parameters

* `source`
  * default: ENST00000531224
* `location`
  * default: join(9981..10084,11071..11982,13339..14004,14146..14315,16660..16765,17768..17852,20235..20410,21508..21685,22672..22818,28370..28582,28732..28734)
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
  options.body = "sources=" + encodeURIComponent('["' + source + '"]') + "&mode=transcript&location=" + encodeURIComponent(location);
  if (target) options.body += "&target=" + target;
  let json = await fetch("./togocoord_convert_ensembl_transcript_ensembl_protein", options).then(r => r.json());
  if (!json[0]) return [];
  let ensembl_protein = json[0].jsonld.id.match(/\/([^\/:]+):/)[1];
  let proteins = await fetch("./togocoord_protein_id_resolver?protein_accession=" + ensembl_protein).then(r => r.json());
  let res = [];
  Object.keys(proteins).forEach(key => {
    for (let id of proteins[key]) {
      if (target && target != id) continue;
      let prefix = key;
      if (key == "isoform") prefix = "uniprot";
      if (key == "insdc_contig") prefix = "insdc";
      let data = JSON.parse(JSON.stringify(json[0]));
      data.jsonld.id = data.jsonld.id.replace(ensembl_protein, id);
      if (data.jsonld.location.begin) data.jsonld.location.begin.reference = prefix + ":" + id;
      if (data.jsonld.location.end) data.jsonld.location.end.reference = prefix + ":" + id;
      if (data.jsonld.location.reference) data.jsonld.location.reference = prefix + ":" + id;
      data.jsonld.location.member?.forEach((m, i) => {
        if (m.begin) data.jsonld.location.member[i].begin.reference = prefix + ":" + id;
        if (m.end) data.jsonld.location.member[i].end.reference = prefix + ":" + id;
        if (m.reference) data.jsonld.location.member[i].reference = prefix + ":" + id;
      });
      data.correspondence.reference = data.correspondence.reference.replace(ensembl_protein, id);
      res.push(data);
    }
  });
  return res;
} 
```