# TogoCoord - convert protein to reference genome

## Description

* req.
  * [./togocoord_convert_refseq_protein_refseq_rna](./togocoord_convert_refseq_protein_refseq_rna)
  * [./togocoord_protein_id_resolver](./togocoord_protein_id_resolver)

## Parameters

* `source`
  * default: NP_001123985
  * example: Q9NYF8, Q9NYF8-2, ENSP00000435210, AT2G33550
* `location`
  * default: 1..800
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)

## `target_align_ranges`
* 
```javascript
async ({source, location, target}) => {
  let ensembl_ids = [];
  if (source.match(/^(?:(?:(?<id1>ENS[A-Z]*P\d{11})(?:\.\d+)?)|(?<id2>FBpp\d{7})|(?:(?<id3>MGP_[A-Za-z0-9]+_P\d{7})(?:\.\d+)?)|(?<id4>LRG_\d+p\d+(?:-\d)?)|(?<id5>(?:[A-Za-z0-9][A-Za-z0-9_]+)\.t?\d+[a-z]?(?:\.\d+)?)|(?<id6>Y[A-Z]{2}\d{3}[A-Z](?:-[A-Z])?)|(?<id7>Q\d{4}))$/)) {
    ensembl_ids.push(source.replace(/\.\d+$/,""));
  } else {
    let json = await fetch("./togocoord_protein_id_resolver?protein_accession=" + source).then(r => r.json());
    //console.log(json);
    Object.keys(json).forEach(key => {
      if (key == "ensembl") {
        ensembl_ids.push(...json[key]);
      }
    });
  }

  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  let body = "location=" + encodeURIComponent(location);
  //if (target) body += "&target=" + target;

  console.log(JSON.stringify(ensembl_ids));
  options.body = body + "&mode=genome&sources=" + encodeURIComponent(JSON.stringify(ensembl_ids));
  console.log(options.body);
  return await fetch("./togocoord_convert_ensembl_protein_ensembl_genome", options).then(r => r.json());
}
```
