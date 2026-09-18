# TogoCoord - convert protein to genome contig

## Description

* req.
  * [./togocoord_convert_refseq_protein_refseq_rna](./togocoord_convert_refseq_protein_refseq_rna)
  * [./togocoord_protein_id_resolver](./togocoord_protein_id_resolver)

## Parameters

* `source`
  * default: NP_055554
  * example: Q9NYF8, Q9NYF8-2, ENSP00000435210
* `location`
  * default: 1..920
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example: CH471051

## `target_align_ranges`
* 
```javascript
async ({source, location, target}) => {
  let proteins = [];
  let json = await fetch("./togocoord_protein_id_resolver?protein_accession=" + source).then(r => r.json());
  Object.keys(json).forEach(key => {
    if (key == "insdc_contig") {
      proteins.push(...json[key]);
    }
  });
   console.log(proteins);

  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  options.body = "sources=" + encodeURIComponent(JSON.stringify(proteins)) + "&location=" + encodeURIComponent(location);
  if (target) options.body += "&target=" + target;
  //console.log(options.body);
  return await fetch("./togocoord_convert_refseq_protein_refseq_rna", options).then(r => r.json());
}
```
