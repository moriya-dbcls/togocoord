# TogoCoord - convert protein to PDB

## Description

* req.
  * [./togocoord_convert_uniprot_pdb](./togocoord_convert_uniprot_pdb)
  * [./togocoord_protein_id_resolver](./togocoord_protein_id_resolver)

## Parameters

* `source`
  * default: Q9NYF8
  * example: ENSP00000435210, NP_055554, AAI32781
* `location`
  * default: join(30..336,400..960)
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example: 7RJN

## `uniprot`
```javascript
async ({source, location, target}) => {
  let accession = source.replace(/\s/g, "").replace(/\.\d+$/, "");
  const json = await fetch("./togocoord_protein_id_resolver?protein_accession=" + accession).then(r => r.json());
  let uniprot = json.uniprot[0] || false;
  
  if (!uniprot) return [];

  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  options.body = "source=" + uniprot + "&location=" + encodeURIComponent(location);
  if (target) options.body += "&target=" + target;
  return await fetch("./togocoord_convert_uniprot_pdb", options).then(r => r.json());
}
```