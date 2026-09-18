# TogoCoord - convert RefSeq (insdc) mRNA (isndc) to protein 

## Description

* req.
  * [./togocoord_location_decoder](./togocoord_location_decoder)
  * [./togocoord_location_encoder](./togocoord_location_encoder)
  * [./togocoord_coord_converter](./togocoord_coord_converter)

## Parameters

* `sources`
  * default: ["NM_014739"]
  * example: ["AF249273"]
* `location`
  * default: 249..3008
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example: NP_055554

## `target_align_ranges`
```javascript
async ({sources, target}) => {
  let proteins = JSON.parse(sources);
  
  let ranges = [];                            
  for (let accession of proteins) {
    const api = "https://togows.org/entry/ncbi-nucleotide/" + accession + ".json";
    const json = await fetch(api).then(r => r.json());
    
    json.forEach(data => {
      data.features.forEach(feature => {
        if (feature.feature != "CDS" || !feature.protein_id || !feature.position) return 0;
        feature.protein_id.forEach(id => {
          const reference = id.replace(/\.\d+$/, "");
          if (target && !target.match(reference)) return 0;
          if (feature.position.match(/^\d+\.\.\d+$/)) {
            const s = feature.position.match(/^(\d+)\.\.(\d+)$/);
            ranges.push({
              reference: reference,
              mode: "nt2aa",
              begin: parseInt(s[1]) * (-1),
              end: parseInt(s[2]) - parseInt(s[1]) + 1 - 3, // - stop codon
              delta: parseInt(s[1]) - 1
            });
          }
        });
      });
    });
  }

  return ranges;
}
```

## `return`
```javascript
async ({location, target_align_ranges}) => {
  if (!target_align_ranges[0]) return [];
  
  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  
  // decode location to coords
  options.body = "location_id=" + encodeURIComponent(location);
  const source_coords = await fetch("./togocoord_location_decoder", options).then(r => r.json());

  // convert coords
  options.body = "source_coords=" + encodeURIComponent(JSON.stringify(source_coords))
        + "&target_align_ranges=" + encodeURIComponent(JSON.stringify(target_align_ranges));
  console.log(JSON.stringify(source_coords));
  console.log(JSON.stringify(target_align_ranges));
  const target_coords = await fetch("./togocoord_coord_converter", options).then(r => r.json());

  // encode json-ld from coords
  const promises = Object.keys(target_coords).map(async (reference) => {
    const new_json = {
      ...(source_coords.assembl && {assembl: source_coords.assembl}),
      reference: reference,
      ...(source_coords.join && {join: source_coords.join}),
      ...((source_coords.complement || target_coords[reference][0]?.joined_complement) && {complement: true}),
      ranges: target_coords[reference]
    };
    //console.log(new_json);
    let target_prefix = "insdc";
    if (reference.match(/^[ANXYWZ]P_\d+/)) target_prefix = "refseq";
    options.body = "data=" + encodeURIComponent(JSON.stringify(new_json)) + "&reference_prefix=" + target_prefix;
    //console.log(options.body);
    return {jsonld: await fetch("./togocoord_location_encoder", options).then(r => r.json()), correspondence: new_json};
  });
  return Promise.all(promises);
}
```