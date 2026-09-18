# TogoCoord - convert ensembl transcript to protein

## Description

* req.
  * [./togocoord_location_decoder](./togocoord_location_decoder)
  * [./togocoord_location_encoder](./togocoord_location_encoder)
  * [./togocoord_coord_converter](./togocoord_coord_converter)

## Parameters

* `sources`
  * default: ["ENST00000252519"]
  * example: ["ENSMUST00000073973"] ["ENST00000252519","ENST00000427411","ENST00000678073","ENST00000680121"] ["ENST00000531224"]
* `location`
  * default: join(233..235,5958..6116,8639..8641)
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example: ENSP00000389326
* `mode`
  * default: transcript
  * example: mrna, transcript

## `target_align_ranges`
```javascript
async ({sources, target, mode}) => {
  const transcripts = JSON.parse(sources);

  const make_range = (reference, begin, end, complement, pre_ranges, cds_begin) => {
    let delta = false;
    let original_end = end;
    if (!complement) {
      if (!pre_ranges[0]) {
        delta = cds_begin - 1;
        end = end - begin + 1 - cds_begin + 1;
        begin = 1 - cds_begin;
      } else {
        delta = pre_ranges.at(-1).delta + (begin - pre_ranges.at(-1).original_end) - 1;
        end = end - begin + 1 + pre_ranges.at(-1).end;
        begin = pre_ranges.at(-1).end + 1;
      }
    } else {
      console.log(false);
    }
    return {
      reference: reference,
      ...(complement && {joined_complement: complement}),
      mode: "nt2aa",
      begin: begin,
      end: end,
      delta: delta,
      original_end: original_end
    };
  }
  
  let ranges = [];                            
  for (let accession of transcripts) {
    const api = "https://rest.ensembl.org/lookup/id/" + accession + "?expand=1&content-type=application/json";
    console.log(api);
    const json = await fetch(api).then(r => r.json());
    
    let range_list = [];
    let reference = false;
    let complement = false;
    if (json.strand == -1) complement = true;
  	let genome_begin = json.Translation.start;
 	let genome_end = json.Translation.end;
    
    reference = json.Translation.id;
    if (target && target != reference) continue;
    let cds_begin = 0;
    let cds_end = 0;
    let begin_flag = true;
    if (!complement) {
      for (let d of json.Exon) {
        if (begin_flag) {
          if (d.end < genome_begin) {
            cds_begin += d.end - d.start + 1;
          } else if (genome_begin <= d.end) {
            cds_begin += genome_begin - d.start + 1;
            begin_flag = false;
          }
        }
        if (d.end < genome_end) {
          cds_end += d.end - d.start + 1;
        } else if (genome_end <= d.end) {
          cds_end += genome_end - d.start + 1;
          break;
        }
      }
    } else {
       for (let d of json.Exon) {
        if (begin_flag) {
          if (genome_end < d.start && begin_flag) {
            cds_begin += d.end - d.start + 1;
          } else if (d.start <= genome_end) {
            cds_begin += d.end - genome_end + 1;
            begin_flag = false;
          }
        }
        if (genome_begin < d.start) {
          cds_end += d.end - d.start + 1;
        } else if (d.start <= genome_begin) {
          cds_end += d.end - genome_begin + 1;
          break;
        }
      }       
    }
    if (mode == "mrna") {
      range_list.push({begin: 1, end: cds_end});
      complement = false;
    } else {
      let offset = json.Translation.start - 1;
      if (complement) {
        offset = json.Translation.end + 1;
      }

      for (let d of json.Exon) {
        if (complement && d.start <= genome_begin && genome_begin <= d.end) {
          range_list.push({begin: genome_begin - offset, end: d.end - offset});
          break;
        } else if (!complement && d.start <= genome_end && genome_end <= d.end) {
          range_list.push({begin: d.start - offset, end: genome_end - offset});
          break;
        } else {
          range_list.push({begin: d.start - offset, end: d.end - offset});
        }
      }
      if (complement) {
        complement = false;
        range_list.forEach((d, i) => {
          let tmp = d.begin;
          range_list[i].begin = d.end * (-1);
          range_list[i].end = tmp * (-1);
        });
      }
    }
   //console.log(range_list); 
    
    let pre_ranges = [];
    for (let d of range_list) {
      pre_ranges.push(make_range(reference, d.begin, d.end, complement, pre_ranges, cds_begin));
    }
        
    if (!complement) {
      pre_ranges.at(-1).end -= 3; // not align stop codon in AAseq
      ranges.push(...pre_ranges);
    } else {
      pre_ranges.at(-1).begin += 3; // not align stop codon in AAseq
      ranges.push(...pre_ranges.reverse());
    }
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
    let target_prefix = "ensembl";
    options.body = "data=" + encodeURIComponent(JSON.stringify(new_json)) + "&reference_prefix=" + target_prefix;
    console.log(JSON.stringify(new_json));
    return {jsonld: await fetch("./togocoord_location_encoder", options).then(r => r.json()), correspondence: new_json};
  });
  return Promise.all(promises);
}
```