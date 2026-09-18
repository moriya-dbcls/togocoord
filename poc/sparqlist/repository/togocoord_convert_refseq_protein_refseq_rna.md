# TogoCoord - convert RefSeq protein (isndc) to RefSeq mRNA (isndc) & insdc contig 

## Description

* req.
  * [./togocoord_location_decoder](./togocoord_location_decoder)
  * [./togocoord_location_encoder](./togocoord_location_encoder)
  * [./togocoord_coord_converter](./togocoord_coord_converter)

## Parameters

* `sources`
  * default: ["AAF99721","AAO25651"]
  * example: ["AAF64304"] ["NP_055554","NP_001373629","NP_001373630"] ["EAW47950"]
* `location`
  * default: 1..920
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example: NM_014739
* `mode`
  * example: mrna

## `target_align_ranges`
```javascript
async ({sources, target, mode}) => {
  let proteins = JSON.parse(sources);

  const make_range = (reference, begin, end, complement, pre_ranges) => {
    let delta = false;
    if (!complement) {
      if (!pre_ranges[0]) delta = 1 - begin;
      else delta = pre_ranges.at(-1).delta - (begin - pre_ranges.at(-1).end - 1);  // gap
    } else {
      if (!pre_ranges[0]) delta = (1 + end) * (-1);
      else delta = pre_ranges.at(-1).delta + (pre_ranges.at(-1).begin - end - 1);  // gap      
    }
    return {
      reference: reference,
      ...(complement && {joined_complement: complement}),
      mode: "aa2nt",
      begin: begin,
      end: end,
      delta: delta
    };
  }

  // UniParc で EMBL タイプでも mRNA じゃないことがあるので、怪しいものを個別チェック
  const check_mrna = async (code) => {
    if (code.match(/join/) || code.match(/order/) || parseInt(code.match(/^[^\d]*(\d+)\.\./)) > 1000) {
      const reference = code.match(/([^:\(]+):-*\d+\.\./)[1].replace(/\.\d+$/, "");
      const api = "https://togows.org/entry/ncbi-nucleotide/" + reference + ".json";
      try {
        const r = await fetch(api);
        const json = await r.json();
        if (json[0].moltype != "mRNA") return false;
        else return true;
      } catch (error) {
        console.log(error);
      }
    }
    else return true;
  }

  let ranges = [];      
  let chk_reference = {};
  for (let accession of proteins) {
    const api = "https://togows.org/entry/ncbi-protein/" + accession + ".json";
    const json = await fetch(api).then(r => r.json());
    
    for (const data of json) {  // await のため for
      for (const feature of data.features) {
        if (feature.coded_by) {
          for (let code of feature.coded_by) {
            //console.log(code);
            code = code.replace(/\s/g,"");
            if (mode == "mrna") {
              const mrna_flag = await check_mrna(code); 
              if (!mrna_flag) continue;
            }
            if (target && !code.match(target)) continue;
            let reference = false;
            let complement = false;
            let join = false;
            let pre_ranges = [];
            let begin;
            let end;
            if (code.match(/^complement\(.+\)$/)) {
              code = code.match(/^complement\((.+)\)$/)[1];
              complement = true;
            }
            if (code.match(/^join/) || code.match(/^order/)) {
              if (!complement) {
                code.match(/^\w+\((.+)\)$/)[1].split(/,/).forEach(d => {
                  [reference, begin, end] = d.match(/^([^:]+):(-*\d+)\.\.(-*\d+)/).slice(1,4);
                  pre_ranges.push(make_range(reference.split(/\./)[0], parseInt(begin), parseInt(end), complement, pre_ranges));
                });
              } else {
                code.match(/^\w+\((.+)\)$/)[1].split(/,/).reverse().forEach(d => {
                  [reference, begin, end] = d.match(/^([^:]+):(-*\d+)\.\.(-*\d+)/).slice(1,4);
                  pre_ranges.push(make_range(reference.split(/\./)[0], parseInt(begin), parseInt(end), complement, pre_ranges));
                });
              }
            } else {
              [reference, begin, end] = code.match(/^([^:]+):(-*\d+)\.\.(-*\d+)/).slice(1,4);
              if (!complement) pre_ranges.push(make_range(reference.split(/\./)[0], parseInt(begin), parseInt(end), complement, pre_ranges));
              else pre_ranges.push(make_range(reference.split(/\./)[0], parseInt(begin), parseInt(end), complement, pre_ranges));
            }

            if (chk_reference[reference]) continue;
            if (!complement) {
              pre_ranges.at(-1).end -= 3; // not align stop codon in AAseq
              ranges.push(...pre_ranges);
            }
            else {
              pre_ranges.at(-1).begin += 3; // not align stop codon in AAseq
           	  ranges.push(...pre_ranges.reverse());
            }
            chk_reference[reference] = true;
          }
        }
      }
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
    let target_prefix = "insdc";
    if (reference.match(/^(?:NM|NR|XM|XR)_\d+/)) target_prefix = "refseq";
    options.body = "data=" + encodeURIComponent(JSON.stringify(new_json)) + "&reference_prefix=" + target_prefix;
    //console.log(options.body);
    return {jsonld: await fetch("./togocoord_location_encoder", options).then(r => r.json()), correspondence: new_json};
  });
  return Promise.all(promises);
}
```