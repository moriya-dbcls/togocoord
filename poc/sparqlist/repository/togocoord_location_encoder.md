# TogoCoord - encode location JSON-LD from pre-JSON

## Description

* JSON から JSON-LD を構築
* Remote-Reference にも一応対応

## Parameters

* `data` (position, range (before, byond, before-beyond, uncertain, between), join, complement)
  * default: {"reference":"2AJF-A","join":true,"ranges":[{"type":"range","beyond":true,"begin":30,"end":210},{"type":"range","before":true,"beyond":true,"begin":211,"end":300},{"type":"range","beyond":true,"begin":400,"end":403},{"type":"range","before":true,"beyond":true,"begin":411,"end":624},{"type":"ramge","before":true,"beyond":true,"begin":645,"end":731},{"type":"range","before":true,"beyond":true,"begin":736,"end":787}]}
* `reference_prefix`
  * default: pdb
* `concat`  1: 連続する隣接 range を連結
  * default: 1

## `parse`
```javascript
({data, reference_prefix, concat}) => {
  const json = JSON.parse(data);

  // concatenate consecutive adjacent ranges strictly
  let concatenated_ranges = [];
  const concat_range = (a, i) => {
    let b = json.ranges[i];
    if (b 
        && ((b.begin - a.end == 1 
             && (!a.before_codon_position || (a.before_codon_position == 3 && b.end_codon_position == 1)))
         || (a.end_codon_position && b.begin == a.end && b.begin_codon_position - a.end_codon_position == 1))
        && !a.beyond && !b.before
        && !a.complement && !b.complement
        && a.reference === b.reference) {
      concat_range({
        type: "range",
        ...(a.reference && {reference: a.reference}),
        ...(a.before && {before: a.before}),
        ...(b.beyond && {beyond: b.beyond}),
        begin: a.begin,
        end: b.end,
        ...(a.before_codon_position && {before_codon_position: a.before_codon_position}),
        ...(b.end_codon_position && {end_codon_position: b.end_codon_position}),
      }, i + 1);
    } else {
      concatenated_ranges.push(a);
      if (b) concat_range(b, i + 1);
    }
  };
  if (concat == "1") concat_range(json.ranges[0], 1);
  else concatenated_ranges = json.ranges;
  console.log(concatenated_ranges);
        
  // construct location ID and location JSON
  let ranges = [];
  let locations = [];
  concatenated_ranges.forEach((d, i) => {
    let reference = json.reference;
    let type = d.type;
    let complement = d.complement;
    let begin = d.begin;
    let end = d.end;
    let r = ""; // location ID
    let l = {}; // location JSON
    if (complement) {
      r += "complement(";
    }
    if (type == "position" && begin == end) {
      l = {
        type: "ExactPosition",
        position: begin,
        reference: reference_prefix + ":" + reference
      };
      r += begin;
      if (d.begin_codon_position) r += "c" + d.begin_codon_position;
    } else {
      if (d.reference) {
        reference = d.reference;
        r += reference + ":";
      }
      l = {
        type: "Region",
        begin: {
          type: "ExactPosition",
          position: begin,
          reference: reference_prefix + ":" + reference  
        },
        end: {
          type: "ExactPosition",
          position: end,
          reference: reference_prefix + ":" + reference  
        }
      }
      if (d.before) { // && i == 0) {
        l.begin.type = "FuzzyPosition";
        r += "<";
      }
      r += begin;
      if (d.begin_codon_position) r += "c" + d.begin_codon_position;
      if (type.match("uncertain")) {
        l.type = "InRangePosition";
        l.begin.type = "Position";
        l.end.type = "Position";
        r += ".";
      } else if (type.match("between")) {
        l.type = "InBetweenPosition";
        l.after = l.begin;
        l.before = l.end;
        delete(l.after.type);
        delete(l.before.type);
        delete(l.begin);
        delete(l.end);
        r += "^";  
      } else r += "..";
      if (d.beyond) { // && i == concatenated_ranges.length - 1) {
        l.end.type = "FuzzyPosition";
        r += ">";
      }
      r += end;
      if (d.end_codon_position) r += "c" + d.end_codon_position;
    }
    if (complement) {
      l.strand = "NegativeStrand";
      r += ")";
    }
    ranges.push(r);
    locations.push(l);
  });
  
  let res = {
    "@context": "http://example.org/context/faldo.jsonld",
    id: ""
  };
  let location_code;
  if (concatenated_ranges.length > 1) {
    location_code = "join(" + ranges.join(",") + ")";
    res.location = {
      type: "ListOfRegions",
      member: []
    };
    locations.forEach((l, i) => {
      l.order = i + 1;
      res.location.member.push(l);
    }); 
  } else {
    location_code = ranges[0];
    res.location = locations[0];
  }
  if (json.complement) {
    location_code = "complement(" + location_code + ")";
    res.location.strand = "NegativeStrand";
  }
  let assembl_reference = json.reference;
  if (json.assembl) assembl_reference = json.assembl + "-" + json.reference;
  res.id = "http://example.org/" + assembl_reference + ":" + location_code
  return res;                  
}
```