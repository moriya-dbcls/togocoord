# TogoCoord - decode 'Location ID' code to JSON

## Description

* Remote-Reference にも一応対応

## Parameters

* `location_id`
  * default: P12883:complement(join(<12..78,134..202))
  * example: 467, 340..565, <345..500, 1..>888, 102.110, 123^124, join (order), complement

## `parse`
```javascript
({location_id}) => {
  const make_range = (code) => {
    let type = "";
    let complement = false;
    let reference = false;
    let before = false;
    let beyond = false;
    if (code.match(/^complement\(.+\)$/)) {
      code = code.match(/^complement\((.+)\)$/)[1];
      complement = true;
    }
    if (code.match(":")) {
      [, reference, code] = code.match(/^([^:]+):(.+)/);
    }
    if (code.match(/^-*\d+$/)) {
      return {
        type: "position",
        ...(complement && { complement: complement}),
        begin: parseInt(code),
        end: parseInt(code)
      };
    } else if (code.match(/-*\d+[^\d]+-*\d+/)) {
      const s = code.match(/(-*\d+)[^\d]+(-*\d+)/);
      const b = parseInt(s[1]);
      const e = parseInt(s[2]);
      let type = "range";
      if (code.match(/-*\d+\.\.-*\d+/)) {
        if (code.match(/</)) before = true;
        if (code.match(/>/)) beyond = true;
      } else if (code.match(/-*\d+\.-*\d+/)) { 
        type = "uncertain";
      } else if (code.match(/-*\d+\^-*\d+/)) {
        if (e - b == 1) type = "between";
        else type = "between-error";
      }
      return {
        ...(reference && {reference: reference}),
        type: type,
        ...(before && {before: before}),
        ...(beyond && {betond: beyond}),
        ...(complement && {complement: complement}),
        begin: b,
        end: e
      };
    }
  }
  
  let assembl = false;
  let reference = false;
  let code = location_id.replace(/\s/g, '');
  if (code.match(":")) {
    [, reference, code] = code.match(/^([^:]+):(.+)$/);
    if (reference.match(/-/)) {
      [, assembl, reference] = reference.match(/^(.+)-(.+)$/);
    }
  }
  let complement = false;
  let ranges = [];
  if (code.match(/^complement\(.+\)$/)) {
    code = code.match(/^complement\((.+)\)$/)[1];
    complement = true;
  }
  if (code.match(/^join/) || code.match(/^order/)) {
    code.match(/^\w+\((.+)\)$/)[1].split(/,/).forEach(d => {
      ranges.push(make_range(d));
    });
  } else {
    ranges.push(make_range(code));
  }
  
  return {
    ...(assembl && {assembl: assembl}),
    reference: reference,
    ...(complement && {complement: complement}),
    ranges: ranges
  }
}
```