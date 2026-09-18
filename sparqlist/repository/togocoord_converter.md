# TogoCoord - converter

## Description

* req.
  * [./togocoord_location_decoder](./togocoord_location_decoder)
  * [./togocoord_location_encoder](./togocoord_location_encoder)
  * [./togocoord_coord_transformer](./togocoord_coord_transformer)

## Parameters

* `source`
  * default: P12883
* `location`
  * default: join(30..300,400..960)
  * example: 46, 100..230, <105..230, 1..>138, 102.110, 123^124, join(12..78,134..202)
* `target` (Opt.)
  * example: 4P7H
* `api`
  * default: togocoord_convert_uniprot_pdb
* `target_prefix`
  * default: pdb

## `return`
```javascript
async ({source, location, target, api, target_prefix}) => {
  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  
  // decode location to coords
  options.body = "location_id=" + encodeURIComponent(source + ":" + location);
  const sourceInput = await fetch("./togocoord_location_decoder", options).then(r => r.json());
  // get mappings
  options.body = "source=" + encodeURIComponent(source);
  const mappings = await fetch("./" + api, options).then(r => r.json());
  if (mappings.length === 0) return [];

  const result = [];
  await Promise.all(
    mappings.map(async (mapping) => {
     if (target && target != mapping.target.reference.replace(/-.+$/, "")) return; 
     // convert coords
   	 options.body = "sourceInput=" + encodeURIComponent(JSON.stringify(sourceInput))
    	    + "&mapping=" + encodeURIComponent(JSON.stringify(mapping));
     const targetRanges = await fetch("./togocoord_coord_transformer", options).then(r => r.json());
  	 if (targetRanges.ranges.length === 0) return;

   	 const targetObj = {
   	   ...(sourceInput.assembl && {assembl: sourceInput.assembl}),
    	  reference: mapping.target.reference,
    	  ...(sourceInput.complement && {complement: sourceInput.complement}),
    	  ranges: targetRanges.ranges.map(r => {return {begin: r.target.begin, end: r.target.end, before: r.before, beyond: r.beyond, type: r.type};})
   	 };
     //onsole.log(targetObj);
   	 options.body = "data=" + encodeURIComponent(JSON.stringify(targetObj)) + "&reference_prefix=" + target_prefix;
   	 //console.log(options.body);
   	 result.push({jsonld: await fetch("./togocoord_location_encoder", options).then(r => r.json()), mapping: targetRanges});
    })
  );
  return result;
}
```