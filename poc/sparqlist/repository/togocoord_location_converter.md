# TogoCoord - convert location ID to JSON-LD

## Description

* req.
  * [./togocoord_location_decoder](./togocoord_location_decoder)
  * [./togocoord_location_encoder](./togocoord_location_encoder)

## Parameters

* `location_id`
  * default: GCA000000000-J00000:join(1..100,J00194.1:100..202)
  * example: 467, 340..565, <345..500, 1..>888, 102.110, 123^124, join (order), complement
* `reference_prefix`
  * default: undefined

## `return`
```javascript
async ({location_id, reference_prefix}) => {
  let options = {
    method: 'POST',
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/x-www-form-urlencoded'
    }
  }
  options.body = "location_id=" + encodeURIComponent(location_id);
  const json = await fetch("./togocoord_location_decoder", options).then(r => r.json());
  options.body = "data=" + encodeURIComponent(JSON.stringify(json)) + "&reference_prefix=" + reference_prefix;
  return await fetch("./togocoord_location_encoder", options).then(r => r.json());
}
```