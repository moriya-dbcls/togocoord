# TogoCoord - convert UniProt to PDB

## Description

* 25/07/09

## Parameters

* `source`
  * default: P12883

## Endpoint

https://rdfportal.org/pdb/sparql

## `pdb`
* UniProt と PDB の align range と座標のズレ
* chian ID
```sparql
PREFIX pdbo: <http://rdf.wwpdb.org/schema/pdbx-v50.owl#>
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX up: <http://purl.uniprot.org/uniprot/>
SELECT DISTINCT ?pdb ?db_align_begin ?auth_align_begin ?auth_align_end (GROUP_CONCAT(DISTINCT ?chain_id ;  separator=',') AS ?chains)
WHERE {
  VALUES ?uniprot { up:{{source}} }
  {{#if target}}
  VALUES ?pdb { "{{target}}" }
  {{/if}}
  [] dct:identifier ?pdb ;
     pdbo:has_pdbx_audit_revision_historyCategory/pdbo:has_pdbx_audit_revision_history [
       pdbo:pdbx_audit_revision_history.ordinal "1" ;
       pdbo:pdbx_audit_revision_history.revision_date ?date 
     ] ;
     pdbo:has_entityCategory/pdbo:has_entity [
       pdbo:referenced_by_struct_ref [
         pdbo:link_to_uniprot ?uniprot ; # UniProt
         pdbo:referenced_by_struct_ref_seq [
           pdbo:struct_ref_seq.db_align_beg ?db_align_begin ;              # 配列 align の UniProt の begin
           pdbo:struct_ref_seq.pdbx_auth_seq_align_beg ?auth_align_begin ; # 配列 align の PDB chain の begin
           pdbo:struct_ref_seq.pdbx_auth_seq_align_end ?auth_align_end ; # 配列 align の PDB chain の end
         ]
       ] ;
       pdbo:referenced_by_struct_asym/pdbo:struct_asym.id ?chain_id
     ] .
  #OPTIONAL {
  #  ?entry pdbo:has_refineCategory/pdbo:has_refine ?refine .
  #  OPTIONAL {
  #    ?refine pdbo:refine.ls_d_res_high ?resolution_high .
  #  }
  #  OPTIONAL {
  #    ?refine pdbo:refine.ls_R_factor_R_free ?rfree .
  #  }
  #  OPTIONAL {
  #    ?refine pdbo:refine.ls_R_factor_R_work ?rwork .
  #  }
  #}
}
ORDER BY ?pdb ?auth_align_begin
```

## `pdb_3d_pos`
* 3D 座標がある position
  * 3D データ不連続性取得のための
```sparql
PREFIX pdbo: <http://rdf.wwpdb.org/schema/pdbx-v50.owl#>
PREFIX dct: <http://purl.org/dc/terms/>
PREFIX up: <http://purl.uniprot.org/uniprot/>
SELECT DISTINCT ?pdb ?pos
WHERE {
  VALUES ?uniprot { up:{{source}} }
  {{#if target}}
  VALUES ?pdb { "{{target}}" }
  {{/if}}
  [] dct:identifier ?pdb ;
     pdbo:has_entityCategory/pdbo:has_entity [
       pdbo:referenced_by_struct_ref/pdbo:link_to_uniprot ?uniprot ;
       pdbo:referenced_by_entity_poly/pdbo:referenced_by_entity_poly_seq/pdbo:referenced_by_pdbx_poly_seq_scheme/pdbo:pdbx_poly_seq_scheme.auth_seq_num ?auth_pos
     ] .
  BIND (xsd:integer(?auth_pos) AS ?pos)
}
ORDER BY ?pdb ?pos
```

## `mappings`
* 不連続な align range に変換
```javascript
({source, pdb, pdb_3d_pos}) => {
  if (!pdb.results.bindings[0] || !pdb_3d_pos.results.bindings[0]) return [];

  let pdb2begin = {};
  let pdb2end   = {};
  let pdb2delta = {};
  let pdb2chain = {};
  pdb.results.bindings.forEach(d => {
    const pdb_id = d.pdb.value;
    console.log(d);
    pdb2begin[pdb_id] = parseInt(d.auth_align_begin.value);
    pdb2end[pdb_id]   = parseInt(d.auth_align_end.value);
    pdb2delta[pdb_id] = parseInt(d.db_align_begin.value) - parseInt(d.auth_align_begin.value);
    pdb2chain[pdb_id] = d.chains.value.split(/,/);
  });

  let res = {};
  let pre_id = false;
  let pre_pos = false;
  let begin = false;
  const addRes = (id, begin, end) => {
    for (let chain of pdb2chain[id]) {
      const reference = id + "-" + chain;
      if (!res[reference]) {
        res[reference] = {
          source: {
            reference: source,
            type: "aa"
          },
          target: {
            reference: reference,
            type: "aa",
            strand: "+"
          },
          ranges: []
        };
      }
      res[reference].ranges.push({
        source: {begin: begin + pdb2delta[pre_id], end: end + pdb2delta[pre_id]},
        target: {begin: begin, end: end}
      });
    }
  };
  
  for (let d of  pdb_3d_pos.results.bindings) {
    let pdb_id = d.pdb.value;
    let pos = parseInt(d.pos.value);
    if (pos < pdb2begin[pdb_id] || pdb2end[pdb_id] < pos) continue;
    if (pre_id != pdb_id || pos - pre_pos != 1) {  // align が不連続な場合
      if (pre_id) {
        addRes(pre_id, begin, pre_pos);
        begin = pos;
      }
      if (pre_id != pdb_id) pre_id = pdb_id;
    }
    if (begin === false) begin = pos;
    pre_pos = pos;
  };
  addRes(pre_id, begin, pre_pos);

  return Object.values(res);
}
```

