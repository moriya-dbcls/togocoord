# TogoCoord - transform coords from source to target

## Description

* 座標変換

## Parameters

* `sourceInput`
  * default: {"reference":"P12883","ranges":[{"type":"range","begin":30,"end":300},{"type":"range","begin":400,"end":960}]}
* `mapping`
  * default: {"source":{"reference":"P12883","type":"aa"},"target":{"reference":"4P7H-A","type":"aa","strand":"+"},"ranges":[{"source":{"begin":3,"end":203},"target":{"begin":3,"end":203}},{"source":{"begin":211,"end":368},"target":{"begin":211,"end":368}},{"source":{"begin":371,"end":403},"target":{"begin":371,"end":403}},{"source":{"begin":411,"end":624},"target":{"begin":411,"end":624}},{"source":{"begin":645,"end":731},"target":{"begin":645,"end":731}},{"source":{"begin":736,"end":787},"target":{"begin":736,"end":787}}]}
  
## `return`
 * transform mapping results through a chain, handling input complement and mRNA/nt/aa conversions
 * @param {Object} mapping
 *   mapping.source: { reference: string, type: 'nt'|'aa' }
 *   mapping.target: { reference: string, type: 'nt'|'aa', strand: '+'|'-' }
 *   mapping.ranges: Array<{ source:{begin:number,end:number}, target:{begin:number,end:number} }>
 * @param {Object} input  // sourceRanges wrapper
 *   { complement?: boolean, ranges: Array<{begin:number,end:number,before?:boolean,beyond?:boolean}> }
 * @returns {Object} result
 *   { reference: string, complement: boolean,
 *     ranges: Array<{
 *       source:{begin:number,end:number},
 *       target:{begin:number,end:number},
 *       before: boolean,  // gap before source or input before
 *       beyond: boolean,  // gap after source or input beyond
 *       codonIndex?:{begin:number,end:number}
 *   }> }
```javascript
({sourceInput, mapping}) => {
  const { source, target, ranges: mapRanges } = JSON.parse(mapping);
  const negStrand = target.strand === '-';
  // determine mapping.source coordinate span
  const sortedMap = mapRanges.slice().sort((a, b) => a.source.begin - b.source.begin);
  const srcMin = sortedMap[0].source.begin;
  const srcMax = sortedMap[sortedMap.length - 1].source.end;
  
  const result = { reference: target.reference, complement: false, ranges: [] };
  const input = JSON.parse(sourceInput);
  const inCompRoot = Boolean(input.complement);
  for (const seg of input.ranges) {
    let { begin: qBeg, end: qEnd, before: inBefore = false, beyond: inBeyond = false } = seg;
    // handle input complement: invert coordinates
    if (inCompRoot) {
      const newBeg = srcMin + (srcMax - qEnd);
      const newEnd = srcMin + (srcMax - qBeg);
      qBeg = newBeg;
      qEnd = newEnd;
      // swap input flags
      [inBefore, inBeyond] = [inBeyond, inBefore];
    }

    // collect mapped source->target intervals
    const mappedInts = [];
    for (const m of sortedMap) {
      const s = m.source;
      const t = m.target;
      const oBeg = Math.max(qBeg, s.begin);
      const oEnd = Math.min(qEnd,   s.end);
      if (oBeg > oEnd) continue;
      // offset in source units
      let offBeg, offEnd;
      if (source.type === 'nt' && target.type === 'aa') {
        offBeg = Math.floor((oBeg - s.begin) / 3);
        offEnd = Math.floor((oEnd - s.begin) / 3);
      } else if (source.type === 'aa' && target.type === 'nt') {
        offBeg = (oBeg - s.begin) * 3;
        offEnd = (oEnd - s.begin) * 3;
      } else {
        offBeg = oBeg - s.begin;
        offEnd = oEnd - s.begin;
      }
      // compute target coords
      const p1 = negStrand ? t.end - offBeg : t.begin + offBeg;
      const p2 = negStrand ? t.end - offEnd : t.begin + offEnd;
      const tgtBeg = Math.min(p1, p2);
      const tgtEnd = Math.max(p1, p2);
      // compute codon indices if nt->aa
      let cBeg, cEnd;
      if (source.type === 'nt' && target.type === 'aa') {
        cBeg = ((oBeg - s.begin) % 3 + 3) % 3 + 1;
        cEnd = ((oEnd - s.begin) % 3 + 3) % 3 + 1;
      }
      let obj = { source: { begin: oBeg, end: oEnd }, target: { begin: tgtBeg, end: tgtEnd }, type: seg.type};
      if (cBeg) obj.targetCodonIndex = { begin: cBeg, end: cEnd };
      mappedInts.push(obj);
    }
    // sort and merge with gaps
    mappedInts.sort((a, b) => a.source.begin - b.source.begin);
    let cursor = qBeg;
    const segs = [];
    for (const mi of mappedInts) {
      if (mi.source.begin > cursor) {
        segs.push({ mapped: false, source: { begin: cursor, end: mi.source.begin - 1 }, target: false });
      }
      let obj = { mapped: true, source: mi.source, target: mi.target, type: mi.type};
      if (mi.targetCodonIndex) obj.targetCodonIndex = mi.targetCodonIndex;
      segs.push(obj);
      cursor = mi.source.end + 1;
    }
    if (cursor <= qEnd) segs.push({ mapped: false, source: { begin: cursor, end: qEnd }, target: false });
    
    // set result.complement first iteration
    result.complement = negStrand !== inCompRoot;
    // emit mapped entries
    for (let i = 0; i < segs.length; i++) {
      const r = segs[i];
      if (!r.mapped) continue;
      const prev = segs[i - 1];
      const next = segs[i + 1];
      const before = (inBefore && i === 0) || (prev && !prev.mapped);
      const beyond = (inBeyond && i === segs.length - 1) || (next && !next.mapped);
      const entry = { source: r.source, target: r.target, before, beyond, type: r.type };
      if (r.targetCodonIndex !== undefined) entry.targetCodonIndex = r.codonIndex;
      result.ranges.push(entry);
    }
  }
  return result;
}
```