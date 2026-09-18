import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { createContext, formatLocationId, mapLocation, parseLocationId, type CoordContext } from "@togocoord/core";
import { ingestGenBank } from "../src/adapter-gbff.ts";
import { ingestGff3 } from "../src/adapter-gff3.ts";
import { edgeMapping, type Edge, type IngestResult } from "../src/model.ts";
import { extract, MemorySequenceSource, translate, type SequenceSource } from "../src/sequence.ts";
import { validateTranscript } from "../src/validate.ts";
import { parseFasta, parseFastaHeaders } from "../src/fasta.ts";
import { assemblyReportAliases, assemblyReportInfo, assemblyReportSeqids, assemblyReportSequences } from "../src/common.ts";
import { chr3Source, fixture, genbankSource } from "./helpers.ts";

const GENOMES = ["NC_045512.2", "NC_012920.1", "NC_001405.1"];
const NUCLEOTIDE_GB = [...GENOMES.map((g) => `${g}.gb`), "NC_002127.1.gb", "NM_000581.4.gb", "NM_014739.3.gb"];

const gb = new Map(NUCLEOTIDE_GB.map((name) => [name, ingestGenBank(fixture(name), { file: name })]));
const cdsEdges = (r: IngestResult) => r.edges.filter((e) => e.attributes.codonStart !== undefined);
const edge = (r: IngestResult, from: string) => r.edges.find((e) => e.from === from)!;

function ctxFor(...results: IngestResult[]): CoordContext {
  const units = new Map(results.flatMap((r) => r.sequences.map((s) => [s.ref, s.unit] as const)));
  return createContext({ units: (ref) => units.get(ref) });
}

describe("GenBank adapter: CDS self-validation on real records", () => {
  it("every CDS translation is reproduced from the genome (68 CDS, 6 nucleotide records)", () => {
    const all = [...gb.values()].flatMap(cdsEdges);
    assert.equal(all.length, 68);
    const failed = all.filter((e) => e.validation.status !== "ok").map((e) => `${e.from}: ${e.validation.detail}`);
    assert.deepEqual(failed, []);
  });

  it("handles alternative starts, selenocysteine and incomplete stop codons", () => {
    assert.match(edge(gb.get("NC_012920.1.gb")!, "refseq:YP_003024027.1").validation.detail!, /alternative start codon ATT/);
    assert.match(edge(gb.get("NM_000581.4.gb")!, "refseq:NP_000572.2").validation.detail!, /transl_except Sec at residue 49/);
    const nd1 = edge(gb.get("NC_012920.1.gb")!, "refseq:YP_003024026.1");
    assert.equal(nd1.attributes.aaLength, "318");
    assert.equal(nd1.location, "refseq:NC_012920.1:3307..4262");
  });

  it("keeps the ribosomal slippage of SARS-CoV-2 ORF1ab", () => {
    const e = edge(gb.get("NC_045512.2.gb")!, "refseq:YP_009724389.1");
    assert.equal(e.location, "refseq:NC_045512.2:join(266..13468,13468..21555)");
    assert.equal(e.attributes.ribosomalSlippage, "true");
  });

  it("records topology, taxon and refget digests", () => {
    const mt = gb.get("NC_012920.1.gb")!.sequences.find((s) => s.ref === "refseq:NC_012920.1")!;
    assert.equal(mt.topology, "circular");
    assert.equal(mt.taxon, 9606);
    assert.match(mt.digest!, /^SQ\.[A-Za-z0-9_-]{32}$/);
    const dloop = gb.get("NC_012920.1.gb")!.annotations.find((a) => a.type === "D-loop")!;
    assert.equal(dloop.location, "refseq:NC_012920.1:complement(join(16024..16569,1..576))");
  });

  it("detects a wrong translation or genetic code (validation is not vacuous)", () => {
    const text = fixture("NM_000581.4.gb");
    const broken = text.replace('/translation="MCAAR', '/translation="MCAAW');
    assert.equal(ingestGenBank(broken).edges[0]!.validation.status, "mismatch");
    const noExcept = text.replace(/\s+\/transl_except=\(pos:220\.\.222,aa:Sec\)/, "");
    // Without /transl_except the UGA reads as a stop; the published U is accepted as a recoded stop and reported.
    const recoded = ingestGenBank(noExcept).edges[0]!.validation;
    assert.equal(recoded.status, "ok");
    assert.match(recoded.detail!, /selenocysteine at residue 49 \(recoded stop codon\)/);
    const wrongTable = fixture("NC_012920.1.gb").replace(/\/transl_table=2/g, "/transl_table=1");
    assert.ok(cdsEdges(ingestGenBank(wrongTable)).some((e) => e.validation.status === "mismatch"));
  });
});

describe("GenPept /coded_by", () => {
  it("gives the same protein -> mRNA mapping as the mRNA record, and validates with its sequence", () => {
    const gp = ingestGenBank(fixture("NP_055554.1.gp"), { source: genbankSource("NM_014739.3.gb") });
    const fromGp = edge(gp, "refseq:NP_055554.1");
    const fromNm = edge(gb.get("NM_014739.3.gb")!, "refseq:NP_055554.1");
    assert.equal(fromGp.location, "refseq:NM_014739.3:249..3011");
    assert.deepEqual(fromGp.blocks, fromNm.blocks);
    assert.equal(fromGp.validation.status, "ok");
  });
});

describe("GFF3 adapter agrees with the GenBank adapter", () => {
  for (const g of GENOMES) {
    it(`${g}: identical CDS mappings, inferred lengths equal /translation lengths, all CDS validate`, () => {
      const fromGb = gb.get(`${g}.gb`)!;
      const fromGff = ingestGff3(fixture(`${g}.gff3`), { source: genbankSource(`${g}.gb`) });
      assert.deepEqual(fromGff.warnings, []);
      const summary = (r: IngestResult) =>
        cdsEdges(r)
          .map((e) => ({ from: e.from, location: e.location, blocks: e.blocks, codonStart: e.attributes.codonStart, aaLength: e.attributes.aaLength }))
          .sort((a, b) => a.from.localeCompare(b.from));
      assert.deepEqual(summary(fromGff), summary(fromGb));
      assert.deepEqual(cdsEdges(fromGff).filter((e) => e.validation.status !== "ok"), []);
    });
  }

  it("wraps origin-spanning features of circular molecules", () => {
    const mt = ingestGff3(fixture("NC_012920.1.gff3"));
    assert.equal(mt.annotations.find((a) => a.type === "D_loop")!.location, "refseq:NC_012920.1:complement(join(16024..16569,1..576))");
    assert.equal(mt.sequences.find((s) => s.ref === "refseq:NC_012920.1")!.topology, "circular");
  });
});

describe("GFF3 cDNA_match (RefSeq transcript vs GRCh38)", () => {
  const result = ingestGff3(fixture("NC_000003.12_cDNA_match.gff3"), { source: chr3Source() });
  const ctx = ctxFor(result);
  const convert = (text: string, e: Edge, inverse = false) => {
    const m = edgeMapping(e);
    const r = mapLocation(parseLocationId(text, ctx), inverse ? m.inverse() : m, ctx);
    return { targets: r.targets.map((t) => formatLocationId(t.location, ctx)), unmapped: r.unmapped && formatLocationId(r.unmapped, ctx) };
  };
  const rybp = edge(result, "refseq:NM_012234.7");
  const tmem = edge(result, "refseq:NR_047574.1");

  it("builds alignment edges and reproduces the reported mismatch counts", () => {
    assert.deepEqual(result.warnings, []);
    assert.equal(result.edges.length, 3);
    assert.equal(rybp.kind, "alignment");
    assert.equal(rybp.validation.status, "ok");
    assert.equal(rybp.validation.detail, "4440/4445 aligned bases identical; num_mismatch=5");
    assert.equal(tmem.validation.detail, "1429/1429 aligned bases identical; num_mismatch=0");
    assert.equal(edge(result, "refseq:NR_047573.1").validation.status, "skipped");
  });

  it("minus strand: a base inserted in the transcript has no genomic counterpart", () => {
    assert.deepEqual(convert("refseq:NM_012234.7:1154", rybp), { targets: [], unmapped: "refseq:NM_012234.7:1154" });
    assert.deepEqual(convert("refseq:NM_012234.7:1153..1155", rybp), {
      targets: ["refseq:NC_000003.12:complement(72378101..72378102)"],
      unmapped: "refseq:NM_012234.7:1154",
    });
    assert.deepEqual(convert("refseq:NC_000003.12:complement(72378101..72378102)", rybp, true), {
      targets: ["refseq:NM_012234.7:join(1153,1155)"],
      unmapped: null,
    });
  });

  it("minus strand: the unaligned 5' end is truncated with a fuzzy boundary", () => {
    assert.deepEqual(convert("refseq:NM_012234.7:1..300", rybp), {
      targets: ["refseq:NC_000003.12:complement(72446537..>72446623)"],
      unmapped: "refseq:NM_012234.7:1..213",
    });
  });

  it("plus strand: a genomic base deleted from the transcript is unmapped", () => {
    assert.deepEqual(convert("refseq:NC_000003.12:194589594", tmem, true), { targets: [], unmapped: "refseq:NC_000003.12:194589594" });
    assert.deepEqual(convert("refseq:NC_000003.12:194589593..194589595", tmem, true), {
      targets: ["refseq:NR_047574.1:419..420"],
      unmapped: "refseq:NC_000003.12:194589594",
    });
    assert.deepEqual(convert("refseq:NR_047574.1:419^420", tmem), { targets: [], unmapped: "refseq:NR_047574.1:419^420" });
  });
});

describe("GFF3 transcripts are located by their exons (RYBP and GPX1, GRCh38)", () => {
  const source = new MemorySequenceSource();
  for (const name of ["NM_012234.7.fa", "NM_000581.4.fa"]) for (const [id, s] of parseFasta(fixture(name))) source.add(`refseq:${id}`, s);
  for (const [name, offset] of [["NC_000003.12_72374597-72446623.fa", 72374596], ["NC_000003.12_49357176-49358353.fa", 49357175]] as const) {
    for (const [, s] of parseFasta(fixture(name))) source.add("refseq:NC_000003.12", s, offset);
  }
  const result = ingestGff3(fixture("NC_000003.12_RYBP_GPX1.gff3"), { source });
  const transcript = (id: string) => result.edges.find((e) => e.from === id && e.attributes.type === "mRNA")!;

  it("builds the spliced model, not the gene-long mRNA row", () => {
    assert.equal(
      transcript("refseq:NM_012234.7").location,
      // exon 1 carries end_range=72446623,. (the transcript continues 5' of what GRCh38 contains)
      "refseq:NC_000003.12:complement(join(72374597..72378637,72379040..72379138,72379251..72379427,72446496..>72446623))",
    );
    assert.equal(result.edges.filter((e) => e.attributes.type === "mRNA").length, 6);
    assert.ok(result.annotations.every((a) => a.type !== "exon"));
    assert.match(result.annotations.find((a) => a.type === "mRNA" && a.attributes.ID?.[0] === "rna-NM_012234.7")!.location, /join\(/);
  });

  it("validates the model against the transcript sequence", () => {
    // RYBP: the RefSeq transcript has 5' sequence absent from GRCh38 and indels -> positions would shift.
    const rybp = transcript("refseq:NM_012234.7").validation;
    assert.deepEqual([rybp.status, rybp.basis], ["mismatch", "full"]);
    assert.match(rybp.detail!, /4445 nt, refseq:NM_012234.7 is 4662 nt/);
    const gpx1 = transcript("refseq:NM_000581.4").validation;
    assert.deepEqual([gpx1.status, gpx1.basis], ["ok", "full"]);
  });
});

describe("Ensembl GFF3 (GRCh38 release 116, GPX1): seqid map, versions, validation against Ensembl proteins", () => {
  const source = new MemorySequenceSource();
  for (const [header, s] of parseFastaHeaders(fixture("ensembl_GRCh38.116_GPX1.pep.fa"))) source.add(`ensembl:${header.split(" ")[0]}`, s);
  for (const [, s] of parseFasta(fixture("NC_000003.12_49357176-49358353.fa"))) source.add("refseq:NC_000003.12", s, 49357175);
  const seqids = assemblyReportSeqids(fixture("GRCh38.p14_assembly_report_chr3.txt"));
  const result = ingestGff3(fixture("ensembl_GRCh38.116_GPX1.gff3"), { source, seqidToRef: (id) => seqids.get(id) });
  const cds = result.edges.filter((e) => e.attributes.codonStart !== undefined);

  it("maps Ensembl seqid 3 to NC_000003.12 through the NCBI assembly report", () => {
    assert.deepEqual([seqids.get("3"), seqids.get("chr3"), seqids.get("CM000665.2")], ["refseq:NC_000003.12", "refseq:NC_000003.12", "refseq:NC_000003.12"]);
    assert.deepEqual(result.warnings, []);
  });

  it("keys proteins and transcripts with their Ensembl versions and validates every CDS against pep.all", () => {
    assert.equal(cds.length, 15);
    assert.ok(cds.every((e) => /^ensembl:ENSP\d{11}\.\d+$/.test(e.from) && e.to === "refseq:NC_000003.12"));
    assert.deepEqual(cds.filter((e) => e.validation.status !== "ok" || e.validation.basis !== "full").map((e) => `${e.from}: ${e.validation.detail}`), []);
    const mane = cds.find((e) => e.from === "ensembl:ENSP00000407375.1")!;
    assert.equal(mane.location, "refseq:NC_000003.12:complement(join(49357388..49357747,49358027..49358278))");
    assert.equal(mane.attributes.aaLengthSource, "protein sequence");
    const tx = result.edges.find((e) => e.from === "ensembl:ENST00000419783.3")!;
    assert.equal(tx.location, "refseq:NC_000003.12:complement(join(49357176..49357747,49358027..49358353))");
  });

  it("records checksums of the Ensembl protein: identical to UniProt P07203 and RefSeq NP_000572.2", () => {
    const seq = result.sequences.find((s) => s.ref === "ensembl:ENSP00000407375.1")!;
    assert.equal(seq.digest, "SQ.qsfE5UDDg5US4PKbR3-mgCvYbsfBsxFH");
  });
});

describe("Ensembl 5'-incomplete CDS: the protein starts with X for the incomplete first codon (ENSP00000349216.4)", () => {
  it("detects the leading X, maps residue 2 to the first complete codon and validates", () => {
    const source = new MemorySequenceSource();
    for (const [header, s] of parseFastaHeaders(fixture("ensembl_GRCh38.116_ENSP00000349216.pep.fa"))) source.add(`ensembl:${header.split(" ")[0]}`, s);
    for (const [, s] of parseFasta(fixture("NC_000001.11_930312-944575.fa"))) source.add("refseq:NC_000001.11", s, 930311);
    const seqids = assemblyReportSeqids(fixture("GRCh38.p14_assembly_report_chr1.txt"));
    const result = ingestGff3(fixture("ensembl_GRCh38.116_ENSP00000349216.gff3"), { source, seqidToRef: (id) => seqids.get(id) });
    const e = result.edges.find((x) => x.from === "ensembl:ENSP00000349216.4")!;
    assert.deepEqual([e.attributes.codonStart, e.attributes.leadingPartialCodon, e.validation.status], ["3", "true", "ok"]);
    assert.match(e.validation.detail!, /first residue is the incomplete first codon/);
    const ctx = createContext({ units: (r) => (r.startsWith("ensembl:ENSP") ? "aa" : undefined) });
    const ids = (t: string) => mapLocation(parseLocationId(t, ctx), edgeMapping(e), ctx).targets.map((x) => formatLocationId(x.location, ctx));
    assert.deepEqual(ids("ensembl:ENSP00000349216.4:1"), ["refseq:NC_000001.11:<930312..930313"]); // incomplete codon
    assert.deepEqual(ids("ensembl:ENSP00000349216.4:2"), ["refseq:NC_000001.11:930314..930316"]);
    assert.equal(translate(source.get("refseq:NC_000001.11", 930313, 930316)!), source.get("ensembl:ENSP00000349216.4", 1, 2));
  });
});

describe("trans-spliced CDS on both strands (Arabidopsis chloroplast rps12, NP_051038.1)", () => {
  it("keeps the GFF3 row order and reproduces the protein", () => {
    const source = new MemorySequenceSource();
    for (const [, s] of parseFasta(fixture("arabidopsis_NC_000932.1.fa"))) source.add("refseq:NC_000932.1", s);
    for (const [, s] of parseFasta(fixture("arabidopsis_NP_051038.1.faa"))) source.add("refseq:NP_051038.1", s);
    const result = ingestGff3(fixture("arabidopsis_NC_000932.1_rps12.gff3"), { source });
    assert.deepEqual(result.warnings, []);
    const e = edge(result, "refseq:NP_051038.1");
    assert.equal(e.location, "refseq:NC_000932.1:join(complement(69611..69724),139856..140087,140625..140650)");
    assert.deepEqual([e.validation.status, e.validation.basis], ["ok", "full"]);
  });
});

describe("transcript validation", () => {
  const ctx = createContext();
  const loc = parseLocationId("refseq:NC_000001.1:join(1..4,9..12)", ctx);
  const check = (transcript: string) =>
    validateTranscript(loc, "refseq:NM_000001.1", new MemorySequenceSource().add("refseq:NC_000001.1", "ACGTnnnnTTGC").add("refseq:NM_000001.1", transcript), ctx);
  it("accepts identity, sparse substitutions and a poly(A) tail; rejects indels", () => {
    assert.match(check("ACGTTTGC").detail!, /^identical/);
    assert.equal(check("ACGTTTGCAAAAAAAAAA").status, "ok");
    assert.match(check("ACGTTTGCAAAAAAAAAA").detail!, /10 nt poly\(A\) tail/);
    assert.equal(check("ACGTTTGCCGTACGTA").status, "mismatch");
    assert.equal(check("ACGTTTG").status, "mismatch");
  });
});

describe("GFF3 match with a reverse-complemented Target (RefSeqGene on the opposite strand)", () => {
  it("maps NG_044083.1 onto the minus strand and validates 304/304", () => {
    const source = new MemorySequenceSource();
    for (const [, s] of parseFasta(fixture("NG_044083.1_1-304.fa"))) source.add("refseq:NG_044083.1", s);
    for (const [, s] of parseFasta(fixture("NC_000003.12_282589-282892.fa"))) source.add("refseq:NC_000003.12", s, 282588);
    const result = ingestGff3(fixture("NC_000003.12_match_reverse_target.gff3"), { source });
    assert.deepEqual(result.warnings, []);
    const e = edge(result, "refseq:NG_044083.1");
    assert.match(e.validation.detail!, /^304\/304 aligned bases identical/);
    const ctx = ctxFor(result);
    const r = mapLocation(parseLocationId("refseq:NG_044083.1:1..10", ctx), edgeMapping(e), ctx);
    assert.deepEqual(r.targets.map((t) => formatLocationId(t.location, ctx)), ["refseq:NC_000003.12:complement(282883..282892)"]);
  });
});

describe("GFF3 gapped match with a reverse-complemented Target", () => {
  it("reads Gap operations along the genome ('+' row) and aligns NG_162589.1 without mismatches", () => {
    const source = new MemorySequenceSource();
    for (const [, s] of parseFasta(fixture("NG_162589.1_1-270.fa"))) source.add("refseq:NG_162589.1", s);
    for (const [, s] of parseFasta(fixture("NW_025791756.1_1186933-1187207.fa"))) source.add("refseq:NW_025791756.1", s, 1186932);
    const result = ingestGff3(fixture("NW_025791756.1_match_reverse_gapped.gff3"), { source });
    assert.deepEqual(result.warnings, []);
    const e = edge(result, "refseq:NG_162589.1");
    assert.equal(e.blocks.length, 2); // Gap=M217 D5 M53
    assert.match(e.validation.detail!, /^270\/270 aligned bases identical/);
  });
});

describe("end to end: every residue maps to a codon that encodes it", () => {
  const check = (results: IngestResult[], source: SequenceSource) => {
    const ctx = ctxFor(...results);
    const failures: string[] = [];
    let residues = 0;
    for (const e of results.flatMap(cdsEdges)) {
      const expected = source.get(e.from, 0, Number(e.attributes.aaLength));
      if (!expected) continue;
      const m = edgeMapping(e);
      const exceptions = new Set([...(e.validation.detail ?? "").matchAll(/at residue (\d+)/g)].map((x) => Number(x[1])));
      if (/alternative start/.test(e.validation.detail ?? "")) exceptions.add(1);
      for (let r = 1; r <= expected.length; r++) {
        residues++;
        const result = mapLocation(parseLocationId(`${e.from}:${r}`, ctx), m, ctx);
        const codon = result.targets.length === 1 ? extract(result.targets[0]!.location, source) : undefined;
        const aa = codon && translate(codon, Number(e.attributes.translTable));
        if (exceptions.has(r)) continue;
        if (aa !== expected[r - 1]) failures.push(`${e.from}:${r} -> ${codon} (${aa}) expected ${expected[r - 1]}`);
      }
    }
    return { failures, residues };
  };

  it("GenBank edges (all residues of all fixture proteins)", () => {
    const source = genbankSource(...NUCLEOTIDE_GB);
    // Expected residues come from /translation, read independently of the adapter.
    const proteins = genbankProteins();
    const { failures, residues } = check([...gb.values()], { get: (ref, s, t) => proteins.get(ref)?.slice(s, t) ?? source.get(ref, s, t) });
    assert.ok(residues > 20000, `checked ${residues} residues`);
    assert.deepEqual(failures.slice(0, 10), []);
  });
});

/** /translation of every CDS in the nucleotide fixtures. */
function genbankProteins(): Map<string, string> {
  const out = new Map<string, string>();
  for (const name of NUCLEOTIDE_GB) {
    const text = fixture(name);
    for (const m of text.matchAll(/\/protein_id="([^"]+)"[\s\S]*?\/translation="([^"]+)"/g)) {
      out.set(`${/^[A-Z]{2}_/.test(m[1]!) ? "refseq" : "insdc"}:${m[1]}`, m[2]!.replace(/\s+/g, ""));
    }
  }
  return out;
}

describe("assembly report as sequences and names (GRCh38.p14 chr13 rows, real data)", () => {
  const text = fixture("GRCh38.p14_assembly_report_chr13.txt");
  it("lists the sequences with length and species", () => {
    const [chr13] = assemblyReportSequences(text, "report.txt").filter((r) => r.ref === "refseq:NC_000013.11");
    assert.deepEqual([chr13!.length, chr13!.moltype, chr13!.taxon, chr13!.provenance.adapter], [114364328, "DNA", 9606, "assembly-report"]);
  });
  it("maps the assembly's own names to RefSeq accessions", () => {
    const aliases = assemblyReportAliases(text);
    assert.deepEqual([aliases["13"], aliases["chr13"], aliases["CM000675.2"]], ["refseq:NC_000013.11", "refseq:NC_000013.11", "refseq:NC_000013.11"]);
  });
  it("uses GenBank accessions for an INSDC-only assembly (MpTak_v7.1, real data)", () => {
    const report = fixture("MpTak_v7.1_assembly_report.txt");
    assert.deepEqual([assemblyReportAliases(report)["1"], assemblyReportAliases(report)["MT"]], ["insdc:AP031342.1", "insdc:AP025456.1"]);
    const seqs = assemblyReportSequences(report);
    assert.deepEqual([seqs.length, seqs.find((r) => r.ref === "insdc:AP025455.1")?.topology, seqs[0]!.taxon], [12, "circular", 1480154]);
    assert.equal(assemblyReportInfo(report).released, "2024-03-27");
  });
});
