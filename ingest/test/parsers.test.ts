import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { parseGenBank, qualifier, qualifiers } from "../src/gbff.ts";
import { parseGff3 } from "../src/gff3.ts";
import { parseFasta } from "../src/fasta.ts";
import { refgetDigest, reverseComplement, translate } from "../src/sequence.ts";
import { fixture } from "./helpers.ts";

const SNIPPET = `LOCUS       AB000001                  30 bp    DNA     circular BCT 01-JAN-2000
DEFINITION  Test record
            spanning two lines.
ACCESSION   AB000001
VERSION     AB000001.1
FEATURES             Location/Qualifiers
     source          1..30
                     /organism="Test organism"
                     /db_xref="taxon:12345"
     CDS             join(1..6,
                     10..18)
                     /codon_start=1
                     /note="a note that is
                     continued"
                     /pseudo
                     /translation="MKK
                     GG"
ORIGIN
        1 atgaaaccca agggcggctg attttttttt
//
`;

describe("GenBank flat file parser", () => {
  const [r] = parseGenBank(SNIPPET);
  it("reads the header", () => {
    assert.equal(r?.version, "AB000001.1");
    assert.equal(r?.length, 30);
    assert.equal(r?.topology, "circular");
    assert.equal(r?.moltype, "DNA");
    assert.equal(r?.definition, "Test record spanning two lines.");
    assert.equal(r?.sequence, "ATGAAACCCAAGGGCGGCTGATTTTTTTTT");
  });
  it("joins location and qualifier continuations", () => {
    const cds = r!.features[1]!;
    assert.equal(cds.location, "join(1..6,10..18)");
    assert.equal(qualifier(cds, "note"), "a note that is continued");
    assert.equal(qualifier(cds, "translation"), "MKKGG");
    assert.equal(qualifier(cds, "pseudo"), "");
    assert.deepEqual(qualifiers(r!.features[0]!, "db_xref"), ["taxon:12345"]);
  });
  it("parses every fixture record with its full sequence", () => {
    for (const name of ["NC_045512.2.gb", "NC_012920.1.gb", "NC_001405.1.gb", "NC_002127.1.gb", "NM_000581.4.gb", "NM_014739.3.gb", "NP_055554.1.gp"]) {
      for (const rec of parseGenBank(fixture(name))) assert.equal(rec.sequence?.length, rec.length, name);
    }
  });
});

describe("GFF3 parser", () => {
  const doc = parseGff3(
    "##gff-version 3\n##sequence-region chr 1 100\n" +
      "chr\tsrc\tCDS\t1\t10\t.\t+\t0\tID=c1;Note=a%3Bb;Dbxref=A:1,B:2\n" +
      "chr\tsrc\tCDS\t20\t31\t.\t+\t2\tID=c1\n" +
      "chr\tsrc\texon\t1\t10\t.\t+\t.\tParent=t1\n" +
      "##FASTA\n>chr desc\nacgt\nAC\n",
  );
  it("groups rows by ID and decodes attributes", () => {
    assert.equal(doc.features.length, 2);
    assert.equal(doc.features[0]!.rows.length, 2);
    assert.deepEqual(doc.features[0]!.attributes.Note, ["a;b"]);
    assert.deepEqual(doc.features[0]!.attributes.Dbxref, ["A:1", "B:2"]);
    assert.equal(doc.features[0]!.rows[1]!.phase, 2);
  });
  it("reads regions and the FASTA section", () => {
    assert.equal(doc.regions.get("chr"), 100);
    assert.equal(doc.sequences.get("chr"), "ACGTAC");
  });
});

describe("sequence helpers", () => {
  it("refget digest matches the GA4GH example", () => {
    assert.equal(refgetDigest("ACGT"), "SQ.aKF498dAxcJAqme6QYQ7EZ07-fiw8Kw2");
    assert.equal(refgetDigest("acgt"), refgetDigest("ACGT"));
  });
  it("translates with NCBI genetic codes", () => {
    assert.equal(translate("ATGTGAAGAATA", 1), "M*RI");
    assert.equal(translate("ATGTGAAGAATA", 2), "MW*M");
    assert.equal(translate("ATGNNN", 1), "MX");
  });
  it("reverse-complements IUPAC codes", () => {
    assert.equal(reverseComplement("ACGTRYN"), "NRYACGT");
  });
  it("parses FASTA", () => {
    assert.deepEqual([...parseFasta(">a x\nAC\ngt\n>b\nNN\n")], [["a", "ACGT"], ["b", "NN"]]);
  });
});
