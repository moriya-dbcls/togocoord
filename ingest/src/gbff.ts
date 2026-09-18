// GenBank / GenPept flat file parser (INSDC feature table layout).

export interface GbFeature {
  key: string;
  /** Location text with whitespace removed. */
  location: string;
  qualifiers: Array<[name: string, value: string | true]>;
}

export interface GbRecord {
  name: string;
  length: number;
  unit: "bp" | "aa";
  moltype?: string;
  topology?: "linear" | "circular";
  division?: string;
  accession?: string;
  /** accession.version from the VERSION line. */
  version?: string;
  definition?: string;
  features: GbFeature[];
  /** Residues from ORIGIN (lower case removed), if present. */
  sequence?: string;
}

/** Parse every record of a (multi-record) GenBank or GenPept flat file. */
export function parseGenBank(text: string): GbRecord[] {
  const records: GbRecord[] = [];
  let lines: string[] = [];
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith("//")) {
      if (lines.some((l) => l.startsWith("LOCUS"))) records.push(parseGenBankRecord(lines));
      lines = [];
    } else {
      lines.push(line);
    }
  }
  if (lines.some((l) => l.startsWith("LOCUS"))) records.push(parseGenBankRecord(lines));
  return records;
}

/** Parse the lines of one record (without the terminating `//`). */
export function parseGenBankRecord(lines: string[]): GbRecord {
  const locus = lines.find((l) => l.startsWith("LOCUS"))!;
  const tokens = locus.split(/\s+/).slice(1);
  const unit = tokens[2] === "aa" ? "aa" : "bp";
  const record: GbRecord = { name: tokens[0]!, length: Number(tokens[1]), unit, features: [] };
  const rest = tokens.slice(3);
  const moltype = rest.find((t) => /(?:DNA|RNA)$/i.test(t));
  if (moltype) record.moltype = moltype;
  const topology = rest.find((t) => t === "linear" || t === "circular");
  if (topology === "linear" || topology === "circular") record.topology = topology;
  const division = rest.find((t) => /^[A-Z]{3}$/.test(t));
  if (division) record.division = division;

  let section = "";
  const definition: string[] = [];
  const featureLines: string[] = [];
  const sequence: string[] = [];
  for (const line of lines) {
    const head = /^[A-Z]/.test(line) ? line.split(/\s+/)[0]! : "";
    if (head) section = head;
    switch (section) {
      case "ACCESSION":
        if (head) record.accession = line.split(/\s+/)[1];
        break;
      case "VERSION":
        if (head) record.version = line.split(/\s+/)[1];
        break;
      case "DEFINITION":
        definition.push(line.slice(12).trim());
        break;
      case "FEATURES":
        if (!head) featureLines.push(line);
        break;
      case "ORIGIN":
        if (!head) sequence.push(line.replace(/[^A-Za-z]/g, ""));
        break;
    }
  }
  if (definition.length) record.definition = definition.join(" ");
  record.features = parseFeatures(featureLines);
  if (sequence.length) record.sequence = sequence.join("").toUpperCase();
  return record;
}

function parseFeatures(lines: string[]): GbFeature[] {
  const features: GbFeature[] = [];
  let current: { key: string; location: string; qualifiers: Array<[string, string | true]>; open?: [string, string] } | undefined;

  const finish = () => {
    if (!current) return;
    if (current.open) current.qualifiers.push([current.open[0], unquote(current.open[1])]);
    features.push({ key: current.key, location: current.location, qualifiers: current.qualifiers });
  };

  for (const line of lines) {
    if (/^ {5}\S/.test(line)) {
      finish();
      current = { key: line.slice(5, 21).trim(), location: line.slice(21).replace(/\s+/g, ""), qualifiers: [] };
      continue;
    }
    if (!current) continue;
    const body = line.slice(21);
    if (current.open) {
      // Continuation of a quoted value: translations are joined without spaces.
      const [name, value] = current.open;
      current.open = [name, value + (name === "translation" ? "" : " ") + body.trim()];
      if (quotesClosed(current.open[1])) {
        current.qualifiers.push([name, unquote(current.open[1])]);
        delete current.open;
      }
    } else if (body.startsWith("/")) {
      const m = /^\/([^=]+)(?:=(.*))?$/.exec(body.trim())!;
      const name = m[1]!;
      if (m[2] === undefined) {
        current.qualifiers.push([name, true]);
      } else if (m[2].startsWith('"') && !quotesClosed(m[2])) {
        current.open = [name, m[2]];
      } else {
        current.qualifiers.push([name, unquote(m[2])]);
      }
    } else if (current.qualifiers.length === 0) {
      current.location += body.replace(/\s+/g, "");
    }
  }
  finish();
  return features;
}

function quotesClosed(value: string): boolean {
  return (value.match(/"/g)?.length ?? 0) % 2 === 0;
}

function unquote(value: string): string {
  return value.startsWith('"') && value.endsWith('"') ? value.slice(1, -1).replace(/""/g, '"') : value;
}

export function qualifier(f: GbFeature, name: string): string | undefined {
  const q = f.qualifiers.find(([n]) => n === name);
  return q === undefined ? undefined : q[1] === true ? "" : q[1];
}

export function qualifiers(f: GbFeature, name: string): string[] {
  return f.qualifiers.filter(([n]) => n === name).map(([, v]) => (v === true ? "" : v));
}
