/** Parse FASTA into id (first word of the header) -> residues (upper case). */
export function parseFasta(text: string): Map<string, string> {
  const out = new Map<string, string>();
  let id: string | undefined;
  let chunks: string[] = [];
  const flush = () => {
    if (id !== undefined) out.set(id, chunks.join("").toUpperCase());
  };
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      flush();
      id = line.slice(1).trim().split(/\s+/)[0] ?? "";
      chunks = [];
    } else if (id !== undefined) {
      chunks.push(line.replace(/\s+/g, ""));
    }
  }
  flush();
  return out;
}
