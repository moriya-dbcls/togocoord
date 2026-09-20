// Namespaces of Location IDs (spec-core §3.4).
import { LocationSemanticError } from "./errors.ts";

export type Unit = "nt" | "aa";

export interface NamespaceDef {
  /** identifiers.org / bioregistry style prefix, lower case. */
  prefix: string;
  /** Whole-accession pattern (including version / isoform / chain suffixes). */
  pattern: RegExp;
  normalize?: (accession: string) => string;
  defaultUnit?: (accession: string) => Unit | undefined;
}

const PDB_RE = /^(?:([0-9][A-Za-z0-9]{3})|(pdb_[0-9]{4}[0-9][A-Za-z0-9]{3}))\.([A-Za-z0-9]{1,4})$/i;

export const DEFAULT_NAMESPACES: readonly NamespaceDef[] = [
  {
    prefix: "insdc",
    pattern: /^[A-Z]{1,6}\d{5,}(?:\.\d+)?$/,
    defaultUnit: (acc) => (/^[A-Z]{3}\d/.test(acc) ? "aa" : "nt"),
  },
  {
    prefix: "refseq",
    pattern: /^[A-Z]{2}_[A-Z0-9]+(?:\.\d+)?$/,
    defaultUnit: (acc) => (/^(?:NP|XP|YP|WP|AP|ZP)_/.test(acc) ? "aa" : "nt"),
  },
  {
    prefix: "ensembl",
    pattern: /^ENS[A-Z]*[EGTPR]\d{11}(?:\.\d+)?$/,
    defaultUnit: (acc) => (/^ENS[A-Z]*P\d{11}/.test(acc) ? "aa" : "nt"),
  },
  {
    prefix: "uniprot",
    pattern: /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$/,
    defaultUnit: () => "aa",
  },
  {
    prefix: "uniparc",
    pattern: /^UPI[0-9A-F]{10}$/,
    defaultUnit: () => "aa",
  },
  {
    prefix: "pdb",
    pattern: PDB_RE,
    normalize: normalizePdb,
    defaultUnit: () => "aa",
  },
  {
    prefix: "refget",
    pattern: /^SQ\.[A-Za-z0-9_-]{32}$/,
  },
];

/** `4hhb.A` -> `4HHB.A`, `pdb_00004hhb.A` -> `4HHB.A`, other extended IDs lower-cased. Chains keep their case. */
function normalizePdb(accession: string): string {
  const m = PDB_RE.exec(accession)!;
  const chain = m[3]!;
  if (m[1] !== undefined) return `${m[1].toUpperCase()}.${chain}`;
  const ext = m[2]!.toLowerCase();
  return ext.startsWith("pdb_0000") ? `${ext.slice(8).toUpperCase()}.${chain}` : `${ext}.${chain}`;
}

export class NamespaceRegistry {
  readonly #defs = new Map<string, NamespaceDef>();

  constructor(defs: Iterable<NamespaceDef> = DEFAULT_NAMESPACES) {
    for (const def of defs) this.register(def);
  }

  register(def: NamespaceDef): this {
    this.#defs.set(def.prefix.toLowerCase(), def);
    return this;
  }

  get(prefix: string): NamespaceDef | undefined {
    return this.#defs.get(prefix.toLowerCase());
  }

  /** Registered prefixes, in registration order. */
  prefixes(): string[] {
    return [...this.#defs.keys()];
  }

  /** Validate and normalise `namespace:accession` into the internal sequence key. */
  refKey(namespace: string, accession: string): string {
    const def = this.get(namespace);
    if (!def) throw new LocationSemanticError(`unknown namespace '${namespace}'`);
    if (!def.pattern.test(accession)) {
      throw new LocationSemanticError(`'${accession}' is not a valid ${def.prefix} accession`);
    }
    return `${def.prefix}:${def.normalize ? def.normalize(accession) : accession}`;
  }

  defaultUnit(ref: string): Unit | undefined {
    const { namespace, accession } = splitRef(ref);
    return this.get(namespace)?.defaultUnit?.(accession);
  }
}

/** Split an internal key `namespace:accession`. */
export function splitRef(ref: string): { namespace: string; accession: string } {
  const i = ref.indexOf(":");
  if (i < 0) throw new LocationSemanticError(`'${ref}' is not a namespace:accession key`);
  return { namespace: ref.slice(0, i), accession: ref.slice(i + 1) };
}
