// Guessing the namespace of an input written without one (spec-service §6): `NP_000572.2:49` is read as
// `refseq:NP_000572.2:49`. Candidates come from TogoID's curated ID syntax (togoid-patterns.ts), which is stricter
// than the registry's own patterns (ours accept `P07203` as an INSDC accession too), plus the registry's patterns for
// the namespaces TogoID does not list. Some ID forms genuinely belong to two databases (`P07203` is both a UniProt
// accession and an INSDC one), so a candidate that we actually hold wins, and a tie is reported rather than guessed.
import type { NamespaceRegistry } from "@togocoord/core";
import { TOGOID_PATTERNS } from "./togoid-patterns.ts";

export interface Guess {
  /** `namespace:accession`, normalised by the registry. */
  ref: string;
  namespace: string;
  /** The sequence is loaded in a store: the strongest evidence, whatever the patterns say. */
  known: boolean;
}

/** Databases an accession could belong to, the ones we hold first. */
export function guessNamespaces(accession: string, registry: NamespaceRegistry, known: (ref: string) => boolean): Guess[] {
  const namespaces = new Set<string>();
  const byTogoid = new Set<string>();
  for (const { namespace, pattern } of TOGOID_PATTERNS) {
    byTogoid.add(namespace);
    if (pattern.test(accession)) namespaces.add(namespace);
  }
  // TogoID identifies a PDB entry (4HHB), we address a chain of it (4HHB.A); everything TogoID does not list at all
  // (UniParc, refget) keeps its registry pattern.
  for (const prefix of registry.prefixes()) {
    const def = registry.get(prefix)!;
    const isPdbChain = prefix === "pdb" && def.pattern.test(accession);
    if (isPdbChain || (!byTogoid.has(prefix) && def.pattern.test(accession))) namespaces.add(prefix);
  }

  const guesses: Guess[] = [];
  for (const namespace of namespaces) {
    let ref;
    try {
      ref = registry.refKey(namespace, accession); // rejects what we cannot address, e.g. a PDB entry without its chain
    } catch {
      continue;
    }
    guesses.push({ ref, namespace, known: known(ref) });
  }
  return guesses.sort((a, b) => Number(b.known) - Number(a.known) || a.namespace.localeCompare(b.namespace));
}

/**
 * The database of an accession: undefined when nothing matches (the caller then reports the input as written), and
 * `ambiguous` when several databases match equally well, for an error that asks the user to write the database.
 */
export function inferNamespace(
  accession: string,
  registry: NamespaceRegistry,
  known: (ref: string) => boolean,
): { ref?: string; namespace?: string; ambiguous?: string[] } {
  const guesses = guessNamespaces(accession, registry, known);
  if (guesses.length === 0) return {};
  const best = guesses.filter((g) => g.known === guesses[0]!.known);
  if (best.length === 1) return { ref: best[0]!.ref, namespace: best[0]!.namespace };
  return { ambiguous: best.map((g) => g.ref) };
}
