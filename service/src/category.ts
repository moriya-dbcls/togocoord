// Coarse sequence categories used as conversion targets ("to genome", "to protein", ...).
import { splitRef, type Unit } from "@togocoord/core";

export type Category = "genome" | "transcript" | "protein" | "gene_region" | "other";

export const CATEGORIES: readonly Category[] = ["genome", "transcript", "protein", "gene_region", "other"];

/**
 * Category of a sequence key from its unit, molecule type and accession conventions:
 * RefSeq NC_/NT_/NW_/NZ_/AC_ genome, NM_/NR_/XM_/XR_ transcript, NG_ gene region; Ensembl ENST/ENSP;
 * INSDC by molecule type.
 */
export function categoryOf(ref: string, info: { unit?: Unit; moltype?: string } = {}): Category {
  if (info.unit === "aa" || info.moltype === "protein") return "protein";
  const { namespace, accession } = splitRef(ref);
  switch (namespace) {
    case "refseq":
      if (/^(?:NM|NR|XM|XR)_/.test(accession)) return "transcript";
      if (/^NG_/.test(accession)) return "gene_region";
      if (/^(?:NC|NT|NW|NZ|AC)_/.test(accession)) return "genome";
      if (/^(?:NP|XP|YP|WP|AP)_/.test(accession)) return "protein";
      return "other";
    case "ensembl":
      if (/^ENS[A-Z]*T\d/.test(accession)) return "transcript";
      if (/^ENS[A-Z]*P\d/.test(accession)) return "protein";
      return "other";
    case "uniprot":
    case "uniparc":
    case "pdb":
      return "protein";
    default:
      if (info.moltype && /RNA/i.test(info.moltype)) return "transcript";
      if (info.moltype === "DNA") return "genome";
      return "other";
  }
}
