export { LocationSemanticError, LocationSyntaxError, MappingError } from "./errors.ts";
export { parseLocationText, type AstNode, type AstPos, type Codon } from "./parser.ts";
export { DEFAULT_NAMESPACES, NamespaceRegistry, splitRef, type NamespaceDef, type Unit } from "./namespace.ts";
export {
  canonicalize,
  createContext,
  formatLocation,
  formatLocationId,
  parseLocationId,
  splitLocationId,
  type CodonMode,
  type ContextOptions,
  type CoordContext,
  type Location,
  type Segment,
} from "./location.ts";
export {
  cdsMapping,
  compose,
  invert,
  Mapping,
  mappingFromLocation,
  projectInterval,
  residueBlock,
  type Block,
  type BlockHit,
  type CdsOptions,
  type ResidueBlockSpec,
} from "./mapping.ts";
export { mapLocation, type MapResult, type Piece, type TargetLocation } from "./convert.ts";
export {
  decodeLocationId,
  DEFAULT_BASE,
  DEFAULT_VOCABULARY,
  encodeLocationId,
  FALDO,
  locationIri,
  RDF,
  toFaldo,
  type FaldoOptions,
} from "./faldo.ts";
