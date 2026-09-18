/** Malformed Location ID text. `position` is the character offset (whitespace removed). */
export class LocationSyntaxError extends Error {
  readonly position: number;
  constructor(message: string, position: number) {
    super(`${message} (at offset ${position})`);
    this.name = "LocationSyntaxError";
    this.position = position;
  }
}

/** Well-formed text that is meaningless for the referenced sequences (spec-core §3.3). */
export class LocationSemanticError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "LocationSemanticError";
  }
}

/** Invalid blocks or mapping construction parameters (spec-core §4). */
export class MappingError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "MappingError";
  }
}
