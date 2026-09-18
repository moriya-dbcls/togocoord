// Parser for the location part of a Location ID (spec-core §3.2).
import { LocationSyntaxError } from "./errors.ts";

export type Codon = 1 | 2 | 3;

export interface AstPos {
  value: number;
  codon?: Codon;
}

export type AstNode =
  | { type: "point"; ref?: string; pos: AstPos }
  | { type: "range"; ref?: string; begin: AstPos; end: AstPos; fuzzyLow: boolean; fuzzyHigh: boolean }
  | { type: "between"; ref?: string; left: AstPos; right: AstPos }
  | { type: "oneof"; ref?: string; begin: number; end: number }
  | { type: "complement"; ref?: string; child: AstNode }
  | { type: "join" | "order"; ref?: string; children: AstNode[] };

const KEYWORDS = ["complement", "join", "order"] as const;
const REMOTE_RE = /^([A-Za-z0-9_.-]+):/;
const POS_RE = /^(\d+)(?:c(\d))?/;

export function parseLocationText(text: string): AstNode {
  const parser = new Parser(text.replace(/\s+/g, ""));
  const node = parser.location();
  if (!parser.atEnd()) parser.fail("unexpected trailing characters");
  return node;
}

class Parser {
  private readonly s: string;
  private i = 0;

  constructor(s: string) {
    this.s = s;
  }

  atEnd(): boolean {
    return this.i >= this.s.length;
  }

  fail(message: string, at: number = this.i): never {
    throw new LocationSyntaxError(message, at);
  }

  private peek(): string | undefined {
    return this.s[this.i];
  }

  private expect(ch: string): void {
    if (this.peek() !== ch) this.fail(`expected '${ch}'`);
    this.i++;
  }

  location(): AstNode {
    const ref = this.remote();
    const node = this.body();
    if (ref !== undefined) {
      if (node.ref !== undefined) this.fail("nested remote references");
      node.ref = ref;
    }
    return node;
  }

  private remote(): string | undefined {
    const m = REMOTE_RE.exec(this.s.slice(this.i));
    if (!m || !/[A-Za-z]/.test(m[1]!)) return undefined;
    this.i += m[0].length;
    return m[1];
  }

  private body(): AstNode {
    for (const kw of KEYWORDS) {
      if (!this.s.startsWith(`${kw}(`, this.i)) continue;
      this.i += kw.length + 1;
      if (kw === "complement") {
        const child = this.location();
        this.expect(")");
        return { type: "complement", child };
      }
      const children = [this.location()];
      while (this.peek() === ",") {
        this.i++;
        children.push(this.location());
      }
      this.expect(")");
      return { type: kw, children };
    }
    return this.span();
  }

  private span(): AstNode {
    const start = this.i;
    let fuzzyLow = false;
    if (this.peek() === "<") {
      fuzzyLow = true;
      this.i++;
    }
    const a = this.pos();
    if (this.s.startsWith("..", this.i)) {
      this.i += 2;
      let fuzzyHigh = false;
      if (this.peek() === ">") {
        fuzzyHigh = true;
        this.i++;
      } else if (this.peek() === "<") {
        this.fail("'<' is only allowed before the begin position");
      }
      const b = this.pos();
      if (orderKey(a, 1) > orderKey(b, 3)) this.fail("range begin is greater than end", start);
      return { type: "range", begin: a, end: b, fuzzyLow, fuzzyHigh };
    }
    if (fuzzyLow) this.fail("fuzzy marker on a single position", start);
    if (this.peek() === "^") {
      this.i++;
      return { type: "between", left: a, right: this.pos() };
    }
    if (this.peek() === ".") {
      this.i++;
      const b = this.pos();
      if (a.codon !== undefined || b.codon !== undefined) this.fail("codon extension in a one-of span", start);
      if (a.value > b.value) this.fail("one-of begin is greater than end", start);
      return { type: "oneof", begin: a.value, end: b.value };
    }
    return { type: "point", pos: a };
  }

  private pos(): AstPos {
    if (this.peek() === ">") this.fail("'>' is only allowed before the end position");
    const m = POS_RE.exec(this.s.slice(this.i));
    if (!m) this.fail("expected a position");
    const value = Number(m[1]);
    if (!Number.isSafeInteger(value) || value < 1) this.fail("positions must be integers >= 1");
    const pos: AstPos = { value };
    if (m[2] !== undefined) {
      const c = Number(m[2]);
      if (c < 1 || c > 3) this.fail("codon position must be 1, 2 or 3");
      pos.codon = c as Codon;
    }
    this.i += m[0].length;
    return pos;
  }
}

/** Sortable key of a position; `defaultCodon` fills in a missing codon extension. */
function orderKey(p: AstPos, defaultCodon: Codon): number {
  return p.value * 3 + (p.codon ?? defaultCodon);
}
