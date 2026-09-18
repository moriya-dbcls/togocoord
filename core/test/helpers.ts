import { createContext, NamespaceRegistry, type ContextOptions, type CoordContext } from "../src/index.ts";

/** Registry with the default namespaces plus a permissive `test` namespace (nt unless overridden). */
export function testRegistry(): NamespaceRegistry {
  return new NamespaceRegistry().register({
    prefix: "test",
    pattern: /^[A-Za-z0-9_-]+$/,
    defaultUnit: () => "nt",
  });
}

export function testContext(options: Omit<ContextOptions, "registry"> = {}): CoordContext {
  return createContext({ ...options, registry: testRegistry() });
}
