# Where things stand

[English](status.md) | [日本語](status.ja.md)

A resume log: what is running, what is waiting for a decision, and what to pick up next. The specs hold the reasoning;
this file holds only what is not yet settled. Last written 2026-10-19 at `6d90ec6`.

## State

- `main` at `6d90ec6`, working tree clean, everything pushed to `github.com/moriya-dbcls/togocoord`.
- Tests pass: core 119, ingest 71, service 62. Node 24 (`npm test`, `npm run typecheck`).
- 51 stores built in `data/stores` (10 GB), listed in `scripts/stores.txt`. Rebuild with `node scripts/data.ts build`
  (docs/data.md); data is never committed.

## Deployment is behind the repository

**The 504 fix is not deployed.** A deployment behind nginx was returning 504s, diagnosed as head-of-line blocking
(spec-service §15) and fixed in `da47887` / `6d90ec6`, but that needs a pull and a restart to take effect:

```sh
git pull
node --no-warnings service/src/serve.ts --host 127.0.0.1 --port 8080 --base https://<host>/<path>/ \
  --workers auto --log --stores scripts/stores.txt --store-dir /path/to/stores
```

`--workers` is what stops one slow conversion from blocking the service; `--log` is how the next such incident gets
diagnosed. Do not raise nginx's `proxy_read_timeout` — a request over a minute is a bug on this side.

Also not deployed, in the same pull: the `/v1/meta` fix (`c8c0a96`, 10.2 s → under 10 ms, which on its own could
cause 504s on every page load), namespace inference, the JSON link, the favicon, `coverage` / `minMatch`.

Four stores built since the last data copy and needed on the server (95 MB): `marchantia_mpv4`,
`marchantia_uniprot_v4`, `mp_v4_to_hub`, `mp_hub_to_v4`.

**Unverified**: which URL the deployment answers on. `https://sparql-support.dbcls.jp/togocoord/v1/meta` returned 301
and `/togocoord_v1/` returned 302 to a 404 from here on 2026-10-08, so neither could be measured. Confirm the URL
before testing a deployment.

## Waiting on a decision

- **Orthologs.** Measured: OMA has 9,668 Arabidopsis↔Marchantia pairs with `sequence_md5`, which matches our own
  sequence keys, and sampled pairs sit around 50% identity. The user's position (2026-09-25) is that a coordinate
  service earns its keep by letting one follow variants, domains and secondary structure — not by restating what
  BLAST or OMA already answer — so aligning ortholog pairs crudely is the wrong move. An ortholog edge should rest on
  a shared coordinate system (an InterPro/Pfam HMM's match states, or a structural superposition) whose basis stays
  visible in the result. Not started; design.md §12 still lists the source as undecided.
- **Residue-level annotation, the stronger candidate.** UniProt is ingested as sequences only, so active sites,
  modified residues, domains, secondary structure and variants are all absent. The suggested order was UniProt
  features → structure-derived (DSSP, AlphaFold pLDDT) → variants (TogoVar/ClinVar), each as annotations on existing
  sequences, plus widening `/v1/annotations` to collect annotations of reachable layers. Not started.
- **T2 for species without `idmapping_selected.tab`** (Marchantia): UniRef90 clusters from the UniProt REST API would
  give candidates for 6,047 of the UP000077202 entries that no identical sequence reaches. Not started; Mp_v4 shrank
  this to roughly 400 entries (cDNA-derived, other species in the genus, and 45 mitochondrial).

## Decided, so do not re-open

- A position with several chain mappings converts to one place; no `-multiple` (spec-service §2).
- Stores never hold residues or sequences — only difference positions (spec-ingest §19).
- `MpTak1_v7.1` stays out: same accessions as the hub, for sequences and CDS proteins alike (spec-ingest §20).
- `minMatch` has no default: partial results are ordinary here, unlike UCSC liftOver (spec-service §4).
- Ordering by chain score was considered for the BRCA1 CRE incident and does not help — the wrong locus scores higher.

## Smaller things left open

- `togocoord-serve` has no rate limiting, authentication or concurrency cap; the reverse proxy should provide them.
- Col-CEN v1.2 (Arabidopsis) has no INSDC accession, so its sequences have no key (spec-ingest §20).
- design.md §12: `faldo:codonPosition` to the FALDO developers, the domain name and identifiers.org prefix, where the
  service runs, and whether to log queries that find no path.
