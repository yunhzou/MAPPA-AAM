# Frozen research source

`campaigns/<report-name>/` contains original drivers previously stored alongside
`reports/<report-name>/`. `manuscript_visuals/` contains replaced movie builders
and early evidence-preparation scripts. Their bytes are unchanged, so recorded
source hashes still identify the original code.

These scripts expect their recorded frozen-engine and staging layouts; this
archive is not a portable current-API entry point. Use the [benchmark guide](../README.md)
for new work. Do not run old evidence-preparation scripts into the current paper.
For exact historical reproduction, use a separate checkout of the manifest's
base revision and the original campaign's input/environment records:

```sh
git worktree add --detach ../mappa-historical 59f55e4
```

The [relocation manifest](../../docs/publication_cleanup_20260916.json) maps every
original source path to its current location and checksum. Numerical evidence
and witnesses stay in `reports/`. Old visual artifacts can be recovered from
the same Git revision; no history was rewritten.
