# Publication cleanup — 16 September 2026

The [manifest](publication_cleanup_20260916.json) records the original/new paths,
checksums, sizes, removal reasons, and preserved Git revision `59f55e4954cb7fb106ec66fe0912f1903550883c`.

- 85 source/contract files relocated: historical report drivers to
  `bench/archive/campaigns/`, replaced visualization builders to
  `bench/archive/manuscript_visuals/`, benchmark utilities to `bench/`,
  and regression contracts to `bench/contracts/`.
- 46 superseded generated media or unused template artifacts
  removed from the current checkout (107.0 MiB).
- Current manuscript figures, verified Golden film, notebooks, library code,
  benchmark measurements, witnesses and original provenance retained.
- Relocated source is byte-identical; no historical measurement or source-hash
  record was rewritten. Current documentation links were repaired.

Every affected file matched its committed version before cleanup and was
independently backed up and checksum-verified. A scan of 204 current library,
test, manuscript and build files found no obsolete-artifact references.
Git history is unchanged. The cleanup reduces current checkout size, not the
historical Git object database. Historical figures/viewers remain accessible
at the manifest's base revision.

Validation: the README's actual Python example ran on the notebook's embedded
molecules; all four relocated profiling CLIs and the final pipeline CLI passed
`--help`; 38 focused tests passed. No newly broken Markdown links were found.
Offline checks passed for the retained comparison, catalogue, trajectory, and
Golden research film, with no JavaScript errors or network requests.
The 20-page paper rebuilt with its numerical/source/reference checks passing;
the final title page was visually reviewed. See the
[validation record](publication_cleanup_validation_20260916.json).

## Repository-wide cleanup — 28 October 2026

Removed 68 retired campaign/visualization files (67 Python modules) from the
current checkout, plus an unused MongoDB Compose file and its optional pymongo
dependency. The exact previous Git revision and checksums are recorded in the
[removal manifest](validation/repository_cleanup_20261028.json). No maintained
library/test/tool/launcher imports the removed archive paths. The archive index
now points to Git history. Existing measurements, witnesses, data, manuscript
figures and uncommitted Desktop work are preserved. The active benchmark tools
and overlapping design documents still need consolidation; this removal does
not claim that the entire repository is clean.
