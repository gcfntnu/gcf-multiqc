# MultiQC migration and acceptance — 2 November 2026

## Decision

Use unmodified upstream MultiQC plus the independently versioned `gcf-multiqc`
extension. The candidate preserves the report commands and input conventions
used by BFQ and the workflows. The old fork remains available for comparison
and rollback; it has not been modified.

The extension belongs in the reporting container. BFQ owns orchestration and
its installed `multiqc` Apptainer wrapper. `gcf-workflows/docker.config` owns the
reporting image default, while explicit project configuration retains precedence.
Do not install the reporting engine into BFQ's command environment.

## Current artifact

- Engine: MultiQC **1.35**; extension: **0.1.0rc4**; Plotly: **6.0.1**.
- Image: `gcfntnu/gcf-multiqc:1.35-gcf0.1.0-rc4`, Linux amd64.
- Registry digest: `sha256:2c356684c03d24c11e63be5724d4f107c253eccd8be98c236952bd3df5435878`.
- Installed plugin source: `4625e75315f5c89475af5ce70313642f4154ce40`.

The image recipe pins its base digest, plugin source archive/checksum and Python
wheel dependencies. Plotly remains below 6.1 because MultiQC 1.35 uses Kaleido
0.2.1 for flat plots. Debian package repositories are not frozen, so rebuilding
creates a separately identified artifact. The workflow selects the accepted
registry digest rather than the moving `latest` tag. Earlier rc1/rc2/rc3 tags retain
their original contents.

## Reviewed source baseline

| Source | Revision |
| --- | --- |
| Legacy `gcfntnu/MultiQC`, `bfq-dev` | `d3fce9dab5b8522ad022f04161a145db121ad23d` |
| Upstream MultiQC v1.35 | `87e504bb77e3d94687dd980638724b7ef8e1e92f` |
| Workflows adoption base, `bfq-dev` | `26fafe406d856b49beead1ecce251cc943e56e63` |
| BFQ adoption base, `bfq-dev` | `2b1fbc51f625c128e9aa5a083f1c5fe1abb06068` |
| Docker recipe base, `master` | `28e66251ec003b0f0d65526d2fa6b8cc03ed92b3` |

The seven copied public workflow templates remain unchanged. Synthetic fixtures
use invented sample values, and their READMEs identify source conventions.
Private operational inputs, reports and detailed comparison evidence remain
outside Git and published images.

## Implemented compatibility

The package retains the selectors used by RNA-seq, microbiome, single-cell,
small-RNA, default, WGS and metagenome reports. GCF-specific parsers are small
standalone modules. Standard parsing and plotting reuse upstream implementations
where possible. Supported hooks select compatibility adapters without competing
entry points; `gcf_compatibility: false` retains an explicit native opt-out.

The custom-content adapter preserves literal sample identities and metadata,
including leading zeros and missing-value-like strings. Existing PCA/UMAP data
are adapted in memory, preserving coordinates, sample tooltips, groups, colors
and component axes. Producer inputs are not rewritten or recalculated.
The separate `bracken` selector keeps estimates distinct from Kraken observations.
The unused old template and `fastv` parser were not copied.

[The README](../README.md#deliberate-corrections) defines intentional changes:
Unitas fractions use each sample's own total; STAR aligned-read counts exclude
reads mapped to too many loci; splice counts describe events rather than a
partition of reads; undefined ratios remain unavailable. Numeric exports and
malformed/duplicate-input handling also have explicit contracts.

Functional acceptance and review found four additional defects, now repaired:

- [Issue #3](https://github.com/gcfntnu/gcf-multiqc/issues/3): plain-text General
  Statistics metadata was coerced during rendering/export. Literal text is now
  preserved and HTML-escaped without changing scientific numeric columns.
- [Issue #4](https://github.com/gcfntnu/gcf-multiqc/issues/4): repeated BCL Convert
  index rows overwrote sample/lane chunks. Counts now accumulate, quality metrics
  retain their weights and all distinct indices remain available. Aggregate
  sample totals and the aggregate sample-detail table were already correct.
- [Issue #5](https://github.com/gcfntnu/gcf-multiqc/issues/5): combined runs that
  reuse lane numbers need consistent labels across all samples. Collision
  detection now uses all runs, preserving short labels for single-run reports.

- [Issue #6](https://github.com/gcfntnu/gcf-multiqc/issues/6): global filename
  cleaning merged paired FastQC mates and retained whichever archive was read
  last. FastQC now retains separate R1/R2 records and plots while other modules
  keep their sample-level metadata joins.

The companion BFQ change also removes duplicate `custom_content` selection for
separate conversion passes. Native demultiplexer selection, report paths and
zero/unavailable assignment semantics remain unchanged.

## Completed validation

The current source passed **119 tests**, Ruff lint/format, dependency checks and
sdist/wheel builds. All **119 tests also passed inside the rc4 image** as non-root,
with networking disabled, read-only source inputs and container filesystem,
and private writable output/temporary directories. Tests exercise installed
MultiQC, strict HTML rendering, JSON/text exports, all seven public profiles,
native sequencing parsers, zero/missing values, aliases and actual flat-image
export. A first sealed-container attempt timed out in the flat renderer with
`/tmp` read-only; the unchanged image passed with a private writable mount.

Eight copied-input reports rendered with the legacy engine and rc4: seven
sequencing reports and one Parse report. All section sets and content gates
passed. **20,594 common exported cells**, 1,170 plotted sample totals, 378
unknown-barcode plot values and 69 Parse metrics across three library summaries
were checked. BCL Convert lane totals agree with independent CSV sums. Quality
means use weighted averages rather than the legacy summed averages. InterOp
unit changes preserve scale and labels. Three representative reports passed
offline browser checks, including the counts-per-lane tab.

Earlier synthetic and rc1/rc2 validation remains historical evidence. The rc4
checks use the current built artifact; they are not inferred from earlier tags.
Original reference-corpus checksums remain unchanged. Full-run Parse HTML and
reduced Parse CSVs from different datasets are not used as matched numeric oracles.

Eight additional copied project-input sets passed with rc4, covering all active
project profiles. All 35 FastQC archive identities, 630 basic/status fields and
171 overrepresentation table values agree with source inputs. General Statistics
and non-FastQC exports/plots remain unchanged from rc2. Eleven reports passed
offline browser checks across the project and earlier sequencing/Parse sets.

The rc2 seven-fixture scientific/operational run passed, but the subsequent
project comparisons identified the FastQC mate loss described above. A fresh
seven-fixture run with final rc4 is being started for complete acceptance. Server acceptance includes scientific outputs, reporting, archives,
checksums, snapshots, restart behavior and configured test SMTP, with report
content checked separately against source inputs and the legacy engine.
Final results will be recorded in the coordinated image validation receipt.

## Review and adoption

| Component | Review | Responsibility |
| --- | --- | --- |
| Extension | [PR #1](https://github.com/gcfntnu/gcf-multiqc/pull/1), [issue #2](https://github.com/gcfntnu/gcf-multiqc/issues/2) | Parsers, compatibility adapters and rendering regressions. |
| Reporting image | [gcf-docker PR #1](https://github.com/gcfntnu/gcf-docker/pull/1) | Pinned recipe, publication and artifact validation. |
| Workflow selection | [gcf-workflows PR #222](https://github.com/gcfntnu/gcf-workflows/pull/222) | Authoritative image mapping and unchanged report invocation contracts. |
| BFQ boundary | [gcf-bfq PR #174](https://github.com/gcfntnu/gcf-bfq/pull/174) | Installed wrapper checks and separate-pass sequencing compatibility. |

The PRs remain review branches. Production promotion and deployment are separate
steps. An installed workflow update does not rewrite retained analysis workdirs
or override curated project image settings; evaluate their actual configuration.
A BFQ rerun from reporting preserves existing analysis reports and therefore
does not by itself demonstrate that project reports used the new engine.

Acceptance prioritizes Parse, current Cell Ranger count and active report
families. QIIME2 and legacy mkfastq are deprecated for this migration; STARsolo
has low priority. Their existing adapters remain installed. Passing the agreed
fixtures does not establish every historical producer version or optional
scientific route. SMTP acceptance is distinct from independently observed
mailbox receipt.

After reviewing the evidence, merge the coordinated source changes, then handle
production promotion and deployment explicitly. Retain the old image and input
references for rollback. A stable plugin/image release can follow that review.
The target is to retire the fork as an active development dependency.

## Maintenance

- Pin upstream, plugin source and dependencies for each reporting image. Run
  actual report rendering and export checks before accepting an upstream bump.
- Add synthetic regressions for new producer formats and retain private original
  examples outside Git. Compare source metrics and visible content, not HTML bytes.
- Use list-form overrides for list-valued `sp` keys such as `bracken` on MultiQC
  1.35; a single dictionary does not reliably replace a list default.
- Keep compatibility adapters while current producers depend on them. Simplify
  them through separately reviewed producer/report changes.
- Submit generally useful fixes upstream when worthwhile; upstream PR acceptance
  is not a prerequisite for GCF to use its independently versioned extension.
