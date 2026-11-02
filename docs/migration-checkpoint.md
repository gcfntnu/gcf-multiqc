# Implementation and suite-readiness checkpoint — 2 November 2026

## Decision and scope

Use unmodified upstream MultiQC plus the independently versioned `gcf-multiqc`
extension. Keep the old fork available during acceptance. This candidate moves
the GCF reporting code into its own repository and container. The follow-on
suite PRs prepare adoption for review; operational configuration, data and the
existing fork remain unchanged.

Package: `gcf-multiqc==0.1.0rc1`. Engine: `multiqc==1.35`. Published candidate:
`gcfntnu/gcf-multiqc:1.35-gcf0.1.0-rc1`, digest
`sha256:5b69b8bdb2b8d055b190386110c55dcddc99e1e30b6ee846477a703b5c8dd39d`.
Its plugin runtime source is `54877b7335635005af113a1a0cd10a9f7fbdfe12` and its
Plotly version is 6.0.1. Later plugin changes add tests and documentation without
changing runtime code or dependencies. The suite-readiness work does not build
or publish another image.

The extension belongs in the reporting container. BFQ continues to orchestrate
reporting, and gcf-workflows continues to own image selection.
Do not install this package into BFQ's command environment.

## Source evidence

| Source | Reviewed baseline |
| --- | --- |
| gcfntnu/MultiQC, `bfq-dev` | `d3fce9dab5b8522ad022f04161a145db121ad23d` |
| Upstream MultiQC, v1.35 | `87e504bb77e3d94687dd980638724b7ef8e1e92f` |
| gcfntnu/gcf-workflows, `bfq-dev` | `26fafe406d856b49beead1ecce251cc943e56e63` |
| gcfntnu/gcf-bfq, `bfq-dev` | `2b1fbc51f625c128e9aa5a083f1c5fe1abb06068` |
| gcfntnu/gcf-docker, original recipe base | `28e66251ec003b0f0d65526d2fa6b8cc03ed92b3` |

The seven copied public workflow templates are unchanged. Tests replace their
logo location and add synthetic metadata. Native sequencing fixtures were
derived from BFQ's existing synthetic tests; their README records the exact
local source revision. No BFQ or scientific workflow was executed.

## What is implemented

The package retains the module selections used by RNA-seq, microbiome,
single-cell, small-RNA, default, WGS and metagenome reports. GCF-specific parsers
are small standalone modules. FastQC, fastp, STAR and Bracken reuse upstream
parsing where possible. The custom-content adapter preserves literal sample
identities and adapts existing PCA/UMAP input conventions without rewriting
producer files or recalculating scientific results.

Testing the existing configurations identified an additional compatibility
requirement: upstream removed the separate `bracken` selector after the old
fork's baseline. The extension restores that selector and keeps its estimates
separate from Kraken observations, including when both use the same sample ID.

The README specifies intentional corrections to Unitas denominators, STAR
splice-event semantics, malformed ordination handling and duplicate identities.
The old template and inactive `fastv` parser were not copied.

Flat rendering exposed an upstream dependency-range defect: MultiQC 1.35 pins
Kaleido 0.2.1 but permits Plotly releases that dropped that exporter. The plugin
adds the same `plotly<6.1` constraint found in subsequent upstream source, and a
real flat-image regression check. This protects large reports as well as
explicit `--flat` invocations.

## Validation and limits

The automated checks exercise installed MultiQC, actual HTML generation under
`--strict`, numeric JSON/text exports and unchanged public workflow profiles.
They cover every newly ported module in its relevant profile, native BCL Convert,
bcl2fastq and InterOp, explicit compatibility disablement, alias coexistence,
literal sample identities, zero totals and malformed inputs. Tests block
network and SMTP access. Package and container checks use isolated paths.

Additional suite-readiness checks use only synthetic data. They render with
the workflow's full relative `--filename` and separate logs/figures input roots,
preserve custom metadata joins, and exercise BFQ's sequencing configuration with
BCL Convert, bcl2fastq or the combined custom report for separate conversion
passes. The last case exposed a repeated `custom_content` entry in BFQ's
`module_order`, which repeats report sections under MultiQC 1.35. BFQ issue #173
removes that duplicate; the plugin
regression uses the corrected selection. A Parse regression covers SplitPipe's
row-oriented producer headers, genome-prefixed extra metrics and negative slope
values using invented values. These checks require no plugin runtime changes.

The suite-readiness run passed all 107 tests in an isolated Python 3.11
environment with MultiQC 1.35 and Plotly 6.0.1. Ruff lint, Ruff format checks,
`pip check` and `git diff --check` passed. This phase used no container execution,
image build/publication or operational report rendering.

This establishes a runnable candidate against synthetic/public examples. It
does **not** establish report parity for every historical producer version or
validate scientific outputs, Apptainer behavior on the facility server, BFQ
orchestration, real project reports, notification delivery or deployment.
No GCF-suite integration test is included in this checkpoint.

Existing reports and approved inputs have been copied into a private local
reference inventory for later acceptance. They are not committed, included in
the container context or rendered by the suite-readiness checks. Header-only
inspection of copied Parse summaries informed the synthetic producer-format
test. The reduced Parse summaries and the original-run HTML represent different
datasets; they must not be compared as a matched input/expected-output pair.
Real-data acceptance remains open. Keep private reference material under ignored
invocation-owned directories outside the build context.

## Review handoff and acceptance

Prepare and review the suite PRs before running operational acceptance:

| Work | Tracking issue | Review boundary |
| --- | --- | --- |
| Plugin regressions and acceptance contract | [gcf-multiqc #2](https://github.com/gcfntnu/gcf-multiqc/issues/2) | Extend [implementation PR #1](https://github.com/gcfntnu/gcf-multiqc/pull/1); preserve the candidate runtime. |
| Workflow image selection | [gcf-workflows #221](https://github.com/gcfntnu/gcf-workflows/issues/221) | Pin the candidate digest; retain existing producer/report contracts. |
| BFQ wrapper and sequencing compatibility | [gcf-bfq #173](https://github.com/gcfntnu/gcf-bfq/issues/173) | Validate the installed wrapper and avoid duplicate custom-content sections for separate conversion passes. |
| Container recipe and rebuild handoff | [gcf-docker #2](https://github.com/gcfntnu/gcf-docker/issues/2) | Record the pinned source/dependencies and reviewable build instructions. |

These are implementation PRs, not approval to merge, promote, deploy or run the
scientific workflows. Operational acceptance remains a later authorized step:

1. Select one matched report-input/reference set for each active workflow family,
   prioritizing Parse, current Cell Ranger count and sequencing reports.
   Include a zero/low-yield case and sample names with leading zeros or suffixes.
   QIIME2 and legacy mkfastq are deprecated for this migration and do not block
   adoption; retain their implemented adapters. STARsolo is lower priority
   because BFQ no longer runs it by default.
2. Render the old and candidate reports from those same copied inputs with the
   existing GCF YAML and module selections. Run reporting only; do not rerun
   demultiplexing, alignment or analysis.
3. Compare sample sets, required sections, raw metrics, units, missing values,
   links, metadata joins and PCA group/color behavior. Classify the documented
   corrections explicitly instead of demanding equality for known old defects.
4. Validate the image through the existing Apptainer wrapper on a non-production
   server using copied inputs. Record image digest and source revisions.
5. After acceptance and explicit merge/deployment decisions, merge the reviewed
   adoption PRs in their documented dependency order. Keep the old image and
   report inputs available for rollback. Promote a stable plugin/image version
   after that review.

The target is to retire the fork as an active development dependency. Keep the
compatibility adapters while existing producers depend on them; simplify them
only through separately reviewed changes to producer/report contracts.

## Maintenance

- Pin the upstream engine and the plugin source revision in the image recipe.
  Pin runtime dependencies in the Docker build directory and record the base
  image digest. Do not build from a moving fork branch.
- For an upstream bump, first run the compatibility suite and inspect report
  exports; then publish another candidate for the acceptance corpus. A passing
  import or version command is insufficient.
- Use list-form overrides for list-valued `sp` keys (for example `bracken`) on
  MultiQC 1.35; its configuration merge does not reliably replace a list default
  with a single dictionary.
- Add regression fixtures when a new producer format is introduced. Prefer
  synthetic reductions; retain private original examples outside Git.
- No plugin has to become an upstream PR before GCF can use it. Submit broadly
  useful changes upstream independently when their maintenance benefit warrants
  that work.
