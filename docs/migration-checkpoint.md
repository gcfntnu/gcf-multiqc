# First implementation checkpoint — 2 November 2026

## Decision and scope

Use unmodified upstream MultiQC plus the independently versioned `gcf-multiqc`
extension. Keep the old fork available during acceptance. This candidate moves
the GCF reporting code into its own repository and container, without changing
the suite's image selection, workflow code, operational configuration or data.

Package: `gcf-multiqc==0.1.0rc1`. Engine: `multiqc==1.35`. Proposed image:
`gcfntnu/gcf-multiqc:1.35-gcf0.1.0-rc1`.

The extension belongs in the reporting container. BFQ continues to orchestrate
reporting, and gcf-workflows continues to own the eventual image selection.
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

This establishes a runnable candidate against synthetic/public examples. It
does **not** establish report parity for every historical producer version or
validate scientific outputs, Apptainer behavior on the facility server, BFQ
orchestration, real project reports, notification delivery or deployment.
No GCF-suite integration test is included in this checkpoint.

The source location referred to as `src` has not yet been identified. Until a
representative operational log corpus is supplied, real-data acceptance remains
open. Copy approved logs into an ignored invocation-owned directory before
testing; never add project data to this repository or the container context.

## Next acceptance step, when authorized

1. Copy one representative report-input set for each active workflow family,
   including distinct Parse layouts, older 10x mkfastq and current count CSVs.
   Include a zero/low-yield case and sample names with leading zeros or suffixes.
2. Render the old and candidate reports from those same copied inputs with the
   existing GCF YAML and module selections. Run reporting only; do not rerun
   demultiplexing, alignment or analysis.
3. Compare sample sets, required sections, raw metrics, units, missing values,
   links, metadata joins and PCA group/color behavior. Classify the documented
   corrections explicitly instead of demanding equality for known old defects.
4. Validate the image through the existing Apptainer wrapper on a non-production
   server using copied inputs. Record image digest and source revisions.
5. Only after acceptance, propose a separate gcf-workflows image-reference
   change. Keep the old image and report inputs available for rollback. Promote
   a stable plugin/image version after that review.

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
