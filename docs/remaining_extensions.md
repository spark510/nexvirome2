# Additional implementation — 2026-09-07

These commands extend the [masking pipeline](masking_pipeline.md). All work is in
NexVirome2. The original NexVirome implementation and real Linux benchmark inputs
are still required for equivalence and natural-data performance validation.

Mixed-feature NMF uses NumPy/SciPy from the existing extended Conda environment
(`environments/extended.yaml`); no new dependency was introduced. Run its commands
inside that environment. The other new adapters use the existing core dependencies.

## Condition-specific mask selection

`mask optimize --reference REF --candidates candidates.json --measurements measurements.tsv
--config configs/mask_optimization.yaml --output NEW_OUTPUT` selects among whole
candidate mask policies on **development data only**. It does not learn arbitrary
new mask boundaries. Every candidate at the requested rank must be measured on
identical independent units and query fingerprints. Reducing wrong assignments
cannot compensate for excessive loss of correct assignments or original-coordinate
retained breadth. Selection uses paired bootstrap intervals over independent units.
Candidate selection is development tuning, not heldout evidence or production approval.

Candidates JSON:

```json
{
  "conditions": {
    "pe150": {
      "rank": "species", "read_length": 150, "library": "genome-derived",
      "database_sha256": "SHA256_OF_ORIGINAL_REFERENCE"
    }
  },
  "policies": {"unmasked": null, "shared": "shared/mask.json"}
}
```

Measurements are headered TSV with
`condition,policy,unit,queries_sha256,wrong,correct,unclassified,retained_breadth`.
Fractions use all queries. Retained breadth is a separate measured quantity on
original coordinates; do not substitute precision or the assigned-read fraction.
Paths in candidates JSON resolve relative to that JSON. Outputs include
`selection.json`, condition-bound mask files, candidate effect intervals and
`development_units.tsv` for the existing policy-freezing interface.

## Mixed-feature NMF

```bash
nexvirome2 multiview fit --features features.tsv --config configs/multiview.yaml --output runs/multiview/fit
nexvirome2 multiview propose --model runs/multiview/fit --coordinates coordinates.tsv --config configs/multiview.yaml --output runs/multiview/proposals
```

Features TSV is long format: `entity,block,feature,value`. Example blocks are
coverage, composition, domains and structural annotations. Values must be finite
and nonnegative; empty/NA cells are missing measurements, distinct from measured
zero. Feature extraction and technical-library aggregation must happen upstream.
This generic feature model does not certify biological sample independence.

Variable columns are scaled by their observed 95th percentile and clipped to
[0,1]. Constant/unobserved columns cannot contribute. Each block's configured
weight is divided across its observed variable features per entity. Weighted
multiplicative updates exclude missing cells from the reconstruction objective;
they do not fit missing entries as zeros. Multiple initializations are aligned
and compared for loading stability. Components are not viral identities or
calibrated probabilities. Normalization is currently fitted for each analysis;
transfer of a frozen normalization/model to new cohorts remains separate work.

Coordinates TSV is `entity,contig,start,end` on the original target, 0-based
half-open. It must cover all entities without overlapping windows. Only stable,
confident component changes between contiguous windows produce proposals.
All-zero informative feature rows are retained as `inactive` entities with zero
loadings and excluded from model fitting. Coordinates are retained, and proposals
never bridge an inactive window. `loadings.tsv` includes an explicit status column.
The `multiview` proposal is an auxiliary gate: existing graph/read support remains
necessary to split a contig.

Set `multiview_features` and `multiview_coordinates` in the extensions workflow
config to add fitting, proposals, correction and a matched random-control arm.
With `truth_sources`, all six correction variants and six random controls are
evaluated alongside the original contigs. Existing five-variant behavior remains
the default when no mixed-feature input is configured.

## Read QC, host evidence and conserved support

`reads prepare --r1 READS1 --r2 READS2 --config configs/read_preparation.yaml
--output NEW_OUTPUT` validates paired four-line Phred+33 FASTQ and applies pair-level
length, mean-quality and N-content filters. Both mates stay together. Retained,
quality-filtered and optionally host-removed pairs are all written as paired gzip
FASTQs; `pairs.tsv` records reasons and original IDs. Reads are not trimmed.

Host subtraction is off by default. To enable it, provide `--host-calls` TSV
(`fragment,status`) and `--host-metadata` JSON containing matching read SHA256s,
method and reference SHA256; set `remove_host: true`. Only `host_supported`
fragments are removed for host evidence. Ambiguous, nonhost, unclassified and
unassessed calls are not removed for that reason, although ordinary quality
filters still apply. Calls must come from independently generated competitive
evidence. This is an import adapter, not automatic host/viral/ERV mapping, and
cannot independently certify the truthfulness or quality of imported evidence.

`conserved reuse --masked MASKED_RUN --unmasked UNMASKED_RUN --config CONFIG
--output NEW_OUTPUT` compares runs with matching queries, taxonomy, reference,
conditions and classifier settings. It checks candidate/prediction hashes and
retains full-reference candidate support, including competitors outside masked
evidence. Fine assignments remain exactly those from the masked run. A conserved
read can support later annotation/detection without becoming a forced fine-rank
assignment. This adapter does not estimate viral origin or recovered genome breadth.

## Near-identity catalog and targeted coassembly

`catalog cluster --contigs CONTIGS --config CONFIG --output NEW_OUTPUT` uses
competitive BLAST or optional precomputed `--alignments` in the existing 11-column
format. Defaults are identity 95% and mutual coverage 95%. Imported alignment
strings must match both original sequences. Identity is recomputed; gap positions
do not count as recovered coverage. Each member needs one qualifying HSP directly
to its representative. No transitive similarity chain, multi-HSP reconstruction
or circular rotation equivalence is assumed. Output retains representative/member
orientation and alignment criteria.

`coassembly prepare --samples samples.tsv --membership membership.tsv --config CONFIG
--output NEW_OUTPUT` pools explicitly selected fragments only. Samples TSV requires
`sample,biological_sample,library,r1,r2`; membership TSV requires
`group,sample,fragment`. Read paths resolve relative to samples TSV. Groups cannot
mix library labels. Sample prefixes prevent read-ID collisions, and lineage TSV
maps every pooled fragment to its original sample/read. Missing selected reads,
duplicate IDs and mismatched pairs fail the stage. Multiple technical libraries
do not increase the count of biological samples.

`workflows/coassembly.smk` takes config keys `samples`, `membership`, `output`,
optional `assemblers` and `threads`. It runs both pooled and per-sample assemblies
using existing metaSPAdes/MEGAHIT wrappers. Its comparison is an assembly-size
summary, not accuracy evidence, and whole per-sample inputs do not constitute a
matched read-budget control. Frozen independent evaluation remains necessary.
`tests/fixtures/coassembly.dryrun.yaml` has deliberately invalid placeholder reads
and is for DAG validation only; actual FASTQ pooling is covered by tests.

## What still requires work

These additions do not complete every research requirement. CheckV native result
integration, automated host/ERV competing alignments, automatic multi-hit coordinate
joins, expanded background/library-artifact simulations and full resource reporting
remain outstanding. Empirical retained-breadth measurements and legacy-classifier
adapters require actual experiment outputs. Natural heldout experiments and a
locked Linux tool environment have not been run. New assemblers, lost-graph-path
recovery, hybrid/long-read and learned graph models remain conditional research.
