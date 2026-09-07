# Reference masking and evidence-driven execution

Implemented on 2026-09-07. These modules are executable research infrastructure;
they do not establish a benefit on natural RefSeq or reproduce an unknown legacy
NexVirome classifier. No new dependencies were added to the Conda environments.

## Installation and runnable fixture

From the checkout root on Linux:

```bash
conda env create -f environment.yml
conda activate nexvirome2
snakemake -s workflows/masking.smk --configfile tests/fixtures/masking.workflow.yaml --cores 1
```

For the production-sized research inputs, replace the paths in
`configs/masking_workflow.yaml` and use that config file instead. Use a fresh
`output` directory for a new experiment. Empty scheduler-created directories are
accepted; any existing artifact prevents a stage from overwriting the directory.
Failed stages retain their manifest and logs and require a new run path.

The fixture is deliberately tiny and synthetic. `A.1`/`B.1` and taxids in that
fixture are test identifiers, not a downloaded natural RefSeq panel.

## Inputs and commands

| Command | Input contract | Main output |
|---|---|---|
| `taxonomy build` | `--mapping` accession/taxid TSV, `--nodes` local nodes.dmp; optional merged/deleted taxdump | frozen rank taxids; unknown ranks remain empty |
| `sharing build` | `--reference`, `--taxonomy`; configured rank and k | positional exact canonical seed atlas, SQLite index |
| `sharing align` | reference, taxonomy; optional 11-column BLAST `--alignments` | approximate shared-region atlas |
| `mask import` | `--reference`, `--source`; format intervals/hard/soft | checksum-bound mask.json |
| `mask propose` | `--reference`, `--atlas` | shared-seed mask and matched random control |
| `mask export` | `--reference`, `--mask`; mode hard/soft/unmasked | FASTA view, intervals TSV, mask.json |
| `classification run` | `--queries` FASTA, reference, taxonomy; optional mask OR competitive hits | assignments plus all candidate hits |
| `classification evaluate` | predictions, classification metadata, truth TSV | all-query and independent-unit metrics |
| `mask evaluate` | same queries/reference/taxonomy/truth, arms JSON | paired arm comparison and raw metrics |
| `evidence domains` | original reference, ORFs, amino-acid evidence intervals | strand-aware nucleotide coordinates |
| `evidence link` | reference, mask, annotation intervals | mask/annotation overlaps and union lengths |
| `evidence joint` | reference, disjoint classification/assembly observations | 2x2 base counts and separate unknowns |
| `discovery summarize` | ORF profiles; optional motif/context assay results | separate sequence/fold/motif/context evidence cards |
| `catalog build` | source/FASTA manifest | exact strand-invariant dereplication and source lineage |
| `catalog split` | accession/cluster/genome_group membership | transitive group-disjoint development/heldout split |
| `design multisample` | existing sources JSON; biological sample count and technical seeds | per-sample source coverage and simulation specs |
| `paths rank` | normalized candidate/path/feature/value TSV and weights | weighted baseline ranking; never authorizes a split |
| `quality genomad` | original contigs and native virus_summary.tsv | original-coordinate candidate/QC table |
| `pipeline freeze/select/execute` | development units, paired evidence, frozen implementations, input recipes | adoption decisions and permitted module execution |

All commands require `--config` and `--output`. CLI help lists required paths.
Paths in catalog manifests, source specifications and arms JSON resolve relative
to the file containing them. Direct CLI and production recipe paths resolve from
the current working directory. No Windows absolute path is embedded in code.

## Mask contract and limits

The original FASTA checksum, accession IDs and lengths must match. The canonical
interval TSV has `accession,start,end,event_id,evidence`; coordinates are 0-based
half-open. Optional event/evidence fields get explicit import labels. A circular
wrap must be imported as two valid intervals with the same event ID. Whole
records are never concatenated. Overlap is merged for rendering and counting;
the original event rows remain in mask.json.

Policies require `version`, `rank`, and `conditions`. A hard-mask FASTA cannot
reveal whether an original N was deliberately masked: that status is recorded as
unknown. Soft-mask import reads case directly instead of the uppercase-only
general FASTA parser. Exporting an unmasked view restores the original nucleotide
sequence, but does not reproduce original header descriptions or FASTA wrapping.

The atlas counts distinct taxids rather than duplicate accessions. Unknown-rank
competitors prevent automatic masking proposals. Low-complexity seeds and
within-reference repetition are separate annotations. Current proposals mask
the union of qualifying seed spans; that is a candidate policy, not a calibrated
read-length-aware policy. Random controls preserve union interval widths and total
masked bases per accession, including separation between distinct intervals.

With `sharing build --metadata`, accession/topology TSV records declared circular
generate junction seeds. The two original-coordinate parts share an evidence ID.
Records shorter than k have no seeds. The masking workflow accepts an optional
`metadata` path for this purpose. Unknown topology defaults to linear indexing.

`sharing align` imports the existing 11-column BLAST format, or runs an all-reference
competitive BLAST search when `--alignments` is omitted. Thresholds are
`min_identity` (default 99) and `min_alignment_bp` (150). It validates alignment
strings against both original references, retains both sides, splits at gaps,
and reports distinct competing taxa along intervals. Approximate atlas files can
be supplied to `mask propose`. BLAST word-size sensitivity remains; exhaustive
homology detection, circular BLAST coordinate projection, and empirical VDS
reconstruction are not claimed. Configured conditions do not turn sharing counts
into a validated discriminability score.

## Classification and evaluation semantics

The built-in backend is explicitly `diagnostic_exact_seed`: it counts distinct
query seeds matching each reference, retains competing best-score references,
then resolves them at the requested rank. Unknown taxonomy never becomes a
confident assignment. It is not an implementation of Kraken, BLAST classification,
or the user's existing NexVirome. Hard-masked bases exclude seeds. External
classifiers can provide a headered `query,target,score` competitive-hit TSV with
one aggregated score per query/target; duplicate entries take their maximum,
not their sum. Record external method/version/view provenance in config.

Queries are FASTA assay units, not directly paired FASTQ input. For an external
fragment classifier, aggregate mates before providing the hit table. Do not
concatenate mates and interpret the artificial junction as a real sequence.
All query IDs must appear exactly once in prediction and truth tables, including
unclassified queries. Truth columns are `query,taxon,unit`; taxon is the answer
at the declared rank and unit is an independent source/pair/sample group.

Correct, wrong, ambiguous and unclassified fractions use **all** queries.
Precision uses assigned queries and is null when none are assigned. Paired
bootstrap intervals resample independent units, not reads. Small unit counts
remain exploratory. No biological performance claim is generated by the masking
fixture workflow. Reads must originate from the unmasked sequence; the pipeline
never generates reads from masked references.

## Conserved-region and graph evidence

`evidence domains` accepts `orf_id,aa_start,aa_end,evidence_id,kind,label`, including
normalized HMM/domain/motif evidence. Amino-acid intervals are 0-based half-open;
ORF `contig,start,end,strand` determines nucleotide projection. Imported evidence
must use the correct genetic code and exclude stop residues. `evidence link`
preserves complete annotation intervals while reporting their overlap with masks;
an overlap by itself does not prove function or classification benefit.

Joint observations use `unit,accession,start,end,classification,assembly`.
Allowed classification calls are ambiguous/informative/unknown, and assembly
calls chimeric/consistent/unknown. Unknowns do not enter the 2x2 denominator.
Overlapping intervals within one unit/accession are rejected. This is an explicit
coordinate adapter, not automatic multi-hit coordinate projection or causality.

Discovery distinguishes no-hit from not-assessed. Set `sequence_assessed` and
`fold_assessed` only when searches covered the supplied ORFs. Motif/context tables
have `orf_id,status` with supported/unsupported/uncertain. No result is proof of
a novel family. The catalog deduplicates exact sequence/reverse-complement only;
near-identity and circular-rotation clustering are not implemented.

Every usable graph correction now saves `fragment_candidates.json` and a ranked
`path_support.json`. Ambiguous pairs retain alternatives but never count toward
unique support. Weighted auxiliary ranking renormalizes over observed features,
records missing weight, and cannot create or force a graph connection.

## Expensive analysis and native QC

Set `representation_cache` in structure settings to enable a database-independent
Foldseek query database cache for ProstT5 searches. Cache identity includes the
query batch, model contents, tool version and prediction threads. Repeated target
DB searches reuse that batch representation; different query batches may still
repeat prediction. Checksums detect corrupt/incomplete caches. This does not yet
cache each individual protein across arbitrary overlapping batches.

The opt-in implementation uses `createdb --prostt5-model`, then `search` and
`convertalis`, following [Foldseek's documentation](https://github.com/steineggerlab/foldseek)
and [conversion example](https://github.com/steineggerlab/foldseek/wiki).
Only command composition/cache behavior has been tested here with a fake runner;
actual Foldseek/ProstT5 execution remains a Linux validation requirement.

The geNomad adapter follows its [native output documentation](https://github.com/apcamargo/genomad#understanding-the-outputs).
It converts provirus coordinates from 1-based inclusive to original 0-based
half-open coordinates, preserves host-contig linkage and does not call genuine
provirus integration a chimera. A missing viral-summary row is not a negative
classification; virus scores are not completeness. CheckV native integration
and cross-tool applicability checks remain outstanding.

## Freezing and operational selection

`configs/adoption.example.yaml` is an example, not an accepted policy. Freeze a
copy with `pipeline freeze --development development.tsv`; its `parameters`
include each module's implementation command/settings, benefit threshold,
maximum recovery loss, resource limit and conditional prerequisites.

Evidence TSV columns:

```text
module unit baseline_error error baseline_recovery recovery resource_ratio
```

The delimiter is TAB. Error and recovery are fractions per independent unit;
resource ratio is positive. Use paired baseline/intervention observations under
the same conditions. Do not put raw wrong-connection counts in a rate column.
`pipeline select` checks parameter hashes, duplicate units, development/heldout
unit leakage, baseline errors, paired uncertainty and resource constraints.
`evidence_kind: natural_heldout` is an explicit caller declaration and cannot
authenticate the biological origin of arbitrary imported measurements.

Outputs distinguish `default`, `conditional`, `research`, and `excluded`.
No baseline errors means reduction is unestablished. Unproven effects remain
research; harmful effects or exceeded resource limits are excluded. Operational
implementation settings must be frozen before a default/conditional decision.

`pipeline execute` takes the production JSON, its decision JSON and a recipes JSON:

```json
{"stages":[{"module":"reference_masking","arguments":{
  "reference":"data/refseq/reference.fasta","mask":"data/policy/mask.json"}}]}
```

It executes only approved operations with the frozen configuration, checks
decision and implementation hashes, and records skipped conditional modules.
Deployment facts are supplied explicitly in `settings.conditions`; absent facts
are not inferred as true. Recipes can declare `depends_on` and refer to a completed
predecessor as `${module}/artifact`. Recipe paths are variable data inputs; they
are not evidence that the biological applicability of a policy has been validated.
Simulation/evaluation operations and config/seed/thread overrides are excluded
from this operational interface. No operation is automatically adopted in this repo.

## Remaining implementation and empirical work

The integrated plan is **not wholly complete**. In particular, legacy classifier
and mask semantics require the actual NexVirome source and example artifacts;
natural full-snapshot experiments and locked tool execution require the intended
Linux environment/data. [Additional extensions](remaining_extensions.md) now provide
development-only condition-specific policy selection, masked mixed-feature NMF,
conserved-support reuse, read QC/imported host-call filtering, near-identity catalog
clustering and selected coassembly pools/workflows. These have documented scope limits.
Automatic coordinate joins, automatic host/ERV mapping, CheckV integration and comprehensive
resource reporting remain outstanding. Lost graph paths, new assemblers, long-read/hybrid and learned
graph models remain evidence-dependent later research, not default modules.
