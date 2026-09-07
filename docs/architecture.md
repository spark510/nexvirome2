# Module boundaries and object composition

The refactor preserves the CLI, output schemas, conservative inference rules and existing
Python entry points. Objects manage execution state, resource ownership and replaceable
behavior. Numeric and sequence calculations remain functions.

| Module | Responsibility | Principal objects/functions |
|---|---|---|
| `domain/` | Coordinates, evidence, decisions and correction settings | `SequenceInterval`, `CutInterval`, `BreakpointProposal`, `SearchHit`, `PairSupport`, `CorrectionDecision`, `CorrectionSettings`, `LocalPaths` |
| `io/` | FASTA/TSV/JSON formats and SAM interpretation | Format adapters, `measure_sam` |
| `runtime/` | Execution lifecycle, commands, artifact discovery | `RunContext`, `CommandRunner` |
| `coverage/` | Matrices, biological sample grouping, composition/window math | `CoverageMatrix`, `SampleSet` |
| `search/` | Backends, SQLite lifetime, search orchestration | `SearchBackend`, `MMseqsBackend`, `FoldseekBackend`, `SearchCache`, `SearchService` |
| `evidence/` | Proposals and structural evidence interpretation | `proposals`, `annotate_orfs`, `summarize_contigs` |
| `correction/` | Local paths, pair support, policy and splitting | `CorrectionPolicy`, `EvidencePolicy`, `CorrectionService` |
| `latent/` | Numerical factorization and output orchestration | `model.factorize`, `service.fit` |
| `evaluation/` | Classification, gap-aware metrics and output orchestration | `classification`, `metrics`, `service` |
| `application/` | Shared search/correction use cases and ablation settings | Request dataclasses, `run_search`, `run_correction`, `ablation_settings` |

## Ownership and composition

`RunContext` owns one fresh output directory and its manifest, including setup/input-checksum
failure records. `CommandRunner` accepts a process executor and preserves command/log/resource
fields. Services accept runners so tests can substitute execution without patching globals.

`SearchCache` owns its SQLite connection and transaction. A missing key (`None`) differs from
a cached no-hit result (`[]`). Connections close on success or error and writes roll back on
an exception. Cache keys retain sequence, settings, tool version and database/model identity.
`SearchBackend` provides tool-specific validation and command arguments through a Protocol.

`CorrectionService` composes local graph construction, mapping, evidence aggregation and a
`CorrectionPolicy`. `EvidencePolicy` retains the existing read thresholds and proposal gates.
`CorrectionDecision` rejects a proposal-created split when the read decision did not permit
one. `CorrectionRequest` has no truth/source input. This is an architectural boundary, not
an access-control mechanism.

## Coordinates and matrices

`SequenceInterval` is half-open `[start, end)` and has positive length. `CutInterval` is an
inclusive range of possible cut positions and permits an exact point. Proposal adapters
preserve extra TSV fields while validating core values. `SearchHit` uses identity/coverage
fractions and remains distinct from BLAST percentage-identity records.

`CoverageMatrix` validates shape, window/sample identity, coordinates and finite nonnegative
values. `SampleSet` groups technical replicates before biological comparisons. Feature
extraction, NMF and binning depend on this shared layer; features and NMF no longer import
each other. Generic proposals and correction policy do not depend on structural annotation.

## Compatibility

CLI search/correction and the extension workflow share application request functions. The
five ablation settings are defined once in `application/analysis.py`. Other stages retain
their existing callable application functions. Original `common`, `features`, `searches`,
`structural`, `correction`, `evaluation` and `latent` import paths remain available through
facades/re-exports. Some stage orchestration remains in facades; not every stage is a class.

The existing `Graph` object retains graph state and sequence spelling. Reverse complements,
interval operations, cosine similarity, factorization, classification and splitting remain
independently testable functions. No new Conda dependencies were introduced.

## Verification

- 65 tests passed: the previous 57 plus eight object/composition contract tests.
- 55 scientific result files are byte-identical to the pre-refactor test outputs, including
  correction FASTA/maps/evidence, NMF outputs, structural proposals, bins and metrics.
- All three Snakemake workflows passed dry-run checks.
- Invalid inputs can fail earlier at typed boundaries. Valid-input thresholds and biological
  algorithms were preserved; manifests naturally differ in code checksums and timing.

The local comparison is recorded in `.validation/refactor-output-comparison.json`. Original
sources were retained in `.validation/refactor-before` during extraction. These checks do
not substitute for actual Linux bioinformatics tool execution.
