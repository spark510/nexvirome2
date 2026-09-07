# Validation record — 2026-09-06

Verified in the Windows development workspace using Python 3.12.12:

- 26 pytest tests passed, covering competitive-alignment calls, graph orientation and
  overlaps, fragment deduplication/multimapping, lossless splitting, mocked end-to-end
  correction and ART coordinate handling, snapshot reconciliation, deterministic panel
  selection, development/evaluation separation, and zero-baseline report behavior.
- Both Snakemake workflows passed dry-run validation with Snakemake 9.26.1 and the supplied
  dry-run fixtures. Reference preparation has four jobs including the aggregate target.
- A platform-independent Python wheel built successfully.
- A live EFetch smoke test downloaded NC_001802.1 (9,181 bp), parsed its metadata, and
  passed accession and checksum validation. This one-record smoke dataset is stored in
  the ignored `.validation/refseq-smoke` directory; it is not the study snapshot.

Not yet executed at that checkpoint: a Linux Conda environment solve/lock, actual ART/SPAdes/MEGAHIT/Bowtie2
integration, a complete NCBI Virus RefSeq snapshot, or the natural/held-out benchmarks.
These require the chosen Linux server and the frozen full accession list. Passing tests
is not evidence that the correction improves biological assembly accuracy.

## Conda installation follow-up — 2026-09-06

- Added root `environment.yml` with automatic editable installation of the CLI and a
  `nexvirome2 doctor` startup check for all required command-line tools.
- 29 tests passed after adding missing-tool, broken-executable, and healthy-tool checks.
- Conda 24.11.3 successfully solved the environment with `--platform linux-64 --dry-run`,
  strict channel priority, Linux 5.10 and glibc 2.17 virtual-package overrides.
- The solve selected Python 3.12.14, SPAdes 4.3.0, MEGAHIT 1.2.9, ART 2016.06.05,
  BLAST 2.17.0, Bowtie2 2.5.5, SAMtools 1.24, and Snakemake 9.26.1.
- Inspected Conda's pip installer path handling: editable paths resolve from the YAML
  directory, so root `environment.yml` uses `-e .` and `environments/linux.yaml` uses `-e ..`.

This validates dependency resolution, not Linux installation or binary execution. The pip
installation phase is not executed by the dry-run. A server-side environment creation and
`nexvirome2 doctor` are still required before the real benchmark. No lockfile is represented
as an installed/validated Linux environment. NMF was still a research option at that checkpoint.

## Multi-sample and structural implementation follow-up — 2026-09-06

- 46 tests passed after implementing the optional extensions. Added checks cover fragment
  coverage and target-digest validation, technical-replicate rejection, NMF identifiability
  and synthetic abundance transitions, ORF strand/coordinate preservation and deduplication,
  search-cache invalidation, competitive structural-label uncertainty, structural lineage
  versus functional-family changes, proposal-gated lossless correction, local neighborhood
  selection, conservative bins, and mocked motif/alignment/quality command interfaces.
- The optional Snakemake workflow passed dry-run with Snakemake 9.26.1. Its local-selection
  fixture contains 16 jobs including baseline/simple/NMF/structural/combined correction.
  Fixture inputs are dependency placeholders, explicitly not real tool inputs.
- The final extended Linux Conda environment solved successfully after adding geNomad and
  CheckV. The solve selected MMseqs2 18.8cc5c, Foldseek 10.941cd33, Folddisco 2.9375a2d,
  FoldMason 4.dd3c235, geNomad 1.12.0, CheckV 1.1.1, scikit-learn 1.9.0 and Biopython 1.88.
  The ignored `.validation/conda-extended-solve.json` stores the resolved dependency list.
- A wheel containing the new Python modules built successfully at
  `.validation/wheels-extensions/nexvirome2-0.1.0-py3-none-any.whl`.

NMF and translation tests execute the actual Python libraries on synthetic data. External
search, coverage-sort, correction-mapping and quality-tool integration tests use mocked
outputs. The Conda solve does not install or start Linux binaries and is not a frozen
server lock. Full database/model downloads, actual structural searches, natural RefSeq
chimera reproduction and held-out correction performance have not been executed.
See [extensions.md](extensions.md) for reproducible installation and execution instructions.

## Bug audit follow-up — 2026-09-06

57 tests passed after reproducing and fixing five defects in gapped-alignment recovery,
TSV data integrity, stage failure provenance, missing structural evidence in binning, and
imported ORF stop handling. See [bug_audit.md](bug_audit.md) for triggers and behavioral changes.

## Modularization follow-up — 2026-09-06

65 tests passed after extracting typed domain values, runtime lifecycle, search backends and
cache, correction policies, shared matrices, numerical inference and evaluation modules.
All three workflows passed dry-run. A byte comparison of 55 scientific output files against
the pre-refactor test run found no differences. See [architecture.md](architecture.md).

## Integrated masking and operational gates — 2026-09-07

- 88 tests passed (`.test-tmp-integrated-5`, 20.22 seconds), including mask round trips,
  taxonomy ranks/redirects, circular seed coordinates, approximate shared alignments,
  all-query evaluation denominators, grouped holdout splits, independent biological
  abundance, ambiguous fragment retention, domain projection, missing evidence,
  frozen adoption/production execution, representation-cache integrity and native
  geNomad provirus coordinates.
- `workflows/masking.smk` completed all five jobs on the synthetic fixture, using
  `.validation/masking-workflow-final`. Its per-stage manifests report completion.
  This was real Snakemake execution on Windows with the Python diagnostic classifier,
  not a natural RefSeq or external-assembler performance experiment.
- Real workflow execution exposed the scheduler-created empty-directory defect in
  RunContext. Empty directories are now accepted with exclusive manifest creation;
  existing artifacts/completed/failed runs remain protected. Both cases have tests.
- `workflows/extensions.smk` passed dry-run with `truth_sources` enabled, connecting
  original contigs, five correction variants and five matched random controls to
  separate source-aware evaluations. Full external-tool execution is still pending.
- Clean source-copy wheel build passed: 70 Python modules, no obsolete latent.py.
  Wheel: `.validation/wheels-integrated/nexvirome2-0.1.0-py3-none-any.whl`.
- No additional Conda dependencies were introduced. The earlier Linux solve is
  unchanged; this does not substitute for installation/locking on the actual server.

The optional Foldseek prediction cache was tested with injected fake executables,
including reuse across two target databases. Actual model inference was not run.
No existing NexVirome source/mask artifact was supplied, so legacy equivalence is
unverified. Detailed contracts and remaining code work are in
[masking_pipeline.md](masking_pipeline.md).

## Resumed remaining extensions — 2026-09-07

- Final regression suite: **101 passed**, 30.46 seconds, `.test-tmp-resume-final`.
- New coverage includes paired development policy comparisons and recall/breadth
  loss gates, masked missing-value NMF, proposal stability, conserved candidate
  retention without forced classification, paired QC/host-call preservation,
  FASTQ identifier validation, representative-based near-identity clustering,
  selected read pooling and technical-library grouping.
- `workflows/coassembly.smk` passed dry-run with
  `tests/fixtures/coassembly.dryrun.yaml`. Actual pooling and source/read lineage
  were exercised by tests with paired FASTQ; external pooled assemblies were not run.
- `workflows/extensions.smk` passed dry-run with mixed-feature input, coordinate
  input and `truth_sources` enabled. This adds a sixth correction variant and its
  matched random control. Placeholder files in that DAG check are not scientific inputs.
- New catalog/coassembly CLI help was checked. Existing source-aware evaluation
  remains the required accuracy test; assembly-size summaries are not accuracy claims.

See [remaining_extensions.md](remaining_extensions.md) for contracts and work still
outstanding. The earlier usage-limit rejection was resolved on resumption; the
previously unapplied catalog/coassembly files are now present and tested.

## Audit fixes — 2026-09-07

All five independently reproduced bugs in the second audit are fixed. The complete
suite passed: **117 tests**, 34.40 seconds, `.test-tmp-audit-fixed`. Regression tests
cover exact CLI options and frozen configuration enforcement, actual aligned anchor
length, unavailable MAPQ, scored secondary competition and inactive NMF entities
without coordinate bridging. These are code-level tests with synthetic inputs,
not a new external-tool or natural-data performance benchmark. See
[the updated audit](bug_audit_20260907.md).

## Installation channel checks - 2026-09-07

- Full regression suite: **117 passed**, 39.64 seconds, `.test-tmp-install-channels`.
- `scripts/install.py` command construction was checked for base, extended and dry-run
  modes. It applies strict priority and a two-channel allowlist during creation, disables
  configured default packages and leaves global Conda configuration untouched.
- The extended environment solved successfully using cached repodata, `--platform linux-64`,
  `--offline --dry-run`, and Linux/glibc virtual-package overrides on the Windows development
  host. Every solved Conda dependency came from `conda-forge` or `bioconda`; the local pip
  editable install was not executed by this check. Output: `.validation/conda-channel-check.json`.
- This checks dependency resolution, not an actual Linux installation, tool startup or
  a natural-panel benchmark. Package versions are not yet locked from a verified server run.
