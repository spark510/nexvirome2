# NexVirome2 implementation and Linux runbook

This file describes currently executable functionality. The complete remaining scope,
including shared-region evaluation feeding back into existing NexVirome reference masks,
is tracked in [implementation_plan.md](implementation_plan.md). Planned CLI names there
are not available commands until implemented and validated.

This is a research prototype: metaSPAdes and MEGAHIT construct assemblies; NexVirome2
evaluates cross-source joins and conservatively splits supported repeat-exit mistakes.
It does not construct a new de Bruijn graph or restore paths removed by an assembler.

## Installation on a Linux server

Run from a Linux checkout of this repository. No Windows path is embedded in the package.
The initial defaults are four threads and 16 GiB for SPAdes; they are configurable and
are not a resource estimate for downloading/indexing the complete viral reference set.

```bash
conda env create -f environment.yml
conda activate nexvirome2
nexvirome2 doctor
python -m pytest -q
conda list --explicit > environments/linux-64.lock.txt
```

Run environment creation from the repository root: the YAML pip section installs `-e .`
automatically after the Conda packages. `environments/linux.yaml` remains an equivalent
compatibility entry point; its editable package path is relative to that YAML's directory.
This is a Linux x86_64
installation; the bioinformatics environment is not intended for native Windows.

The YAML describes the first environment solve; the explicit file records the exact
resolved Linux packages after installation. Keep that lock with the experiment.
An explicit Conda lock does not contain the editable pip package; keep the code revision too.
An existing exact environment can be reconstructed using `conda create --name nexvirome2
--file environments/linux-64.lock.txt` followed by the editable package installation.

## Reference source and selection

Export the **entire nucleotide accession list**, with versions, from the NCBI Virus view
linked in `configs/reference.yaml`. Save it to `data/refseq_accessions.txt`, one accession
per line. This is an explicit input: the code never substitutes a different database's
accession set. The downloader uses EFetch GBSeq XML in batches of 100, stores raw responses,
and requires returned accession.version sets to match exactly. No length, completeness,
or segment filters are applied at acquisition. GBSeq XML includes annotation and metadata;
unknown completeness stays unknown. The current metadata completeness flag is inferred
from the GenBank definition, not represented as a curated NCBI Virus completeness field.

Sources: [NCBI Virus download help](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/)
and [EFetch formats](https://www.ncbi.nlm.nih.gov/sites/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/).

Set dated snapshot, panel, design, run, and report paths in `configs/workflow.yaml`.
Use a new namespace for every changed input/configuration. Then:

```bash
snakemake -s workflows/reference.smk --cores 4 --dry-run
snakemake -s workflows/reference.smk --cores 4
nexvirome2 reference validate --snapshot data/references/refseq_viral/snapshot
snakemake -s workflows/Snakefile --cores 4 --dry-run
snakemake -s workflows/Snakefile --cores 4
```

The example validate path must match the snapshot path you chose. Run all commands from
the repository root. The initial workflow uses the activated Conda environment for every
rule. The workflow passes allocated threads to subprocesses. The genome FASTA files are
real dependencies, not just filenames embedded in the source specifications.

Panel selection keeps a disk-backed canonical 31-mer index. It ranks pairs by shared
non-low-complexity k-mers, verifies shared regions and distinguishable flanks with BLASTn,
and selects two risk pairs plus a disjoint low-sharing control. It reports insufficient
candidates rather than silently inventing a suitable pair. RNA/DNA/retrovirus records are
all eligible; a category is not guaranteed to appear among the top-scoring pairs.
Exact repeats can make the index and pair joins large; preserve `candidate_pairs.tsv`
for auditing. Panel `exclude_accessions` supports disjoint evaluation panels.

## Experiment design

The pilot is three pairs, three seeds, and mixture/single-A/single-B controls: 27 read
sets and 54 baseline assemblies. The expanded configuration uses coverage 10/50/100 and
ratios 1/10/100, producing 243 read sets for a three-pair panel. Coverage is defined for
the first source; the second is divided by the ratio. Singles retain the same coverage.

ART produces PE150 with HS25 error profiles, fragment mean 300 and SD 30. Its SAM output
provides private origin coordinates. FASTQ IDs are neutral and shuffled deterministically.
Circular references are padded to sample across the origin, padded-only fragments are
discarded, and coordinates are mapped back modulo the original length. Actual coverage
is recorded; it can differ from the requested coverage. Segments are never concatenated.

Change `experiment_config` to the expanded YAML and give design/runs/report new paths to
run the predefined expansion if no natural pilot errors appear. For independent evaluation,
exclude every development accession during panel construction and fill the same exclusion
list in `configs/experiments/heldout.yaml`. Freeze correction settings before running it.

## Interfaces and formats

`python -m nexvirome2 --help` lists all subcommands. Every producing command requires
`--config` (YAML or JSON) and `--output`; reference validation is read-only. `prepare`
expands a panel into explicit per-sample source JSON files and `experiments.json`.

Source JSON has `sources: [{accession, fasta, coverage, topology, seed_offset}]`. The seed
offset preserves identical source simulation seeds between mixtures and single-source controls.
Relative FASTA paths
are resolved against that JSON's directory. Generated design paths are absolute on the
execution server; generate the design on that server, not before moving the checkout.

Evaluation accepts raw sources only; correction accepts only contigs, GFA, SPAdes paths,
and the original paired reads. Precomputed evaluation alignments use this BLAST format:

```text
qseqid sseqid pident length qstart qend sstart send bitscore qseq sseq
```

BLAST input coordinates are 1-based inclusive. All reported breakpoint intervals and
fragment maps are 0-based. Fragment maps are half-open `[start,end)`. A breakpoint report
`start=end` denotes one exact cut; otherwise the possible cut lies in the closed range
`[start,end]`. Shared reference regions remain ambiguous if there is no unique evidence.

## Correction behavior and limitations

SPAdes numeric segment paths and sequence-containing GFA with exact `nM` overlaps are
supported. Missing/gapped paths and sequence mismatches remain unresolved. The prototype
examines a merge into a shared path followed by an exit branch. Competing local templates
include alternate incoming and outgoing anchors. Bowtie2 `-a` retains competing alignments;
equal-scoring paths or multiple placements do not count as unique evidence. Both ends must
reach the flanking anchors. Fragment endpoints are deduplicated.

Split requires at least three independent fragments for one alternate exit and at least
90% alternate support relative to the current exit. Absent spanning reads alone do not
trigger a split. Repeats exceeding the fragment span and oversized candidate sets remain
unresolved. Coverage is recorded as supporting context, not used as a standalone split rule.
Local mapping currently scans the original read library per candidate; this prioritizes
auditability in small panels and is not optimized for large real metagenomes.

All corrected fragments, including short pieces, are retained. A random control splits
the same original contigs the same number of times with a fixed seed. Full assembler
directories remain available, including logs, effective k values, and SPAdes intermediates.
See [SPAdes output](https://ablab.github.io/spades/output.html) and
[MEGAHIT usage](https://github.com/voutcn/megahit).

## Metrics and interpretation

Reports compare MEGAHIT, metaSPAdes, corrected metaSPAdes, and matched random splitting,
both without a length cutoff and with 500 bp minimum. They contain per-source compatible
alignment coverage, confirmed cross-source connections, ambiguous contigs, N50, longest
aligned blocks, and longest blocks from consistently assigned contigs. The latter summed
over source lengths is `correctly_reconstructed_fraction`; it is a conservative block
metric and not proof of complete genome reconstruction. Shared sequence coverage is
compatible with multiple sources and is not itself evidence of source-specific recovery.
Circular origin joins may split BLAST HSPs, making longest-block measures conservative.
Alignment gaps terminate continuous blocks; subject bases aligned to query gaps are excluded
from recovery. See [bug_audit.md](bug_audit.md) for the correction to earlier metric behavior.

Stage manifests record elapsed time. On Linux, `/usr/bin/time -v` records maximum resident
memory for external commands; missing resource measurements remain null, not zero.
Correction resource measurements cover subprocesses, not the Python parent's peak memory.

The target is 50% fewer confirmed wrong connections with <=0.02 absolute recovery loss.
When baseline errors are zero, reduction is undefined and the target is not claimed met.
Reports give per-sample targets, not a significance test or automatic generalization claim.

## Restart and verification

A stage creates a fresh output directory and refuses existing directories. Failures retain
logs and a failed manifest. To retry, archive the entire failed stage directory to a chosen
new location, then rerun. Snakemake skips intact completed outputs. Do not use forced reruns
against an existing output namespace. A manifest includes code checksums even when changes
have not been committed. Download/network errors retry with bounded backoff.

```bash
python -m pytest -q
snakemake -s workflows/reference.smk --configfile tests/fixtures/dryrun.yaml \
  --config design=.validation/design --cores 4 --dry-run
snakemake -s workflows/Snakefile --configfile tests/fixtures/dryrun.yaml --cores 4 --dry-run
```

The dry-run fixtures are for DAG validation only. They are deliberately too small for a
biological benchmark. Unit/integration tests mock ART and mapping subprocess output; this
checks wiring and invariants, and does not replace a real Linux ART/assembler/Bowtie2 run.

## Optional multi-sample and structural stages

Implemented extensions and their Conda installation are documented in
[extensions.md](extensions.md). They add feature extraction, NMF, coordinate-preserving ORFs,
sequence/structural searches with caching, structural candidate annotation, conservative
bin associations, geNomad/CheckV wrappers and five correction ablations. Candidate signals
remain separate from the evaluator's source truth. These stages are not enabled by the
base benchmark or required by its Conda environment.
