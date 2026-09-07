# Multi-sample and structural extensions

These are executable research prototypes. They extend the existing graph/read-pair
correction; they do not construct a new de Bruijn graph or establish biological performance.
All examples run from the repository root on Linux x86_64.

## Installation

```bash
conda env create -f environments/extended.yaml
conda activate nexvirome2-extended
nexvirome2 doctor --extended
python -m pytest -q
conda list --explicit > environments/extended.linux-64.explicit.txt
```

The YAML installs the local CLI during Conda's pip phase. It includes the base assembly
tools, scikit-learn, Biopython, MMseqs2, Foldseek, Folddisco, FoldMason, geNomad and CheckV.
The explicit export should be made after successful installation on the actual Linux server;
a Windows cross-platform solver check is not proof that Linux executables start correctly.
Large databases and model weights are downloaded separately, not during Conda installation.
For Python-only development use `python -m pip install -e '.[extensions,test]'`.

## Common-target features and NMF

Create `samples.tsv` with `sample`, `biological_sample`, `r1`, `r2` columns, following
`configs/samples.example.tsv`. Read paths are relative to the TSV. Every library is mapped
to exactly the same contig FASTA, including competing contigs. There must be at least three
independent biological abundance profiles. Different seeds of one mixture are technical
replicates and must share a biological_sample value; their normalized coverage is averaged.
At least three biological samples must have nonzero coverage for factorization.

Alternatively supply `sam` and `target_sha256` columns. The digest must be the SHA256 of the
exact common target FASTA. SAM headers must also match its sequence names and lengths.
SAMtools name-sorts each import; the digest is a provenance assertion by the caller.

```bash
nexvirome2 features neighborhoods --contigs data/common/contigs.fasta \
  --graph data/common/assembly_graph_with_scaffolds.gfa --paths data/common/contigs.paths \
  --config configs/features.yaml --output runs/local
nexvirome2 features coverage --contigs data/common/contigs.fasta --samples configs/samples.tsv \
  --regions runs/local/regions.tsv --config configs/features.yaml --output runs/features
nexvirome2 features propose --coverage runs/features/normalized_coverage.tsv \
  --windows runs/features/windows.tsv --libraries runs/features/libraries.tsv \
  --composition runs/features/composition.json --config configs/features.yaml --output runs/simple
nexvirome2 latent fit --coverage runs/features/normalized_coverage.tsv \
  --windows runs/features/windows.tsv --libraries runs/features/libraries.tsv \
  --config configs/latent.yaml --output runs/nmf
```

Omit `--regions` for full-target windows. The local selector uses verified merge/repeat/exit
paths recognized by the current correction prototype; it is not an exhaustive bubble finder.
An empty selected neighborhood stops feature extraction explicitly. This prevents silently
expanding to the whole dataset. Mapping still scans each full library once; window selection
reduces feature storage and downstream model/search work, not the mapping input size.

Coverage counts uniquely mapped proper fragments, unions overlapping mates, rejects tied
best scores and low MAPQ, and normalizes depth by accepted aligned bases (per million).
This is a relative coverage feature, not an absolute molecule-count estimator. Composition
uses canonical 4-mers. The simple baseline reports coverage and composition cosine changes
separately. NMF uses nonnegative window-by-biological-sample coverage, row L1 normalization,
multiple initializations, rank/error diagnostics, convergence and aligned-component stability.
`rank: auto` uses an explicit error-plus-rank penalty heuristic, not a validated rank estimator.
Constant abundance profiles fail as unidentifiable. Factors are not virus identities.

`proposals.tsv` contains original-contig coordinates and a method-specific score, not a
probability. A window transition reports the union of the two adjacent windows because it
does not localize a nucleotide breakpoint. Gaps between selected regions are not bridged.

## ORFs, sequence and structural evidence

```bash
nexvirome2 orfs predict --contigs data/common/contigs.fasta --imported data/common/orfs.tsv \
  --config configs/structure.yaml --output runs/orfs
nexvirome2 structure database --name ProstT5 --config configs/structure.yaml --output data/databases/prostt5
nexvirome2 structure database --name PDB --config configs/structure.yaml --output data/databases/pdb
nexvirome2 search sequence --proteins runs/orfs/proteins.faa --database data/databases/proteins/db \
  --cache data/cache/sequence.sqlite --config configs/structure.yaml --output runs/sequence
nexvirome2 search structure --proteins runs/orfs/proteins.faa --database data/databases/pdb/db \
  --model data/databases/prostt5/db --cache data/cache/structure.sqlite \
  --config configs/structure.yaml --output runs/structure
nexvirome2 structure annotate --orfs runs/orfs/orfs.tsv --metadata configs/targets.tsv \
  --sequence-hits runs/sequence/hits.tsv --structure-hits runs/structure/hits.tsv \
  --config configs/structure.yaml --output runs/annotation
```

Supply an existing MMseqs2 protein database (or compatible protein FASTA) for sequence search.
PDB can be replaced with an explicitly selected database; `--name BFVD` is also supported.
No database is assumed to be equivalent to the original RefSeq nucleotide snapshot.
The target metadata TSV requires exact search `target` IDs and `origin` values
`viral`, `cellular`, or `uncertain`; optional fields are `family`, `lineage`, `hallmark`.
The example TSV is deliberately an uncertain placeholder. Curated metadata is required to
interpret hits; the program does not infer viral lineage from a protein family name.
Sequence and structure DB target IDs must use a consistent label namespace. Unlabelled or
conflicting competitive hits retain uncertainty. Sequence and structural bit scores are
not compared across tools; each search is summarized separately before combining labels.

Imported ORFs use zero-based, half-open `contig,start,end,strand` coordinates with in-frame
lengths. A terminal stop is removed from the peptide; coordinates retain the supplied interval.
`--regions` selects ORFs intersecting original-coordinate regions, without trimming proteins.
Without `--imported`, six-frame stop-to-stop translation generates candidate ORFs, including
partial ORFs. This is not a validated gene caller: overlapping/spurious candidates can obscure
structural change points. Prefer trusted viral gene calls for actual experiments. Translation
table and minimum length are configurable; no single genetic code covers every virus.

Identical peptides share a content-hash protein ID while every ORF retains its own coordinates.
The SQLite search cache includes sequence, tool version, search settings, database/model
artifact checksums, and coordinate-file identity when relevant. Hits and no-hit results are
cached. It reuses complete searches; it does not export a separately reusable per-protein 3Di
prediction cache. Hashing a large database costs I/O. Keep DB prefixes in dedicated folders.

Foldseek supports amino-acid queries through ProstT5 predicted 3Di, or coordinate input via
`--structure-map` in place of `--model`. The latter TSV has `protein_id,query_id,path`; query_id
must match the exact Foldseek chain identifier, paths are relative to that TSV, and each file
should contain one selected chain with a unique basename. The map must cover the selected
protein FASTA exactly. Predicted 3Di is not atomic coordinates: TM-score, LDDT and local
geometric motif claims are not produced by that route. See the [Foldseek documentation](https://github.com/steineggerlab/foldseek).

Contig annotation reports viral/cellular/uncertain evidence and sequence-close versus remote
structural candidates. These are discovery priorities, not proof of a new viral species.
Change points require consistent labels on multiple ORFs on both sides; functional-family
changes within one lineage do not trigger lineage proposals. Real mosaic genomes still need
read/graph validation. Search E-values and coverage cutoffs remain experimental settings.

## Coordinate-based motif and comparison tools

```bash
nexvirome2 structure motif --structures data/structures --query data/query.pdb \
  --residues B57,B102,C195 --config configs/structure.yaml --output runs/motifs
nexvirome2 structure align --structures data/selected_structures \
  --config configs/structure.yaml --output runs/structure_alignment
```

Folddisco requires actual structures and an explicitly chosen chain/residue motif; use
`--index` to reuse an existing index. Its raw motif hits are retained for review and do not
automatically become functional labels or chimera calls. FoldMason aligns at least two
coordinate files and retains its native report. These wrappers follow the official
[Folddisco](https://github.com/steineggerlab/folddisco) and
[FoldMason](https://github.com/steineggerlab/foldmason) command interfaces.

## Correction ablations and bins

Copy/edit `configs/extensions_workflow.yaml` with actual contigs, graph, paths, mapping
libraries, correction read pair, database/model paths and curated labels. Then run:

```bash
snakemake -s workflows/extensions.smk --configfile configs/extensions_workflow.yaml --cores 16 -n
snakemake -s workflows/extensions.smk --configfile configs/extensions_workflow.yaml --cores 16
```

`local_only: true` restricts coverage windows and ORF searches to graph neighborhoods.
`orf_coordinates` optionally imports trusted calls. The workflow runs baseline, simple,
NMF, structural and combined correction variants. Every variant requires the original
discriminating read-pair support before splitting. The simple/NMF/structural variants also
require a matching candidate interval. Combined requires both `nmf` and `structural_lineage`;
an origin-only structural proposal does not satisfy this particular ablation.

For individual commands, `correct --proposals file.tsv ...` defaults to annotation only.
Set `proposal_policy: require_any` or `require_all` (with explicit `required_methods`) in
the correction config to gate cuts. Candidate evidence can downgrade a split to unresolved;
it cannot promote an unresolved/retained read-supported connection to a split. The audit
lists unmatched proposals, which remain unresolved. Truth is never passed to correction.

The workflow's `summary/ablations.tsv` compares split and unresolved counts only. Evaluate
each corrected FASTA and corresponding random control using the existing `evaluate` command
and held-out source truth to establish correctness/recovery. Fewer cuts alone are not success.

`bins build` takes the same coverage/windows/libraries/composition files and optional
`--profiles runs/annotation/orf_profiles.tsv`. It conservatively associates contigs by
coverage and composition, with optional structural-family Jaccard evidence. Internally
inconsistent contigs remain unresolved. Connected components are candidate associations,
not ordered genomes, species calls or sequence joins; the module emits no concatenated FASTA.
Default pairwise binning is limited to 2,000 selected contigs. Structural weighting is a
heuristic and should be compared with `structural_weight: 0` on a frozen evaluation panel.

## Independent classification and quality checks

```bash
nexvirome2 qc database --tool genomad --config configs/structure.yaml --output data/databases/genomad
nexvirome2 qc database --tool checkv --config configs/structure.yaml --output data/databases/checkv
nexvirome2 qc run --tool genomad --contigs data/common/contigs.fasta \
  --database data/databases/genomad/database/genomad_db --config configs/structure.yaml --output runs/genomad
nexvirome2 qc run --tool checkv --contigs data/common/contigs.fasta \
  --database PATH_TO_EXTRACTED_CHECKV_DB --config configs/structure.yaml --output runs/checkv
```

Use the actual extracted database directory reported by the download. Native outputs,
commands, versions and checksums are preserved. These checks are separate from the optional
workflow and do not alter contigs or supply chimera truth. Completeness/classification scope
depends on the tool and database; it does not cover every RNA/DNA/retroviral genome equally.
Interfaces follow [geNomad](https://portal.nersc.gov/genomad/quickstart.html) and
[CheckV](https://pypi.org/project/checkv/).

## Remaining research scope

The current extension splits supported bad connections; it does not reconstruct lost graph
paths, join contigs, train sequence embeddings, resolve every strain, run Phold as a separate
backend, or perform long-read/COBRA extension. ProstT5 through Foldseek is the implemented
fast structural representation route. Geometric motif outputs need curated interpretation.
Natural RefSeq chimera reproduction and held-out improvement targets remain unvalidated.

Each stage refuses an existing output directory and records failure explicitly. Successful
Snakemake outputs are reused; a failed stage must be archived under a new name, or the output
root changed, before retrying. Do not use `--rerun-incomplete` expecting silent overwrites.
