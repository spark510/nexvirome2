# Bounded graph reconstruction prototype

`reconstruct run` spells a supported alternative path into the output contig,
including length changes. It preserves original inputs and uses existing GFA
segments only; it does not fill gaps or reconstruct paths removed by assembly.

```bash
nexvirome2 reconstruct run --config configs/reconstruction.yaml \
  --contigs contigs.fasta --graph assembly_graph_with_scaffolds.gfa \
  --paths contigs.paths --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz \
  --output runs/reconstruction-inspect
```

Default `apply: false` generates proposals and leaves `reconstructed.fasta`
identical to the input. Set `apply: true` in a separate config and use a fresh
output directory to emit replacements. `proposed.fasta` contains eligible
alternative sequences; `decisions.json` records support, reasons and paths.
`graph_coordinate_map.tsv` uses zero-based half-open output and forward-node
coordinates with explicit strand. New interior sequence has no original-contig
coordinate. Original coordinates of retained suffixes account for length changes.

The first implementation supports simple paths with two common boundary nodes.
This includes a direct original edge replaced by a supported multi-node detour.
Cycles, gapped/mismatched original paths, incomplete enumeration and regions too
long for configured fragments are unresolved. Multiple supported changes on one
contig are all withheld pending joint validation. Paths ending at an unrelated
exit cannot be substituted for the original suffix.

Bowtie2 maps all input pairs competitively against full candidate contigs and
other original contigs. A counted pair must span both outside anchors with at
least 20 aligned bases on each side, connecting the changed region. Default
support is at least 3 endpoint-deduplicated fragments and 90% of discriminating
support, with a positive alignment score margin. All enumerated alternatives
compete; ties and QC/duplicate-marked fragments cannot support reconstruction.
Identical candidate sequences therefore remain ambiguous, rather than creating
evidence from duplicate graph paths.

Alignment competition uses both mates' positions, CIGARs and orientations;
endpoint-based fragment deduplication is applied only after uniqueness is established.
Equal-scoring different placements cannot become unique support merely because
their outer endpoints coincide. Fragment audits retain the mate placements.

Every non-null entry in `output_paths.json` must spell its output FASTA exactly.
Invalid original paths are exported as null with reasons in
`output_path_validation.json`; the original sequence is retained. A mismatch
in an applied reconstruction fails the run instead of exporting an invalid path.

Limitations: this is PE-only reconstruction, not integrated NMF/biological path
scoring or strain phasing. Other original contigs are competitors, but absent
sources and graph paths outside the enumerated regions need not be represented.
Mapping is heuristic, not a proof of global uniqueness. Each candidate still
requires a whole-read mapping pass. Applied candidates are remapped together;
if any replacement loses minimum spanning support, the entire provisional batch
is withheld to avoid unverified interactions after rollback. This is an input-read
consistency check, not independent validation. Biological heldout validation
remains required before production adoption. No natural-data
improvement or Linux external-tool execution is claimed by synthetic tests.
