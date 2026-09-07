# Study design

The complete forward implementation scope and requirement audit are maintained in
[implementation_plan.md](implementation_plan.md). This study includes reference masking
refinement, assembly-ambiguity evaluation/correction, and conserved-region reuse for discovery.
Classification ambiguity and assembly chimera are measured separately and then related;
one is not assumed to imply the other. The sequence below remains the initial assembly pilot.

The reference source is the NCBI Virus VSSI nucleotide view filtered to `SourceDB=RefSeq`.
The browser URL and exported accession.version list define provenance. The implementation
retrieves GenBank XML records using EFetch and reconciles the returned set against that list.
Other download sources must not silently replace the frozen accession set. See
[implementation.md](implementation.md) for the executable workflow and current limitations.

The first milestone is a controlled benchmark with NCBI RefSeq Viral genomes:

1. Freeze a dated reference snapshot and record accession versions and checksums.
2. Select related virus pairs and simulate paired-end reads with known read origins.
3. Run existing assemblers (initially MEGAHIT and metaSPAdes).
4. Align contigs back to the source genomes and quantify true chimeras.
5. Inspect the de Bruijn graph around validated breakpoints.
6. Prototype ambiguity-aware contig breaking or path scoring.
7. Compare chimera rate, genome recovery, completeness, and correctly reconstructed fraction.

The simulation truth must remain separate from the inputs used by the correction method.
The user's NMF/multi-sample graph-resolution discussion is retained as a follow-up research
direction in [latent_resolution.md](latent_resolution.md); it is not a mandatory assembly step.
The structural discovery/binning/QC discussion is recorded in
[structural_evidence.md](structural_evidence.md). Its direct assembly extension is to propose
breakpoint candidates from ORF evidence and validate them with the existing read/graph framework.
The first claim should be whether shared/conserved sequence causes assembly chimeras under
specified similarity, abundance, coverage, insert-size, and k-mer conditions.
