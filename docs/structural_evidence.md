# Research note: structural evidence for viral contigs and assembly QC

Remaining discovery, motif integration, reference-mask reuse, and real-data requirements are
tracked under R07/R16–R21 in the [integrated plan](implementation_plan.md). Tool wrappers
alone do not constitute those complete analyses.

Source: user-provided conversation on structural discovery, viral binning, chimera
detection, and assembly strategies, received 2026-09-06. This note preserves proposals;
it does not establish their effectiveness. The implemented optional stages are described in
[extensions.md](extensions.md), including their remaining experimental limitations.

## Three distinct research outputs

1. Discovery: prioritize sequence-divergent candidates using sequence, structural,
   functional-motif, and genomic-context evidence. No hit is not a viral identity call.
2. Classification/binning: investigate whether structural annotations improve viral,
   cellular, and uncertain assignments and association of contigs. A bin is an association
   of sequences, not a reconstructed ordered genome; shared folds alone must not join them.
3. Assembly QC: use changes in ORF-level evidence to propose breakpoint intervals, then
   validate them against reads and assembly graph paths. This is the direct extension of
   the current NexVirome2 objective.

The current implementation retains metaSPAdes/MEGAHIT plus read/graph-based correction.
The optional structural and [NMF](latent_resolution.md) modules feed candidate/evidence
tables and support-gated correction ablations, preserving the baseline and truth boundary.

## Proposed structural-QC experiment

Preserve original contig IDs, ORF nucleotide coordinates, strand, translation settings,
and links to graph paths. Start with ambiguous graph neighborhoods and uncertain contigs.
Compare the frozen read/graph baseline with sequence/composition features, structural
features, and their combination on the same independent evaluation set.

The pasted discussion proposes MMseqs2/DIAMOND/HMMER for sequence evidence, Foldseek for
structural search, Folddisco for local-motif evidence, and FoldMason for follow-up comparison
of candidate structures. ProstT5/Phold are candidate routes for reducing prediction cost.
For every real experiment, verify each tool's supported input representation, organism scope,
database provenance, confidence measures, and computing requirements. A representation
used for a global search must not be assumed sufficient for coordinate-based motif tests.

Cache evidence by protein sequence and model/database version, retaining uncertainty and
mapping it back to each contig. Reuse the evidence for QC and later classification instead
of repeatedly predicting the same protein. Investigate ORF/window change points as proposed
breakpoint intervals; do not pretend protein-level resolution gives an exact nucleotide cut.

## Interpretation and evaluation safeguards

- Structural family, functional similarity, taxonomic lineage, and genome membership are
  different labels; do not substitute one for another without independent evidence.
- Include genuine recombinant/mosaic source genomes as negative controls for artificial
  assembly joins. A structural transition or a viral/cellular boundary is not sufficient
  to establish an assembly artifact. Integrated viral sequence is also a competing explanation.
- Keep host-flanking annotation, host-virus assembly errors, and virus-virus assembly errors
  as separate outcomes. Proposed geNomad/CheckV comparisons need task-matched endpoints;
  neither an annotation label nor agreement between assemblers defines ground truth.
- Assess raw-read spanning support, alternative graph paths, and mapping ambiguity. Lack
  of spanning pairs by itself remains insufficient for a split, as in the existing resolver.
- Evaluate cross-sample abundance only after defining a shared mapping target and accounting
  for repeats, ambiguous mapping, source abundance, and library effects. Correlation is
  supporting evidence, not proof of one genome or proof that a mosaic is biological.
- Novelty assessment must record reference and structural database versions and potential
  leakage from held-out sources. Retain an uncertain category rather than automatically
  deleting sequences with no confident annotation.

For real data, host subtraction, per-sample assembly, targeted coassembly, and assembly
ensemble/dereplication are experimental choices needing the actual sample and host metadata.
The human-specific examples in the pasted text do not establish that this project's samples
are human. No blanket host subtraction is added to the current synthetic benchmark.

## Deferred decisions

Long-read/hybrid assemblers and COBRA are recorded as candidates for later comparison if
the available data and the reconstruction objective warrant it. The pasted database sizes,
performance multipliers, model scores, distances, and claims about the latest methods are
unverified examples; they are not adopted as thresholds, benchmark facts, or tool rankings.

The pasted preference for MEGAHIT as a primary assembler does not supersede the agreed
metaSPAdes graph-development baseline. Both remain independent assembly comparison arms.
Structural tools, models, and databases are not added to the default Conda installation.
The implemented extension has a separate optional environment and explicit model/database
downloads; additional backends must preserve that separation and installation procedure.
