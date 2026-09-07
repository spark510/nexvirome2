# Research note: NMF and local graph resolution

Implementation gaps and subsequent experiments are tracked under R10–R14/R18 in the
[integrated plan](implementation_plan.md). Coverage-only NMF is implemented; multi-feature
path resolution and a biologically varied multi-sample benchmark remain separate tasks.

Source: user-provided conversation about NMF, contig assembly, and computational cost,
received 2026-09-06. This records candidate research directions, not validated results.

The accepted direction is to retain existing graph construction and investigate additional
evidence at ambiguous junctions. The current implementation uses graph paths and paired-end
support. The optional extension now adds cross-sample coverage and composition and compares
NMF with those simpler measurements. See [extensions.md](extensions.md) for executable steps.

## Proposed sequence of experiments

1. Freeze the current read-pair baseline and independent evaluation panel.
2. Map several genuinely different samples to one common contig/graph coordinate system;
   quantify coverage and mapping ambiguity for the candidate neighborhood only.
3. Test coverage-pattern similarity as an additional feature, with composition as a separate
   ablation. Measure error rate, recovery, runtime, peak memory, and unresolved junctions.
4. Only then compare a nonnegative factorization of the neighborhood-by-sample coverage
   matrix. Record normalization, rank, initialization seeds, convergence, and stability.
5. Use a change in factor loading as candidate evidence, subject to the same read/graph
   verification. Do not turn a factor change directly into an unconditional contig split.

Random simulation seeds at the same abundances are technical replicates, not evidence of
biological abundance covariance. Multi-sample experiments need independently varied source
abundances and a fixed mapping target. Shared nodes can combine coverage from several sources,
and correlated sources may be inseparable by coverage alone.

NMF components are latent numerical factors: they do not automatically identify a virus,
host, plasmid, or contamination. Those interpretations need independent evidence. Likewise,
within-genome composition changes or real recombination can mimic a factor transition.
Fixed 2–5 kb windows from the pasted discussion should not become an unexamined default for
short viral genomes; candidate windows must be evaluated for length and informative content.

## Computational boundary

Avoid a reads-by-all-kmers matrix. Begin with local graph nodes/windows by samples, retain
mapping uncertainty, and benchmark feature extraction as well as factorization. The current
per-candidate mapping scans the whole library; an efficient reusable mapping/index should
be evaluated before adding an expensive model. Threshold/rank tuning must not use held-out
source labels or the simulation truth available only to evaluation.

NMF is implemented in the optional extended Conda environment. It remains disabled in the
base benchmark. Convergence, rank diagnostics, seed stability, technical-replicate grouping
and read/graph gating are implemented; actual multi-sample biological performance is untested.

See [structural_evidence.md](structural_evidence.md) for the complementary proposal to use
ORF-level structural annotations as candidate breakpoint evidence, with the same validation
and uncertainty requirements.
