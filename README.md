# NexVirome2

The [integrated implementation plan](docs/implementation_plan.md) tracks reference masking
refinement, assembly correction, conserved-region reuse, and all audited remaining work.
It distinguishes implemented prototypes from missing functionality and conditional research.

The [masking and evidence-driven execution guide](docs/masking_pipeline.md) covers the
new taxonomy/sharing atlas, reversible reference masks, same-query comparisons,
conserved-region evidence, grouped holdout splits, and frozen production decisions.
`workflows/masking.smk` includes a runnable synthetic fixture. Real legacy NexVirome
equivalence and natural-panel performance remain unvalidated.

Further implemented modules—condition-specific policy selection, mixed-feature NMF,
read QC, conserved-support reuse, near-identity catalog and targeted coassembly—are
documented in [remaining extensions](docs/remaining_extensions.md).

NexVirome2 is a reproducible benchmark and prototype for ambiguity-aware short-read viral
assembly. The initial study uses a dated NCBI RefSeq Viral snapshot, controlled mixed-virus
read simulation, existing assemblers, and source-aware chimera evaluation.

The repository is intentionally split into stable code/configuration and generated data. Large
references, reads, assemblies, and run outputs belong under `data/` and `runs/` and are ignored
by Git. See [docs/study_design.md](docs/study_design.md) for the experimental sequence and
[docs/chimera_definition.md](docs/chimera_definition.md) for the evaluation rule.

The first implementation milestone is:

```text
RefSeq snapshot → virus-pair panel → simulated PE150 reads → MEGAHIT/metaSPAdes
→ source-aware chimera calls → graph inspection → ambiguity-aware correction prototype
```

Implementation, Linux installation, commands, input formats, and limitations are documented
in [docs/implementation.md](docs/implementation.md).

The CLI supports `reference fetch/validate`, `panel build`, `prepare`, `simulate`, `assemble`,
`evaluate`, `correct`, and `report`. The first workflow freezes a RefSeq accession list and
creates a panel/design; the second runs simulation, both assemblers, correction, and evaluation.

```bash
python scripts/install.py
conda activate nexvirome2-extended
nexvirome2 doctor --extended
python -m pytest -q
```

Run these commands from the repository root on **Linux x86_64**. Conda installs the external
tools and uses its pip phase to install the local NexVirome2 CLI automatically. No separate
pip command is needed. This is installation from a checkout, not a published Conda-channel
package. `doctor` checks PATH and tool startup; it does not run an assembly benchmark.

The installer defaults to the extended environment, which includes dependencies for the
full test suite and implemented optional modules. Use `python scripts/install.py --base`
for the core assembly tools only, or add `--dry-run` to check dependency resolution without
installing. Existing environments are not overwritten.

Only `conda-forge` and `bioconda` are allowed, in that order with strict channel priority.
All environment YAML files use `nodefaults` (an exclusion marker, not a download channel).
The installer also sets a channel allowlist during creation and persists strict priority
in the new environment. It does not change global Conda settings. For later package additions,
use `conda install --override-channels -c conda-forge -c bioconda PACKAGE`: without this flag,
Conda can merge channels from your global configuration. Recreate environments through the
installer to retain its installation restrictions.
The local Python package is installed through pip; Conda channel restrictions apply to Conda
packages, not pip build dependencies.

Installation prepares software, not experiment inputs. You must still supply the NCBI Virus
accession list, fetch and validate the RefSeq snapshot, prepare taxonomy/panels and run paths,
and supply reads or simulation settings. Structural/quality workflows also need their
databases and any required model weights. Natural-panel benchmarks and a full external-tool
Linux run remain necessary before treating this as a validated production pipeline.

Bioinformatics commands require this Linux environment. Set the
snapshot/run paths in `configs/workflow.yaml` and supply the versioned accession list exported
from NCBI Virus. Tests use synthetic inputs; biological performance has not yet been validated.

Optional multi-sample NMF, structural evidence, candidate binning, and correction ablations
are implemented as research prototypes. See the [extension runbook](docs/extensions.md)
for input tables, database preparation, interpretation limits, and the third workflow.

```bash
python scripts/install.py
conda activate nexvirome2-extended
nexvirome2 doctor --extended
```

The extension CLI adds `features neighborhoods/coverage/propose`, `latent fit`, `orfs predict`,
`search sequence/structure`, `structure database/annotate/motif/align`, `bins build`, and
`qc database/run`. NMF and structural proposals require read/graph validation before a cut.
The [NMF note](docs/latent_resolution.md) and [structural note](docs/structural_evidence.md)
retain the research rationale; numerical tests do not establish biological effectiveness.

See [architecture.md](docs/architecture.md) for module boundaries, injectable execution/search
backends, typed evidence, correction policies, and compatibility guarantees.
