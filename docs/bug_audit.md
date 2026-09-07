# Bug audit — 2026-09-06

Five defects were reproduced with failing tests before correction. The audit focused on
data integrity, extension evidence handling, and evaluation metrics; it is not an exhaustive
review or a Linux external-tool benchmark.

| Defect | Trigger and effect | Fix |
|---|---|---|
| Recovery includes missing sequence | A 300-base contig aligns across a one-base source deletion; the old bounding interval counts 301 recovered bases. Gaps also inflate continuous-block length. | Accumulate ungapped blocks continuous in both query and subject coordinates. Test insertions/deletions in both subject orientations. |
| Duplicate TSV headers overwrite data | Two sample columns have the same name; DictReader silently retains only the last value. | Reject duplicate/empty column names and malformed column counts; accept UTF-8 BOM input. |
| Input failure loses provenance | A stage creates its output directory, then fails to hash a missing input before creating its manifest. | Initialize the manifest first and include setup/checksumming in the failure handler. |
| Missing structure becomes negative evidence | Two contigs have identical coverage/composition, but only one has a structural family; a zero Jaccard score prevents association. | Structural similarity is unavailable unless both contigs have annotated families; absent evidence has zero weight. |
| Imported ORF hides internal stop | Translation ends in two stop codons; stripping all trailing stops incorrectly accepts the first, internal stop. | Remove at most one terminal stop and reject any remaining stop, preserving valid terminal-stop coordinates. |

The full suite passes **57 tests**, including the new reproductions and boundary cases in
`tests/test_bug_regressions.py`. Existing correction, NMF, feature extraction, search-cache,
reference, and workflow-interface tests remain passing. No external Linux tool execution
was performed during this audit.

Previously generated evaluation metrics are immutable artifacts and are not rewritten.
Re-run `evaluate` into a fresh output directory before comparing gap-sensitive recovery or
continuous-block results. Existing bin outputs with missing structural annotations likewise
need a fresh run to reflect the corrected evidence handling.
