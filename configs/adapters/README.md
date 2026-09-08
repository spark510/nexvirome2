# Known adapter read-overlap screening panel

Prepared 2026-09-08 from [Illumina's trimming sequence guide](https://knowledge.illumina.com/library-preparation/general/library-preparation-general-reference_material-list/000001314).
`catalog.json` records sequences, applicable mates, orientation and source.
`screen_R1.fasta` and `screen_R2.fasta` contain four candidates each: TruSeq,
Nextera/DNA Prep, stranded RNA ligation and TruSeq Small RNA. Shared sequences
do not uniquely identify a library kit. These are read-facing trimming sequences,
not full oligos, index sequences or Trimmomatic palindrome prefixes.

Use this broad panel for diagnosis only. Select the actual kit before trimming.
Stranded RNA protocols may also need kit-specific 5-prime processing; it is not
encoded in this 3-prime panel. PCR-free tagmentation dual-sequence processing is
excluded pending explicit kit configuration. Reverse complements are not silently
added to 3-prime read-through detection.

## Inspect real paired FASTQ without trimming

### Download or refresh the panel

Run from the repository root using standard Python, with no additional packages:

```bash
python scripts/fetch_adapters.py --output data/adapters/illumina-20260908
```

This downloads the official Markdown document and extracts the five supported
entries with kit/mate validation. Outputs are `source.md`, `catalog.json`,
`screen_R1.fasta`, `screen_R2.fasta` and a checksum/timestamp `manifest.json`.
Download/parser failures do not fall back to hardcoded sequences.
Existing output directories are refused; use a fresh path for each refresh.
The committed panel is not overwritten. To use the download for screening,
replace the `configs/adapters/` paths below with your output directory.

For offline replay add `--source-file PATH/TO/source.md`. Local snapshots are
explicitly labeled unverified, with processing time separate from retrieval time.

### Run screening

Cutadapt is not yet in the project environment. On Linux, install it explicitly:

```bash
conda install --override-channels -c conda-forge -c bioconda cutadapt
```

Run from the checkout root, with real R1/R2 paths and a fresh output directory:

```bash
mkdir adapter_screen
cutadapt --action=none --no-indels -e 0.1 -O 8 \
  -a file:configs/adapters/screen_R1.fasta \
  -A file:configs/adapters/screen_R2.fasta \
  --json adapter_screen/report.json \
  -o adapter_screen/inspected_R1.fastq.gz \
  -p adapter_screen/inspected_R2.fastq.gz \
  sample_R1.fastq.gz sample_R2.fastq.gz > adapter_screen/report.txt
```

No quality/length filtering, sequence trimming or read renaming is requested.
Output FASTQs are diagnostic copies, not cleaned reads. See [Cutadapt's guide](https://cutadapt.readthedocs.io/en/stable/guide.html).

Repeat independently with minimum overlap 12 and 16 and fresh output paths;
run `-e 0 --no-indels` separately for exact matches. The above error-tolerant
screen is not an exact-match measurement. Regular 3-prime matching allows a full
adapter within a read and a partial adapter prefix at the read end. It does not
measure arbitrary internal partial matches or just read/read PE overlap.

Report R1/R2 hit fractions, candidate-specific counts and reported length/error
distributions at each threshold. Cutadapt chooses a best adapter per read here:
counts are not an exhaustive all-candidate overlap table. Nextera and stranded
RNA share sequence; resolve ambiguous kit attribution using kit metadata and
separate candidate runs if necessary. Do not sum R1/R2 counts to infer a unique
fragment contamination fraction. Do not interpret short partial hits alone as
confirmed contamination. A full per-read competing-hit/coordinate report is a
separate implementation task.

No real FASTQ has been analyzed as part of preparing this panel.
