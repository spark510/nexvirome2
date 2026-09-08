#!/usr/bin/env python3
"""Download the diagnostic adapter panel from the official Illumina guide."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

SOURCE = ('https://knowledge.illumina.com/library-preparation/general/'
          'library-preparation-general-reference_material-list/000001314.md')


def parse_adapters(document):
    """Associate each sequence with its kit paragraph and explicit mate label."""
    paragraphs = re.split(r'\n\s*\n', document.replace('\r\n', '\n'))
    definitions = [
        ('truseq', 'TruSeq single index', [('truseq_r1', 'R1'), ('truseq_r2', 'R2')]),
        ('nextera_dna_prep', 'AmpliSeq for Illumina;', [('nextera_dna_prep', 'both')]),
        ('stranded_rna_ligation', 'Illumina Stranded mRNA Prep, Ligation', [('stranded_rna_ligation', 'both')]),
        ('truseq_small_rna', 'TruSeq Small RNA:', [('truseq_small_rna', 'both')]),
    ]
    entries = []
    for family, heading, mates in definitions:
        positions = [i for i, p in enumerate(paragraphs) if p.strip().startswith(heading)]
        if len(positions) != 1 or positions[0]+1 >= len(paragraphs):
            raise ValueError(f'Official document layout changed: {family}')
        block = paragraphs[positions[0]+1].strip().splitlines()
        if len(block) != len(mates):
            raise ValueError(f'Unexpected adapter count: {family}')
        for line, (identifier, mate) in zip(block, mates):
            label = f'Read {mate[-1]}: ' if mate != 'both' else ''
            match = re.fullmatch(r'\* ' + re.escape(label) + r'([ACGT]{12,100})', line.strip())
            if match is None:
                raise ValueError(f'Unexpected sequence or mate label: {identifier}')
            entries.append(dict(id=identifier, family=family, mate=mate,
                                sequence=match[1], orientation='read_5to3',
                                expected_location='3prime_readthrough', source_url=SOURCE))
    return entries


def fetch(output, source_file=None):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f'Output already exists; choose a new --output: {output}')
    if source_file is None:
        with urlopen(Request(SOURCE, headers={'User-Agent':'NexVirome2/0.1 adapter-fetch'}), timeout=45) as response:
            raw = response.read(2_000_001)
            resolved_url = response.geturl()
        mode = 'official_download'
    else:
        with open(source_file, 'rb') as handle:
            raw = handle.read(2_000_001)
        resolved_url = None
        mode = 'local_snapshot_unverified'
    if len(raw) > 2_000_000:
        raise ValueError('Unexpected document size')
    entries = parse_adapters(raw.decode('utf-8-sig'))
    now = datetime.now(timezone.utc).isoformat()
    for entry in entries:
        entry['retrieved_date'] = now[:10] if source_file is None else None
    catalog = dict(schema=1, purpose='diagnostic_read_overlap_not_automatic_trimming',
                   entries=entries, excluded=[
                       'PCR-Free tagmentation dual-sequence policy requires kit-specific configuration',
                       'index-specific full oligos', 'PCR primers and Trimmomatic palindrome prefixes'])
    files = {'source.md':raw, 'catalog.json':(json.dumps(catalog, indent=2)+'\n').encode('utf-8')}
    for mate in ('R1','R2'):
        files[f'screen_{mate}.fasta'] = ''.join(
            f">{r['id']}\n{r['sequence']}\n" for r in entries
            if r['mate'] in (mate,'both')).encode('ascii')
    # Exclusive creation prevents simultaneous runs from overwriting one another.
    # A failed write cannot leave a completed manifest.
    output.mkdir(parents=True, exist_ok=False)
    for name, data in files.items():
        (output/name).write_bytes(data)
    manifest = dict(schema=1, status='complete', prepared_at_utc=now,
                    source_url=SOURCE, resolved_url=resolved_url, source_mode=mode,
                    source_file=str(Path(source_file).resolve()) if source_file else None,
                    files={name:hashlib.sha256(data).hexdigest() for name,data in files.items()},
                    real_reads_analyzed=False)
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='New output directory; never overwritten.')
    parser.add_argument('--source-file', type=Path, help='Replay a saved Markdown document offline (unverified snapshot).')
    args = parser.parse_args()
    fetch(args.output, args.source_file)
    print(f'Prepared 5 adapter entries and R1/R2 FASTAs in {args.output}')


if __name__ == '__main__':
    main()
