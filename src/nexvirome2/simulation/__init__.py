"""ART wrappers with private coordinate truth and neutral paired read IDs."""
import gzip
import random
from pathlib import Path

from ..common import command, dump, fasta, stage, table, version, write_fasta


def fastq(path):
    with open(path) as handle:
        while True:
            name = handle.readline().rstrip()
            if not name:
                break
            seq, plus, quality = (handle.readline().rstrip() for _ in range(3))
            if not name.startswith('@') or not plus.startswith('+') or len(seq) != len(quality):
                raise ValueError(f'Malformed FASTQ: {path}')
            yield name[1:].split()[0], seq, quality


def sam_truth(path):
    result = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            flag = int(f[1])
            if flag & 4:
                continue
            mate = 1 if flag & 64 else 2
            key = f[0].removesuffix('/1').removesuffix('/2')
            result[(key, mate)] = (int(f[3])-1, '-' if flag & 16 else '+')
    return result


def simulate(settings, source_file, output):
    """sources JSON: list of {accession, fasta, coverage, topology}."""
    import json
    source_file = Path(source_file)
    spec = json.loads(source_file.read_text())
    sources = spec['sources']
    paths = [(source_file.parent / s['fasta']).resolve() for s in sources]
    with stage(output, settings, [source_file, *paths]) as (out, manifest):
        truth = out / 'truth'
        truth.mkdir()
        working = truth / 'art'
        working.mkdir()
        version('art_illumina', '-h', out, manifest)
        seed = int(settings.get('seed', 20260905))
        rng = random.Random(seed)
        read_length = int(settings.get('read_length', 150))
        records, origin, all_sources = [], [], {}
        for index, (source, path) in enumerate(zip(sources, paths)):
            sequence = fasta(path)[source['accession']]
            if source['accession'] in all_sources:
                raise ValueError('Duplicate simulation source')
            all_sources[source['accession']] = sequence
            circular = source.get('topology', 'linear') == 'circular'
            padding = int(settings.get('insert_mean', 300) + 8 * settings.get('insert_sd', 30)) if circular else 0
            if padding >= len(sequence):
                raise ValueError('Source too short for configured circular fragment range')
            simulated = sequence + sequence[:padding]
            input_fa = working / f'source_{index}.fasta'
            write_fasta(input_fa, {'template': simulated})
            prefix = working / f'source_{index}_'
            command(['art_illumina', '-ss', settings.get('profile', 'HS25'), '-p', '-sam', '-na',
                     '-i', input_fa, '-l', read_length, '-f', source.get('coverage', 50),
                     '-m', settings.get('insert_mean', 300), '-s', settings.get('insert_sd', 30),
                     '-rs', seed + source.get('seed_offset', index), '-o', prefix], out, manifest, f'art_{index}')
            coordinates = sam_truth(str(prefix)+'.sam')
            left, right = list(fastq(str(prefix)+'1.fq')), list(fastq(str(prefix)+'2.fq'))
            if len(left) != len(right):
                raise ValueError('ART mate count mismatch')
            for r1, r2 in zip(left, right):
                key1 = r1[0].removesuffix('/1').removesuffix('/2')
                key2 = r2[0].removesuffix('/1').removesuffix('/2')
                if key1 != key2:
                    raise ValueError('ART mate names differ')
                c1, c2 = coordinates[(key1, 1)], coordinates[(key2, 2)]
                if circular and min(c1[0], c2[0]) >= len(sequence):
                    continue
                records.append((r1, r2, source['accession'], c1, c2, len(sequence)))
        rng.shuffle(records)
        with gzip.open(out / 'reads_R1.fastq.gz', 'wt') as r1out, gzip.open(out / 'reads_R2.fastq.gz', 'wt') as r2out:
            for i, (r1, r2, accession, c1, c2, length) in enumerate(records):
                neutral = f'fragment_{i:012d}'
                for mate, read, coordinate, handle in ((1, r1, c1, r1out), (2, r2, c2, r2out)):
                    handle.write(f'@{neutral}/{mate}\n{read[1]}\n+\n{read[2]}\n')
                    origin.append({'read_id': f'{neutral}/{mate}', 'accession': accession,
                                   'start': coordinate[0] % length, 'strand': coordinate[1],
                                   'wraps_origin': coordinate[0] % length + len(read[1]) > length})
        write_fasta(truth / 'sources.fasta', all_sources)
        table(truth / 'read_origin.tsv', origin, ['read_id', 'accession', 'start', 'strand', 'wraps_origin'])
        dump(truth / 'sources.json', spec)
        manifest['read_pairs'] = len(records)
        manifest['realized_coverage'] = {s: sum(2*read_length for r in records if r[2] == s)/len(seq) for s, seq in all_sources.items()}
