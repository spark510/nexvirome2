"""Explicitly selected paired-read pools with sample lineage."""
from pathlib import Path
from collections import defaultdict, Counter
from contextlib import ExitStack
from itertools import zip_longest
import csv
import gzip
import re
from .common import rows, stage, table, dump
from .read_preparation import records, fragment_id


def prepare(settings, samples_file, membership_file, output):
    """Samples: sample,biological_sample,library,r1,r2. Membership: group,sample,fragment."""
    samples = {}
    for row in rows(samples_file):
        name = row['sample']
        if name in samples or not re.fullmatch('[A-Za-z0-9_-]+', name):
            raise ValueError('Duplicate/unsafe sample ID')
        if not row['biological_sample'] or not row['library']:
            raise ValueError('Samples require biological_sample and library')
        samples[name] = {**row, **{key: (Path(samples_file).resolve().parent/row[key]).resolve() for key in ('r1','r2')}}
    selected, group_samples, seen = defaultdict(lambda: defaultdict(set)), defaultdict(set), set()
    for row in rows(membership_file):
        group, sample, fragment = row['group'], row['sample'], row['fragment']
        if not re.fullmatch('[A-Za-z0-9_-]+', group) or sample not in samples or not fragment:
            raise ValueError('Invalid coassembly membership')
        key = (group, sample, fragment)
        if key in seen:
            raise ValueError('Duplicate coassembly membership')
        seen.add(key)
        selected[sample][fragment].add(group)
        group_samples[group].add(sample)
    if not group_samples or len(group_samples) > int(settings.get('max_groups', 100)):
        raise ValueError('Empty or excessive coassembly groups')
    for names in group_samples.values():
        if len({samples[name]['library'] for name in names}) != 1:
            raise ValueError('Coassembly group mixes incompatible library types')
    with stage(output, settings, [samples_file, membership_file, *[s[k] for s in samples.values() for k in ('r1','r2')]]) as (out, _):
        counts = Counter()
        with ExitStack() as stack:
            pools = {}
            for group in sorted(group_samples):
                (out/group).mkdir()
                pools[group] = [stack.enter_context(gzip.open(out/group/f'reads_R{i}.fastq.gz', 'wt', newline='\n')) for i in (1,2)]
            lineage = csv.DictWriter(stack.enter_context(open(out/'lineage.tsv','w',encoding='utf-8',newline='')),
                fieldnames=['group','pooled_fragment','sample','biological_sample','original_fragment'], delimiter='\t')
            lineage.writeheader()
            for sample, source in sorted(samples.items()):
                visited = set()
                for first, second in zip_longest(records(source['r1']), records(source['r2'])):
                    if first is None or second is None or fragment_id(first[0]) != fragment_id(second[0]):
                        raise ValueError('Coassembly reads have missing/out-of-order mates')
                    if first[0].endswith('/2') or second[0].endswith('/1'):
                        raise ValueError('Coassembly mate suffix conflicts with R1/R2 input')
                    original = fragment_id(first[0])
                    if original in visited:
                        raise ValueError('Duplicate coassembly fragment ID')
                    visited.add(original)
                    pooled = f'{sample}:{original}'
                    for group in sorted(selected[sample].get(original, ())):
                        for mate, (handle, read) in enumerate(zip(pools[group], (first[1], second[1])), 1):
                            handle.write(f'@{pooled}/{mate}\n{read[1]}\n+\n{read[3]}\n')
                        lineage.writerow(dict(group=group, pooled_fragment=pooled, sample=sample,
                            biological_sample=source['biological_sample'], original_fragment=original))
                        counts[group] += 1
                if set(selected[sample])-visited:
                    raise ValueError('Selected coassembly fragments are absent from reads')
        table(out/'pools.tsv', [dict(group=group, pairs=counts[group],
            biological_samples=len({samples[n]['biological_sample'] for n in names}),
            r1=f'{group}/reads_R1.fastq.gz', r2=f'{group}/reads_R2.fastq.gz') for group,names in sorted(group_samples.items())],
            ['group','pairs','biological_samples','r1','r2'])
        dump(out/'summary.json', dict(groups=len(group_samples), pairs=dict(counts), sequence_joined=False,
            next_step='Run existing assemble on each pool and compare against per-sample assemblies'))
