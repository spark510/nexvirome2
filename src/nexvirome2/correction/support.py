import re
import math
import json
import sqlite3
import tempfile
import itertools
from pathlib import Path
from contextlib import closing, contextmanager, nullcontext
from collections import defaultdict
from ..io.sam import blocks, sam_records

def cigar_length(cigar):
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar) if op in 'MDN=X')


def validate_score_margin(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError('score_margin must be finite and greater than zero')
    return value


@contextmanager
def fragment_records(sam):
    """Group arbitrary-order SAM on disk; retain only one fragment's hits in RAM."""
    with tempfile.TemporaryDirectory(prefix='nexvirome2-pairs-') as directory:
        with closing(sqlite3.connect(str(Path(directory)/'records.sqlite'))) as database:
            database.execute('PRAGMA temp_store=FILE')
            database.execute('PRAGMA cache_size=-4096')
            database.execute('CREATE TABLE records (name TEXT, record TEXT)')
            database.executemany('INSERT INTO records VALUES (?, ?)',
                                 ((f[0], json.dumps(f)) for f in sam_records(sam)))
            database.execute('CREATE INDEX fragment_name ON records(name)')
            cursor = database.execute('SELECT name, record FROM records ORDER BY name')
            yield ((name, [json.loads(row[1]) for row in group])
                   for name, group in itertools.groupby(cursor, key=lambda row:row[0]))


def pair_evidence(sam, bounds, score_margin=1, audit_file=None):
    """Count independent, uniquely spanning fragments from arbitrary-order SAM.

    Use temporary disk storage for grouping and stream the optional JSON audit.
    Endpoint sets remain in memory to deduplicate supporting fragments.
    """
    validate_score_margin(score_margin)
    with fragment_records(sam) as groups, (open(audit_file, 'w', encoding='utf-8') if audit_file is not None else nullcontext()) as audit:
        if audit is not None:
            audit.write('{"fragments":[')
        first_item = True
        def record(item):
            nonlocal first_item
            if audit is not None:
                if not first_item:
                    audit.write(',')
                json.dump(item, audit)
                first_item = False
        result = _pair_evidence(groups, bounds, score_margin, record)
        if audit is not None:
            audit.write('],"interpretation":"Ambiguous candidates are retained, not counted as unique support"}\n')
        return result


def _pair_evidence(groups, bounds, score_margin, record):
    """Count distinct fragment endpoints only when a unique best template is spanned.

    bounds: template -> (last upstream-anchor base exclusive, first downstream-anchor base).
    Alignments to common repeat sequence alone cannot support an entrance/exit choice.
    """
    support = defaultdict(set)
    ambiguous = 0
    for read, records in groups:
        if any(int(f[1]) & (512 | 1024) for f in records):
            record({'fragment':read,'status':'excluded_qc_or_duplicate','candidates':[]})
            continue
        templates = defaultdict(list)
        for f in records:
            flag = int(f[1])
            if flag & (4 | 8 | 2048) or not flag & 2 or f[2] not in bounds:
                continue
            tags = {v.split(':', 2)[0]: v.split(':', 2)[2] for v in f[11:]}
            if 'AS' not in tags:
                raise ValueError('Bowtie2 alignment lacks AS score')
            if (flag & (64 | 128)) not in (64, 128):
                continue
            templates[f[2]].append({'mate': 1 if flag & 64 else 2,
                                         'pos': int(f[3])-1, 'mate_pos': int(f[7])-1,
                                         'cigar': f[5],
                                         'blocks': list(blocks(int(f[3])-1,f[5])),
                                         'end': int(f[3])-1+cigar_length(f[5]), 'score': int(tags['AS']),
                                         'reverse': bool(flag & 16)})
        scored = []
        for template, hits in templates.items():
            left_bound, right_bound = bounds[template]
            for first in [h for h in hits if h['mate'] == 1]:
                for second in [h for h in hits if h['mate'] == 2]:
                    if first['mate_pos'] != second['pos'] or second['mate_pos'] != first['pos']:
                        continue
                    left, right = sorted((first, second), key=lambda h: h['pos'])
                    # Require anchor overlap by >= 20 bp on both sides.
                    upstream = sum(max(0,min(end,left_bound)-start) for start,end in left['blocks'])
                    downstream = sum(max(0,end-max(start,right_bound)) for start,end in right['blocks'])
                    spans = upstream >= 20 and downstream >= 20
                    # Alignment identity is distinct from endpoint-based PCR
                    # deduplication. Different CIGARs/mate placements compete.
                    placement = tuple((h['pos'], h['cigar'], h['reverse']) for h in (first, second))
                    scored.append((first['score']+second['score'], template,
                                   (left['pos'], right['end'], first['reverse']), spans, placement))
        if not scored:
            record({'fragment':read,'status':'no_valid_pair','candidates':[]})
            continue
        scored = sorted(set(scored), reverse=True)
        best = scored[0]
        competitors = [s for s in scored[1:] if (s[1],s[4]) != (best[1],best[4])]
        item = {'fragment':read,'status':'unique_spanning' if best[3] else 'nonspanning',
                'candidates':[{'score':s[0],'path':s[1],'endpoints':s[2],'spanning':s[3],
                               'mate_placements':s[4]} for s in scored]}
        if competitors and best[0]-competitors[0][0] < score_margin:
            item['status'] = 'ambiguous'
            ambiguous += 1
            record(item)
            continue
        record(item)
        if best[3]:
            support[best[1]].add(best[2])
    return {name: len(support[name]) for name in bounds}, ambiguous


