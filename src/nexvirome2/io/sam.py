import re
import itertools
import math
from collections import defaultdict
from ..coverage.math import merged

def blocks(position, cigar):
    cursor = position
    if ''.join(n+op for n,op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)) != cigar:
        raise ValueError('Invalid CIGAR')
    for count, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(count)
        if op in 'M=X':
            yield cursor, cursor+length
        if op in 'MDN=X':
            cursor += length


def sam_records(path):
    with open(path) as handle:
        for line in handle:
            if not line.startswith('@'):
                f = line.rstrip().split('\t')
                if len(f) < 11:
                    raise ValueError('Invalid SAM record')
                yield f


def alignment_score(record, tag='AS'):
    tags = {v.split(':',2)[0]:v.split(':',2)[2] for v in record[11:]}
    try:
        value = float(tags[tag]) if tag in tags else None
    except (ValueError, TypeError):
        return None
    return value if value is not None and math.isfinite(value) else None


def competing_alignment(primary, records):
    """Require explicit lower scores when another placement is present.

    Per-mate comparison is conservative: a tied mate is not unique evidence even
    if its partner favors the primary placement. Supplementary mappings also
    remain outside this unique-pair coverage calculation.
    """
    originals = {int(f[1]) & (64|128): f for f in primary}
    for record in records:
        flag = int(record[1])
        if flag & 4:
            continue
        if flag & 2048:
            return True
        if not flag & 256:
            continue
        original = originals.get(flag & (64|128))
        if original is None:
            return True
        placement = lambda f: (f[2],f[3],f[5],int(f[1]) & 16)
        if placement(record) == placement(original):
            continue
        best, alternative = alignment_score(original), alignment_score(record)
        if best is None or alternative is None or alternative >= best:
            return True
    return False


def measure_sam(path, target, regions, min_mapq=20):
    """Input is name-sorted SAM. Count each accepted pair's overlapping bases once."""
    header = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith('@'):
                break
            if line.startswith('@SQ\t'):
                tags = dict(field.split(':',1) for field in line.strip().split('\t')[1:])
                header[tags['SN']] = int(tags['LN'])
    if header != {a:len(s) for a,s in target.items()}:
        raise ValueError('SAM reference names/lengths differ from common target')
    by_contig = defaultdict(list)
    for i,row in enumerate(regions):
        by_contig[row['contig']].append((i,row['start'],row['end']))
    counts = [0]*len(regions)
    stats = {'accepted_fragments':0, 'excluded_fragments':0, 'aligned_bases':0}
    for _,group in itertools.groupby(sam_records(path), key=lambda f:f[0]):
        records = list(group)
        primary = [f for f in records if not int(f[1]) & (256|2048)]
        valid = len(primary) == 2 and all(int(f[1]) & 2 and not int(f[1]) & (4|8|512|1024)
                                        and min_mapq <= int(f[4]) < 255 for f in primary)
        if valid:
            valid = primary[0][2] == primary[1][2] and {int(f[1]) & (64|128) for f in primary} == {64,128}
            for f in primary:
                tags = {v.split(':',2)[0]:v.split(':',2)[2] for v in f[11:]}
                if 'XS' in tags:
                    best, alternative = alignment_score(f), alignment_score(f, 'XS')
                    if best is None or alternative is None or alternative >= best:
                        valid = False
            if valid and competing_alignment(primary, records):
                valid = False
        if not valid:
            stats['excluded_fragments'] += 1
            continue
        contig = primary[0][2]
        intervals = merged([b for f in primary for b in blocks(int(f[3])-1,f[5])])
        if any(a < 0 or b > len(target[contig]) for a,b in intervals):
            raise ValueError('Alignment outside target')
        stats['accepted_fragments'] += 1
        stats['aligned_bases'] += sum(b-a for a,b in intervals)
        for i,left,right in by_contig[contig]:
            counts[i] += sum(max(0,min(right,b)-max(left,a)) for a,b in intervals)
    return [n/(r['end']-r['start']) for n,r in zip(counts,regions)], stats


