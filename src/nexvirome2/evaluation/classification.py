from collections import defaultdict

def aligned_bases(hit):
    q, s = hit['qstart']-1, hit['sstart']-1
    dq = 1 if hit['qend'] >= hit['qstart'] else -1
    ds = 1 if hit['send'] >= hit['sstart'] else -1
    for qb, sb in zip(hit['qseq'], hit['sseq']):
        if qb != '-' and sb != '-':
            yield q, s, qb.upper() == sb.upper() and qb.upper() in 'ACGT'
        q += dq if qb != '-' else 0
        s += ds if sb != '-' else 0


def aligned_runs(hit):
    """Blocks continuous in both coordinates; alignment gaps are not recovered bases."""
    dq = 1 if hit['qend'] >= hit['qstart'] else -1
    ds = 1 if hit['send'] >= hit['sstart'] else -1
    run = []
    for q,s,match in aligned_bases(hit):
        if run and (q != run[-1][0]+dq or s != run[-1][1]+ds):
            yield run
            run = []
        run.append((q,s,match))
    if run:
        yield run


def classify(sequence, hits, minimum=150, identity=99.0, explained=0.99):
    if minimum < 1 or not 0 <= identity <= 100:
        raise ValueError('Invalid alignment thresholds')
    support = [set() for _ in sequence]
    exact = [set() for _ in sequence]
    valid = [h for h in hits if h['pident'] >= identity and h['length'] >= minimum]
    by_source = defaultdict(set)
    for h in valid:
        for q, _, match in aligned_bases(h):
            if not 0 <= q < len(sequence):
                raise ValueError('Alignment outside contig')
            support[q].add(h['sseqid'])
            by_source[h['sseqid']].add(q)
            if match:
                exact[q].add(h['sseqid'])
    if not by_source:
        return {'status': 'unaligned', 'breakpoints': [], 'valid_hits': valid}
    complete = [s for s, positions in by_source.items() if len(positions)/len(sequence) >= explained]
    unique = [next(iter(x)) if len(x) == 1 else None for x in exact]
    if len(complete) > 1:
        return {'status': 'ambiguous', 'breakpoints': [], 'valid_hits': valid}
    if len(complete) == 1:
        # One full explanation wins over a spurious split best-hit chain.
        status = 'consistent' if complete[0] in unique else 'ambiguous'
        return {'status': status, 'source': complete[0], 'breakpoints': [], 'valid_hits': valid}
    blocks = []
    for source in sorted(by_source):
        positions = sorted(by_source[source])
        start = previous = positions[0]
        for p in positions[1:] + [len(sequence)+1]:
            if p != previous+1:
                identifiers = [i for i in range(start, previous+1) if unique[i] == source]
                if previous+1-start >= minimum and identifiers:
                    blocks.append((start, previous+1, source, min(identifiers), max(identifiers)))
                start = p
            previous = p
    blocks.sort(key=lambda b: (b[3], b[4], b[2]))
    breakpoints = []
    for left, right in zip(blocks, blocks[1:]):
        if left[2] != right[2] and left[4] < right[3] and left[1] >= right[0]:
            breakpoints.append({'left_source': left[2], 'right_source': right[2],
                                'start': left[4]+1, 'end': right[3]})
        elif left[2] != right[2] and left[4] < right[3] and right[0]-left[1] <= 10:
            breakpoints.append({'left_source': left[2], 'right_source': right[2],
                                'start': left[4]+1, 'end': right[3]})
    return {'status': 'chimeric' if breakpoints else 'ambiguous', 'breakpoints': breakpoints, 'valid_hits': valid}


