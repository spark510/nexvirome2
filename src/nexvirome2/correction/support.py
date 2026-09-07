import re
from collections import defaultdict
from ..io.sam import blocks

def cigar_length(cigar):
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar) if op in 'MDN=X')


def pair_evidence(sam, bounds, score_margin=1, audit_file=None):
    """Count distinct fragment endpoints only when a unique best template is spanned.

    bounds: template -> (last upstream-anchor base exclusive, first downstream-anchor base).
    Alignments to common repeat sequence alone cannot support an entrance/exit choice.
    """
    alignments = defaultdict(lambda: defaultdict(list))
    with open(sam) as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            f = line.rstrip().split('\t')
            flag = int(f[1])
            if flag & (4 | 8 | 2048) or not flag & 2 or f[2] not in bounds:
                continue
            tags = {v.split(':', 2)[0]: v.split(':', 2)[2] for v in f[11:]}
            if 'AS' not in tags:
                raise ValueError('Bowtie2 alignment lacks AS score')
            alignments[f[0]][f[2]].append({'mate': 1 if flag & 64 else 2,
                                         'pos': int(f[3])-1, 'mate_pos': int(f[7])-1,
                                         'blocks': list(blocks(int(f[3])-1,f[5])),
                                         'end': int(f[3])-1+cigar_length(f[5]), 'score': int(tags['AS']),
                                         'reverse': bool(flag & 16)})
    support = defaultdict(set)
    ambiguous = 0
    audit = []
    for read, templates in alignments.items():
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
                    scored.append((first['score']+second['score'], template,
                                   (left['pos'], right['end'], first['reverse']), spans))
        if not scored:
            audit.append({'fragment':read,'status':'no_valid_pair','candidates':[]})
            continue
        scored = sorted(set(scored), reverse=True)
        best = scored[0]
        competitors = [s for s in scored[1:] if s[1:3] != best[1:3]]
        item = {'fragment':read,'status':'unique_spanning' if best[3] else 'nonspanning',
                'candidates':[{'score':s[0],'path':s[1],'endpoints':s[2],'spanning':s[3]} for s in scored]}
        audit.append(item)
        if competitors and best[0]-competitors[0][0] < score_margin:
            item['status'] = 'ambiguous'
            ambiguous += 1
            continue
        if best[3]:
            support[best[1]].add(best[2])
    if audit_file is not None:
        from ..common import dump
        dump(audit_file, {'fragments':audit,'interpretation':'Ambiguous candidates are retained, not counted as unique support'})
    return {name: len(support[name]) for name in bounds}, ambiguous


