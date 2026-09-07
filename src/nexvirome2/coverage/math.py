from collections import Counter
from ..io.formats import revcomp

def composition(sequence, k=4):
    counts = Counter()
    for i in range(len(sequence)-k+1):
        word = sequence[i:i+k]
        if set(word) <= set('ACGT'):
            counts[min(word, revcomp(word))] += 1
    total = sum(counts.values())
    return {key: value/total for key,value in counts.items()} if total else {}


def cosine(a, b):
    import numpy as np
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    norm = float(np.linalg.norm(a)*np.linalg.norm(b))
    return float(np.clip(a@b/norm, -1, 1)) if norm else None


def merged(intervals):
    result = []
    for a,b in sorted(intervals):
        if result and a <= result[-1][1]:
            result[-1] = (result[-1][0], max(b,result[-1][1]))
        else:
            result.append((a,b))
    return result


def windows(sequences, size=500, step=500, regions=None):
    if size < 1 or step < 1 or step > size:
        raise ValueError('Require 0 < step <= window size')
    result = []
    for name, seq in sequences.items():
        for start in range(0, len(seq), step):
            end = min(start+size, len(seq))
            if regions is not None and not any(start < b and end > a for a,b in regions.get(name, [])):
                continue
            result.append({'window_id': f'w{len(result):09d}', 'contig': name, 'start': start, 'end': end})
            if end == len(seq):
                break
    return result


