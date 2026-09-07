"""Oriented GFA links and SPAdes contig paths with checked overlap spelling."""
import re
from collections import defaultdict

from ..common import revcomp


def flip(node):
    return node[:-1] + ('-' if node[-1] == '+' else '+')


class Graph:
    def __init__(self):
        self.sequences = {}
        self.links = defaultdict(dict)
        self.coverage = {}

    def sequence(self, node):
        seq = self.sequences[node[:-1]]
        return seq if node[-1] == '+' else revcomp(seq)

    @classmethod
    def load(cls, path):
        graph = cls()
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip().split('\t')
                if fields[0] == 'S':
                    if fields[2] == '*':
                        raise ValueError('GFA must include segment sequences')
                    graph.sequences[fields[1]] = fields[2].upper()
                    for tag in fields[3:]:
                        if tag.startswith('DP:f:') or tag.startswith('DP:i:'):
                            graph.coverage[fields[1]] = float(tag[5:])
                elif fields[0] == 'L':
                    match = re.fullmatch(r'(\d+)M', fields[5])
                    if not match:
                        raise ValueError('Only exact GFA overlaps (nM) are supported')
                    a, b, overlap = fields[1]+fields[2], fields[3]+fields[4], int(match[1])
                    graph.links[a][b] = overlap
                    graph.links[flip(b)][flip(a)] = overlap
        for a, edges in list(graph.links.items()):
            for b, overlap in edges.items():
                if overlap > min(len(graph.sequence(a)), len(graph.sequence(b))):
                    raise ValueError('GFA overlap exceeds segment length')
                if overlap and graph.sequence(a)[-overlap:] != graph.sequence(b)[:overlap]:
                    raise ValueError(f'Inconsistent GFA overlap: {a} -> {b}')
        return graph

    def spell(self, path):
        seq = self.sequence(path[0])
        offsets = [0]
        for a, b in zip(path, path[1:]):
            overlap = self.links[a][b]
            offsets.append(len(seq)-overlap)
            seq += self.sequence(b)[overlap:]
        return seq, offsets

    def incoming(self, node):
        return sorted(flip(n) for n in self.links.get(flip(node), {}))


def contig_paths(path):
    result = {}
    name = None
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if not re.fullmatch(r'[0-9]+[+-](?:[,;][0-9]+[+-])*;?', line):
                name = line
                result[name] = []
            elif name is None:
                raise ValueError('Path before contig name')
            else:
                # Semicolons describe gaps; retain a sentinel, never invent a link.
                chunks = line.rstrip(';').split(';')
                if result[name]:
                    result[name].append(None)
                for i, chunk in enumerate(chunks):
                    if i:
                        result[name].append(None)
                    result[name].extend(chunk.split(','))
    return result


def candidates(graph, sequences, paths, max_repeat_nodes=20):
    """Find an entrance merge and downstream exit branch on each contig path."""
    result = []
    for contig, sequence in sequences.items():
        path = paths.get(contig)
        if not path or None in path:
            result.append({'contig': contig, 'reason': 'missing_or_gapped_path', 'usable': False})
            continue
        try:
            spelled, offsets = graph.spell(path)
        except (KeyError, ValueError):
            result.append({'contig': contig, 'reason': 'invalid_path', 'usable': False})
            continue
        if spelled != sequence:
            result.append({'contig': contig, 'reason': 'path_sequence_mismatch', 'usable': False})
            continue
        seen = set()
        for end in range(1, len(path)-1):
            if len(graph.links[path[end]]) < 2:
                continue
            entrances = [start for start in range(max(1, end-max_repeat_nodes+1), end+1)
                         if len(graph.incoming(path[start])) > 1]
            if not entrances:
                continue
            start = entrances[-1]
            current = path[start-1:end+2]
            if len(set(current)) < len(current):
                continue
            alternatives = [[*current[:-1], other] for other in sorted(graph.links[path[end]]) if other != current[-1]]
            cut = offsets[end+1] + graph.links[path[end]][path[end+1]]
            if (contig, cut) in seen or not 0 < cut < len(sequence):
                continue
            seen.add((contig, cut))
            repeat, _ = graph.spell(current[1:-1])
            result.append({'contig': contig, 'cut': cut, 'current': current, 'alternatives': alternatives,
                           'repeat_length': len(repeat), 'usable': True, 'reason': 'candidate'})
    return result
