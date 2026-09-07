"""Conservative near-identity clustering: every member must match its representative."""
from collections import defaultdict
from .common import fasta, stage, table, dump, write_fasta, revcomp
from .alignment import blast, read_hits


def cluster(settings, contigs, output, alignments=None):
    identity = float(settings.get('min_identity', 95))
    coverage = float(settings.get('min_coverage', .95))
    if not 0 <= identity <= 100 or not 0 < coverage <= 1:
        raise ValueError('Invalid clustering identity/coverage thresholds')
    with stage(output, settings, [contigs] + ([alignments] if alignments else [])) as (out, manifest):
        sequences = fasta(contigs)
        if not sequences or any(not seq for seq in sequences.values()):
            raise ValueError('Catalog sequences must be nonempty')
        edges = defaultdict(dict)
        hits = read_hits(alignments or blast(contigs, contigs, out, manifest, threads=settings.get('threads', 1)))
        for hit in hits:
            query, target = hit['qseqid'], hit['sseqid']
            if query not in sequences or target not in sequences:
                raise ValueError('Unknown clustering alignment ID')
            if len(hit['qseq']) != len(hit['sseq']) or len(hit['qseq']) != hit['length']:
                raise ValueError('Inconsistent clustering alignment strings')
            for accession, prefix in ((query, 'q'), (target, 's')):
                start, end = hit[prefix+'start'], hit[prefix+'end']
                if not 1 <= min(start, end) <= max(start, end) <= len(sequences[accession]):
                    raise ValueError('Clustering alignment outside contig')
                expected = sequences[accession][min(start, end)-1:max(start, end)]
                if start > end:
                    expected = revcomp(expected)
                if expected != hit[prefix+'seq'].replace('-', '').upper():
                    raise ValueError('Clustering alignment differs from original sequence')
            pairs = list(zip(hit['qseq'].upper(), hit['sseq'].upper()))
            aligned = sum(a != '-' and b != '-' for a, b in pairs)
            matches = sum(a == b and a in 'ACGT' for a, b in pairs)
            measured_identity = 100*matches/len(pairs) if pairs else 0
            mutual = min(aligned/len(sequences[query]), aligned/len(sequences[target]))
            if query == target or measured_identity < identity or mutual < coverage:
                continue
            direction = '+' if (hit['qend'] >= hit['qstart']) == (hit['send'] >= hit['sstart']) else '-'
            evidence = dict(identity=measured_identity, coverage=mutual, orientation=direction)
            for a, b in ((query, target), (target, query)):
                if b not in edges[a] or (measured_identity, mutual, direction) > tuple(edges[a][b][k] for k in ('identity','coverage','orientation')):
                    edges[a][b] = evidence
        remaining, representatives, lineage = set(sequences), {}, []
        for representative in sorted(sequences, key=lambda name: (-len(sequences[name]), name)):
            if representative not in remaining:
                continue
            representatives[representative] = sequences[representative]
            members = {representative} | (set(edges[representative]) & remaining)
            for member in sorted(members):
                evidence = edges[representative].get(member, dict(identity=100, coverage=1, orientation='+'))
                lineage.append(dict(representative=representative, member=member, **evidence))
            remaining -= members
        write_fasta(out/'representatives.fasta', representatives)
        table(out/'clusters.tsv', lineage, ['representative','member','identity','coverage','orientation'])
        dump(out/'summary.json', dict(input_contigs=len(sequences), representatives=len(representatives),
             method='greedy_representative_single_HSP_mutual_coverage',
             limitations=['No transitive chain merging', 'No circular rotation equivalence', 'No multi-HSP chaining']))
