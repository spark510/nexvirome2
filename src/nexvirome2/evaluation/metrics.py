from collections import defaultdict
from ..io.formats import union_length
from .classification import aligned_runs

def summarize_metrics(calls,outcomes,references,settings):
    metrics = {}
    for cutoff in sorted({0, settings.get('min_contig_length', 500)}):
        included = [r for r in calls if r['length'] >= cutoff]
        coverage, longest, consistent_longest = defaultdict(list), defaultdict(int), defaultdict(int)
        for row in included:
            for h in outcomes[row['contig']]['valid_hits']:
                for run in aligned_runs(h):
                    positions = [s for _,s,_ in run]
                    if min(positions) < 0 or max(positions) >= len(references[h['sseqid']]):
                        raise ValueError('Subject alignment outside source')
                    coverage[h['sseqid']].append((min(positions), max(positions)+1))
                    longest[h['sseqid']] = max(longest[h['sseqid']], len(positions))
                    call = outcomes[row['contig']]
                    if call['status'] == 'consistent' and call.get('source') == h['sseqid']:
                        consistent_longest[h['sseqid']] = max(consistent_longest[h['sseqid']], len(positions))
        per_source = {s: {'recovered_bp': union_length(coverage[s]), 'recovery': union_length(coverage[s])/len(seq),
                          'longest_aligned_block_bp': longest[s],
                          'longest_consistent_block_bp': consistent_longest[s]} for s, seq in references.items()}
        total = sum(len(s) for s in references.values())
        lengths = sorted((r['length'] for r in included), reverse=True)
        accumulated, n50 = 0, 0
        for length in lengths:
            accumulated += length
            if accumulated >= sum(lengths)/2:
                n50 = length
                break
        metrics[str(cutoff)] = {'contigs': len(included), 'chimeric_contigs': sum(r['status']=='chimeric' for r in included),
                               'wrong_connections': sum(r['wrong_connections'] for r in included),
                               'chimera_rate': sum(r['status']=='chimeric' for r in included)/len(included) if included else 0,
                               'genome_recovery': sum(v['recovered_bp'] for v in per_source.values())/total if total else 0,
                               'per_source': per_source, 'total_contig_bp': sum(lengths), 'n50': n50,
                               'longest_aligned_block_bp': max(longest.values(), default=0),
                               'correctly_reconstructed_fraction': sum(min(consistent_longest[s], len(seq)) for s, seq in references.items())/total if total else 0,
                               'ambiguous_contigs': sum(r['status']=='ambiguous' for r in included)}
    return metrics
