from collections import defaultdict
from ..alignment import blast, read_hits
from ..common import dump, fasta, stage, table
from .classification import classify
from .metrics import summarize_metrics

def evaluate(settings, contigs, sources, output, alignment_file=None):
    inputs = [contigs, sources] + ([alignment_file] if alignment_file else [])
    with stage(output, settings, inputs) as (out, manifest):
        sequences, references = fasta(contigs), fasta(sources)
        hits = read_hits(alignment_file or blast(contigs, sources, out, manifest, threads=settings.get('threads', 1)))
        grouped = defaultdict(list)
        for hit in hits:
            if hit['qseqid'] not in sequences or hit['sseqid'] not in references:
                raise ValueError('Unknown alignment ID')
            grouped[hit['qseqid']].append(hit)
        calls, breakpoints, outcomes = [], [], {}
        for name, seq in sequences.items():
            result = classify(seq, grouped[name], settings.get('min_segment_bp', 150), settings.get('min_identity', 99))
            outcomes[name] = result
            calls.append({'contig': name, 'length': len(seq), 'status': result['status'], 'wrong_connections': len(result['breakpoints'])})
            breakpoints.extend({'contig': name, **b} for b in result['breakpoints'])
        table(out / 'contig_calls.tsv', calls, ['contig', 'length', 'status', 'wrong_connections'])
        table(out / 'breakpoints.tsv', breakpoints, ['contig', 'left_source', 'right_source', 'start', 'end'])
        metrics = summarize_metrics(calls,outcomes,references,settings)
        dump(out / 'metrics.json', metrics)


