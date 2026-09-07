"""Explicit diagnostic classifier and external competitive-hit adapter.

This is not a replacement claimed to reproduce legacy NexVirome semantics.
Query IDs are the evaluation unit; paired fragments must be supplied as such
by the external adapter. Never simulate from a masked reference view.
"""
from collections import defaultdict, Counter
from dataclasses import dataclass
from typing import Protocol
import math
from .common import fasta, rows, stage, table, dump, sha
from .taxonomy import TaxonomyResolver
from .masking import MaskSet
from .sharing import seeds


class ClassifierBackend(Protocol):
    def hits(self, queries, reference): ...


@dataclass
class ExactSeedBackend:
    k: int = 31

    def hits(self, queries, reference):
        if self.k < 3:
            raise ValueError('k must be at least 3')
        index = defaultdict(set)
        for accession,sequence in reference.items():
            for _,word in seeds(sequence,self.k):
                index[word].add(accession)
        result = []
        for query,sequence in queries.items():
            counts = Counter()
            for word in {word for _,word in seeds(sequence,self.k)}:
                counts.update(index.get(word, ()))
            result.extend(dict(query=query,target=target,score=float(score)) for target,score in counts.items())
        return result


def classify(settings, queries, reference, taxonomy, output, mask_file=None, hits_file=None):
    from .decisions import fingerprint
    rank = settings.get('rank','species')
    minimum, fraction = float(settings.get('min_score',1)), float(settings.get('competitive_fraction',1))
    if not math.isfinite(minimum) or minimum<=0 or not 0<fraction<=1:
        raise ValueError('Invalid classifier score thresholds')
    with stage(output,settings,[p for p in (queries,reference,taxonomy,mask_file,hits_file) if p]) as (out,_):
        query_sequences, reference_sequences = fasta(queries), fasta(reference)
        resolver = TaxonomyResolver.load(taxonomy)
        for acc in reference_sequences:
            resolver.taxon(acc,rank)
        if hits_file and mask_file:
            raise ValueError('External hits must declare their masked view via config; do not apply a second mask')
        if mask_file:
            mask = MaskSet.load(mask_file,reference)
            if mask.policy['rank'] != rank:
                raise ValueError('Classifier and mask ranks differ')
            reference_sequences = mask.render(reference,'hard')
        records = rows(hits_file) if hits_file else ExactSeedBackend(int(settings.get('k',31))).hits(query_sequences,reference_sequences)
        per_query = defaultdict(dict)
        for row in records:
            query,target,score = row['query'],row['target'],float(row['score'])
            if query not in query_sequences or target not in reference_sequences:
                raise ValueError('Classifier hit has an unknown query/target')
            if not math.isfinite(score) or score<0:
                raise ValueError('Classifier scores must be finite and nonnegative')
            per_query[query][target] = max(score,per_query[query].get(target,0))
        predictions, candidates = [], []
        for query in query_sequences:
            scores = per_query[query]
            best = max(scores.values(),default=0)
            selected = sorted(acc for acc,score in scores.items() if score>=minimum and score>=best*fraction)
            known = sorted({resolver.taxon(acc,rank) for acc in selected if resolver.taxon(acc,rank) is not None})
            unknown = any(resolver.taxon(acc,rank) is None for acc in selected)
            status = 'unclassified' if not selected else 'assigned' if len(known)==1 and not unknown else 'ambiguous'
            predictions.append(dict(query=query,status=status,taxon=known[0] if status=='assigned' else '',
                                    candidate_taxa=';'.join(known),unknown_taxonomy=str(unknown).lower()))
            candidates.extend(dict(query=query,target=acc,score=score,selected=str(acc in selected).lower())
                              for acc,score in sorted(scores.items()))
        table(out/'predictions.tsv',predictions,['query','status','taxon','candidate_taxa','unknown_taxonomy'])
        table(out/'candidates.tsv',candidates,['query','target','score','selected'])
        dump(out/'classification.json', dict(schema=1,rank=rank,queries_sha256=sha(queries),
             reference_sha256=sha(reference),taxonomy_sha256=sha(taxonomy),
             mask_sha256=sha(mask_file) if mask_file else None,
             backend='external_competitive_hits' if hits_file else 'diagnostic_exact_seed',
             unit=settings.get('unit','query'),conditions=settings.get('conditions',{}),
             predictions_sha256=sha(out/'predictions.tsv'),candidates_sha256=sha(out/'candidates.tsv'),
             settings_sha256=fingerprint(settings)))


def evaluate(settings, predictions, truth, output, metadata):
    """Truth TSV: query, taxon, unit (independent source/pair/sample cluster)."""
    from .common import config
    with stage(output,settings,[predictions,truth,metadata]) as (out,_):
        info = config(metadata)
        if info['rank']!=settings['rank'] or info['predictions_sha256']!=sha(predictions):
            raise ValueError('Classification metadata rank/checksum mismatch')
        def unique(path):
            result = {}
            for row in rows(path):
                if row['query'] in result:
                    raise ValueError('Duplicate query ID')
                result[row['query']] = row
            return result
        observed, expected = unique(predictions), unique(truth)
        if observed.keys()!=expected.keys() or not expected:
            raise ValueError('Predictions must include every truth query exactly once')
        groups = defaultdict(Counter)
        for query,row in expected.items():
            if not row['taxon'] or not row['unit']:
                raise ValueError('Truth needs an explicit taxon and independent unit')
            pred = observed[query]
            if pred['status'] not in ('assigned','ambiguous','unclassified'):
                raise ValueError('Unknown classification status')
            if (pred['status']=='assigned') != bool(pred['taxon']):
                raise ValueError('Assigned status and taxon disagree')
            label = ('correct' if pred['taxon']==row['taxon'] else 'wrong') if pred['status']=='assigned' else pred['status']
            groups[row['unit']]['total'] += 1
            groups[row['unit']][label] += 1
        totals = sum(groups.values(),Counter())
        def metrics(counts):
            n = counts['total']
            assigned = counts['correct']+counts['wrong']
            return dict(total=n,**{key:counts[key]/n for key in ('correct','wrong','ambiguous','unclassified')},
                        precision=counts['correct']/assigned if assigned else None)
        table(out/'units.tsv',[dict(unit=unit,**metrics(counts)) for unit,counts in sorted(groups.items())],
              ['unit','total','correct','wrong','ambiguous','unclassified','precision'])
        dump(out/'metrics.json', dict(schema=1,rank=info['rank'],**metrics(totals),
             independent_units=len(groups),classification=info,truth_sha256=sha(truth),
             interpretation='Correct/wrong/ambiguous/unclassified fractions use all queries; precision uses assigned queries'))
