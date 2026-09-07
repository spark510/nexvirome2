"""Retain full-reference support without forcing masked-read fine taxonomy."""
from collections import defaultdict
from .common import config,rows,stage,table,dump,sha


def link(settings,masked_run,unmasked_run,output):
    from pathlib import Path
    masked,full=Path(masked_run),Path(unmasked_run)
    paths=[folder/name for folder in (masked,full) for name in ('classification.json','predictions.tsv','candidates.tsv')]
    with stage(output,settings,paths) as (out,_):
        a,b=config(masked/'classification.json'),config(full/'classification.json')
        for field in ('rank','queries_sha256','reference_sha256','taxonomy_sha256','unit','conditions','backend','settings_sha256'):
            if a.get(field)!=b.get(field):
                raise ValueError('Conserved reuse requires matched queries, database and conditions')
        if b.get('mask_sha256') is not None or not a.get('mask_sha256'):
            raise ValueError('Expected masked and unmasked classification runs')
        for folder,info in ((masked,a),(full,b)):
            if info['predictions_sha256']!=sha(folder/'predictions.tsv'):
                raise ValueError('Classification predictions checksum mismatch')
            if info.get('candidates_sha256')!=sha(folder/'candidates.tsv'):
                raise ValueError('Classification candidates checksum mismatch')
        def unique(path):
            source=rows(path)
            result={r['query']:r for r in source}
            if len(result)!=len(source):
                raise ValueError('Duplicate classification query')
            return result
        predictions=unique(masked/'predictions.tsv')
        original=unique(full/'predictions.tsv')
        if predictions.keys()!=original.keys():
            raise ValueError('Conserved support query IDs disagree')
        unique=defaultdict(set)
        for row in rows(masked/'candidates.tsv'):
            if row['selected']=='true':
                unique[row['query']].add(row['target'])
        support=[]
        for row in rows(full/'candidates.tsv'):
            query=row['query']
            if query not in predictions:
                raise ValueError('Unknown conserved-support query')
            # Candidates outside the masked evidence are preserved, never silently removed.
            support.append(dict(query=query,target=row['target'],full_score=row['score'],
                 unmasked_competitor=row['selected'],masked_supported=str(row['target'] in unique[query]).lower(),
                 fine_status=predictions[query]['status'],fine_taxon=predictions[query]['taxon'],
                 interpretation='annotation_or_detection_support_only'))
        table(out/'reference_support.tsv',support,['query','target','full_score','unmasked_competitor',
              'masked_supported','fine_status','fine_taxon','interpretation'])
        table(out/'predictions.tsv',predictions.values(),['query','status','taxon','candidate_taxa','unknown_taxonomy'])
        dump(out/'summary.json',dict(fine_assignments_changed=False,all_full_reference_candidates_retained=True,
             support_rows=len(support),conserved_only_queries=sum(original[q]['status']!='unclassified' and
             r['status']=='unclassified' for q,r in predictions.items()),
             note='Read support is not a count of independent taxa or proof of viral origin'))
