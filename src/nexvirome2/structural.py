"""Structural annotation application facade and compatibility exports."""
from collections import defaultdict
from .common import rows,stage,table
from .search.formats import read_search
from .evidence.structural import summarize,annotate_orfs,summarize_contigs
from .evidence.proposals import proposals
from .correction.policy import gate

PROPOSAL_FIELDS = ['contig','start','end','method','score','detail']

def annotation(settings,orf_file,metadata,output,sequence_hits=None,structure_hits=None):
    if not sequence_hits and not structure_hits:
        raise ValueError('Provide sequence and/or structural search evidence')
    inputs = [orf_file,metadata]+[p for p in (sequence_hits,structure_hits) if p]
    with stage(output,settings,inputs) as (out,manifest):
        orfs,meta = rows(orf_file),rows(metadata)
        labels = {r['target']:r for r in meta}
        if len(labels)!=len(meta) or any(r.get('origin') not in ('viral','cellular','uncertain') for r in meta):
            raise ValueError('Target labels must be unique and origin viral/cellular/uncertain')
        ids = {r['protein_id'] for r in orfs}
        evidence = []
        for path in (sequence_hits,structure_hits):
            grouped = defaultdict(list)
            if path:
                for h in read_search(path):
                    if h['query'] not in ids:
                        raise ValueError('Search protein not present in ORF map')
                    grouped[h['query']].append(h)
            evidence.append(grouped)
        profiles,by_contig = annotate_orfs(orfs,evidence,labels,settings)
        fields = ['orf_id','protein_id','contig','start','end','strand','origin','lineage','sequence_family','structure_family',
                  'sequence_related','sequence_close','structural_related','structural_only','hallmark','confidence']
        table(out/'orf_profiles.tsv',profiles,fields)
        changes,novelty = summarize_contigs(by_contig,settings)
        table(out/'proposals.tsv',changes,PROPOSAL_FIELDS)
        table(out/'contigs.tsv',novelty,['contig','classification','viral_orfs','cellular_orfs','hallmark_orfs','sequence_close_orfs',
                                       'structural_rescue_orfs','novelty_class','priority_score'])
        manifest['interpretation'] = 'Curated-label-based candidate annotations; priority score is not a probability or taxonomic novelty proof'


