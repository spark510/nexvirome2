import json
from collections import defaultdict

def summarize(hits,labels,settings):
    accepted = [h for h in hits if h['evalue'] <= settings.get('evalue',1e-3)
                and min(h['qcov'],h['tcov']) >= settings.get('min_coverage',0.5)]
    if not accepted:
        return {'origin':'uncertain','lineage':'','family':'','hallmark':False,'confidence':0.0,'related':False,'close':False}
    best = max(h['bits'] for h in accepted)
    competing = [h for h in accepted if h['bits'] >= best*settings.get('competitive_fraction',0.95)]
    annotations = [labels.get(h['target'],{}) for h in competing]
    origins = {a.get('origin','uncertain') for a in annotations}
    origin = next(iter(origins)) if len(origins)==1 and origins <= {'viral','cellular'} else 'uncertain'
    lineages = {a.get('lineage','') for a in annotations}
    families = {a.get('family','') for a in annotations}
    return {'origin':origin,'lineage':next(iter(lineages)) if len(lineages)==1 else '',
            'family':next(iter(families)) if len(families)==1 else '',
            'hallmark':origin=='viral' and all(a.get('hallmark','').lower() in ('true','1','yes') for a in annotations),
            'confidence':min(h['qcov'] for h in competing) if origin!='uncertain' else 0.0,
            'related':True,'close':any(h['fident'] >= settings.get('close_identity',0.8) and h['qcov'] >= 0.8 for h in accepted)}


def annotate_orfs(orfs,evidence,labels,settings):
    profiles = []
    by_contig = defaultdict(list)
    for r in orfs:
        a,b = int(r['start']),int(r['end'])
        if not 0 <= a < b:
            raise ValueError('Invalid ORF coordinates')
        seq = summarize(evidence[0][r['protein_id']],labels,settings)
        struct = summarize(evidence[1][r['protein_id']],labels,settings)
        informative = [s for s in (seq,struct) if s['origin']!='uncertain']
        origins = {s['origin'] for s in informative}
        origin = next(iter(origins)) if len(origins)==1 else 'uncertain'
        if any(s['related'] and s['origin']=='uncertain' for s in (seq,struct)):
            origin = 'uncertain'
        lineages = {s['lineage'] for s in informative if s['lineage']}
        lineage = next(iter(lineages)) if len(lineages)==1 and origin!='uncertain' else ''
        row = {**r,'start':a,'end':b,'origin':origin,'lineage':lineage,
               'sequence_family':seq['family'],'structure_family':struct['family'],
               'sequence_related':seq['related'],'sequence_close':seq['close'],
               'structural_related':struct['related'],'structural_only':struct['related'] and not seq['related'],
               'hallmark':origin=='viral' and any(s['hallmark'] for s in informative),
               'confidence':min((s['confidence'] for s in informative),default=0.0)}
        profiles.append(row)
        by_contig[r['contig']].append(row)
    return profiles,by_contig

def summarize_contigs(by_contig,settings):
    changes,novelty = [],[]
    n = int(settings.get('flank_orfs',2))
    if n < 1:
        raise ValueError('flank_orfs must be positive')
    for contig,records in by_contig.items():
        records.sort(key=lambda r:(r['start'],r['end'],r['orf_id']))
        for i in range(n,len(records)-n+1):
            left,right = records[i-n:i],records[i:i+n]
            start,end = max(r['end'] for r in left),min(r['start'] for r in right)
            if start>end:
                continue
            for key,method in [('origin','structural_origin'),('lineage','structural_lineage')]:
                la,lb = {r[key] for r in left},{r[key] for r in right}
                if len(la)!=1 or len(lb)!=1 or la==lb or '' in la|lb or 'uncertain' in la|lb:
                    continue
                if not any(r['structural_related'] for r in left) or not any(r['structural_related'] for r in right):
                    continue
                confidence = min(r['confidence'] for r in left+right)
                if confidence < settings.get('min_confidence',0.5):
                    continue
                changes.append({'contig':contig,'start':start,'end':end,'method':method,'score':confidence,
                                'detail':json.dumps({'left':sorted(la),'right':sorted(lb),'orfs':[r['orf_id'] for r in left+right]})})
        viral = sum(r['origin']=='viral' for r in records)
        cellular = sum(r['origin']=='cellular' for r in records)
        hallmark = sum(r['hallmark'] for r in records)
        category = 'viral' if viral>=settings.get('min_viral_orfs',2) and hallmark and not cellular else 'cellular' if cellular>=2 and not viral else 'uncertain'
        close = sum(r['sequence_close'] for r in records)
        rescue = sum(r['structural_only'] and r['origin']=='viral' for r in records)
        novelty.append({'contig':contig,'classification':category,'viral_orfs':viral,'cellular_orfs':cellular,
                        'hallmark_orfs':hallmark,'sequence_close_orfs':close,'structural_rescue_orfs':rescue,
                        'novelty_class':'known_sequence' if close else 'remote_structure_candidate' if rescue else 'sequence_supported' if viral else 'unresolved',
                        'priority_score':rescue+hallmark-close})
    return changes,novelty
