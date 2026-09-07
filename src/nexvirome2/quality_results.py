"""geNomad native summary adapter with explicit original provirus coordinates."""
import math
import re
from .common import fasta,rows,stage,table,dump


def genomad(settings,contigs,summary,output):
    with stage(output,settings,[contigs,summary]) as (out,_):
        sequences=fasta(contigs)
        results=[]
        reported=set()
        identifiers=set()
        for row in rows(summary):
            name=row['seq_name']
            if name in identifiers:
                raise ValueError('Duplicate geNomad sequence ID')
            identifiers.add(name)
            coordinates=row['coordinates']
            if coordinates not in ('','NA'):
                match=re.fullmatch(r'(.+)\|provirus_(\d+)_(\d+)',name)
                if not match or coordinates!=f'{match[2]}-{match[3]}':
                    raise ValueError('geNomad provirus ID/coordinates disagree')
                original,start,end=match[1],int(match[2])-1,int(match[3])
            else:
                original,start,end=name,0,len(sequences.get(name,''))
            if original not in sequences or not 0<=start<end<=len(sequences[original]) or end-start!=int(row['length']):
                raise ValueError('geNomad result does not match original contig coordinates/length')
            score=float(row['virus_score'])
            if not math.isfinite(score) or not 0<=score<=1:
                raise ValueError('Invalid geNomad virus score')
            reported.add(original)
            results.append(dict(contig=original,result_id=name,start=start,end=end,status='reported_viral_candidate',
                 virus_score=score,topology=row['topology'],taxonomy=row['taxonomy'],
                 completeness=None,provirus_is_chimera=False))
        for name,seq in sequences.items():
            if name not in reported:
                results.append(dict(contig=name,result_id='',start=0,end=len(seq),status='not_reported',
                                    virus_score=None,topology='',taxonomy='',completeness=None,provirus_is_chimera=False))
        table(out/'quality.tsv',results,['contig','result_id','start','end','status','virus_score','topology','taxonomy',
                                        'completeness','provirus_is_chimera'])
        dump(out/'interpretation.json',dict(tool='genomad',missing_result_is_negative=False,
              coordinate_system='0-based-half-open',
              note='Native viral candidate scores are not completeness estimates or cross-source chimera calls'))
