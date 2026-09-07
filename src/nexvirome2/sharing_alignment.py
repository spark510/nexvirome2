"""Approximate nucleotide sharing from competitive BLAST alignments."""
from collections import defaultdict,Counter
import math
from .common import fasta,stage,table,sha,dump,revcomp
from .taxonomy import TaxonomyResolver
from .alignment import blast,read_hits
from .evaluation.classification import aligned_runs
from .sharing import low_complexity


def build(settings,reference,taxonomy,output,alignments=None):
    rank=settings.get('rank','species')
    identity=float(settings.get('min_identity',99))
    minimum=int(settings.get('min_alignment_bp',150))
    if not 0<=identity<=100 or minimum<1:
        raise ValueError('Invalid shared-alignment thresholds')
    with stage(output,settings,[reference,taxonomy]+([alignments] if alignments else [])) as (out,manifest):
        sequences=fasta(reference)
        resolver=TaxonomyResolver.load(taxonomy)
        taxa={acc:resolver.taxon(acc,rank) for acc in sequences}
        digest=sha(reference)
        events=defaultdict(lambda:defaultdict(Counter))
        hits=read_hits(alignments or blast(reference,reference,out,manifest,threads=settings.get('threads',1)))
        for hit in hits:
            query,target=hit['qseqid'],hit['sseqid']
            if query not in sequences or target not in sequences:
                raise ValueError('Unknown shared-alignment accession')
            if not math.isfinite(hit['pident']) or not 0<=hit['pident']<=100:
                raise ValueError('Invalid alignment identity')
            if len(hit['qseq'])!=len(hit['sseq']) or len(hit['qseq'])!=hit['length']:
                raise ValueError('Inconsistent alignment strings/length')
            if len(hit['qseq'].replace('-',''))!=abs(hit['qend']-hit['qstart'])+1 or len(hit['sseq'].replace('-',''))!=abs(hit['send']-hit['sstart'])+1:
                raise ValueError('Alignment coordinates disagree with strings')
            for acc,prefix in ((query,'q'),(target,'s')):
                first,last=hit[prefix+'start'],hit[prefix+'end']
                if not 1<=min(first,last)<=max(first,last)<=len(sequences[acc]):
                    raise ValueError('Shared alignment outside original reference')
                expected=sequences[acc][min(first,last)-1:max(first,last)]
                if first>last:
                    expected=revcomp(expected)
                if hit[prefix+'seq'].replace('-','').upper()!=expected:
                    raise ValueError('Alignment sequence differs from original reference')
            if query==target or hit['pident']<identity or hit['length']<minimum:
                continue
            for run in aligned_runs(hit):
                for acc,other,column in ((query,target,0),(target,query,1)):
                    positions=[p[column] for p in run]
                    a,b=min(positions),max(positions)+1
                    if not 0<=a<b<=len(sequences[acc]):
                        raise ValueError('Shared alignment outside original reference')
                    events[acc][a][other]+=1
                    events[acc][b][other]-=1
        def atlas():
            for acc,boundaries in sorted(events.items()):
                active=Counter()
                points=sorted(boundaries)
                for index,a in enumerate(points[:-1]):
                    active.update(boundaries[a])
                    b=points[index+1]
                    competitors=sorted({acc,*[other for other,count in active.items() if count>0]})
                    if len(competitors)<2:
                        continue
                    distinct=sorted({taxa[x] for x in competitors if taxa[x] is not None})
                    yield dict(accession=acc,start=a,end=b,rank=rank,taxa_count=len(distinct),
                         taxa=';'.join(distinct),accessions=';'.join(competitors),
                         unknown_taxonomy=str(any(taxa[x] is None for x in competitors)).lower(),
                         repeat='unknown',low_complexity=str(low_complexity(sequences[acc][a:b])).lower(),
                         evidence_id=f'blast:{acc}:{a}:{b}',reference_sha256=digest,method='competitive_blast',k='')
        table(out/'atlas.tsv',atlas(),['accession','start','end','rank','taxa_count','taxa','accessions',
              'unknown_taxonomy','repeat','low_complexity','evidence_id','reference_sha256','method','k'])
        dump(out/'summary.json',dict(method='competitive_blast',reference_sha256=digest,taxonomy_sha256=sha(taxonomy),
             min_identity=identity,min_alignment_bp=minimum,
             limitations=['BLAST word-size sensitivity remains; no exhaustive homology guarantee',
                          'Circular junction alignments require a separately projected input',
                          'Repeat multiplicity is not inferred from aggregate HSPs']))
