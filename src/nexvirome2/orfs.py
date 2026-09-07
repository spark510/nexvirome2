"""Coordinate-preserving candidate ORFs, with optional validated imported calls."""
import hashlib
from .common import dump, fasta, revcomp, rows, sha, stage, table, write_fasta
from .domain import SequenceInterval

ORF_FIELDS = ['orf_id','protein_id','contig','start','end','strand','frame','partial','method','genetic_code']


def six_frame(sequence, minimum=75, genetic_code=1, max_unknown=0.1):
    from Bio.Seq import Seq
    if minimum < 1:
        raise ValueError('Minimum protein length must be positive')
    for strand,dna in [('+',sequence),('-',revcomp(sequence))]:
        for frame in range(3):
            length = (len(dna)-frame)//3*3
            if length <= 0:
                continue
            peptide = str(Seq(dna[frame:frame+length]).translate(table=genetic_code))
            offset = 0
            for segment in peptide.split('*'):
                left,right = frame+offset*3,frame+(offset+len(segment))*3
                if len(segment) >= minimum and segment.count('X')/len(segment) <= max_unknown:
                    a,b = (left,right) if strand == '+' else (len(sequence)-right,len(sequence)-left)
                    yield {'start':a,'end':b,'strand':strand,'frame':frame,
                           'partial':offset == 0 or offset+len(segment) == len(peptide),'protein':segment}
                offset += len(segment)+1


def predict(settings,contigs,output,imported=None,regions=None):
    from Bio.Seq import Seq
    with stage(output,settings,[contigs]+([imported] if imported else [])+([regions] if regions else [])) as (out,manifest):
        sequences = fasta(contigs)
        code = settings.get('genetic_code',1)
        calls = []
        if imported:
            for r in rows(imported):
                name = r['contig']
                a,b = int(r['start']),int(r['end'])
                strand = r['strand']
                interval = SequenceInterval(a,b)
                if name not in sequences or interval.end > len(sequences[name]) or strand not in ('+','-') or (b-a)%3:
                    raise ValueError('Imported ORFs require valid in-frame half-open coordinates')
                dna = sequences[name][a:b]
                peptide = str(Seq(dna if strand == '+' else revcomp(dna)).translate(table=code))
                if peptide.endswith('*'):
                    peptide = peptide[:-1]
                if '*' in peptide or not peptide:
                    raise ValueError('Imported ORF contains internal stop or is empty')
                calls.append({'contig':name,'start':a,'end':b,'strand':strand,'frame':-1,'partial':r.get('partial','unknown'),
                              'protein':peptide,'method':'imported_coordinates'})
        else:
            for name,seq in sequences.items():
                calls.extend({'contig':name,**r,'method':'six_frame_stop_to_stop'} for r in six_frame(
                    seq,settings.get('min_aa',75),code,settings.get('max_unknown_fraction',0.1)))
        if regions:
            selected = []
            for r in rows(regions):
                a,b = int(r['start']),int(r['end'])
                if r['contig'] not in sequences or not 0 <= a < b <= len(sequences[r['contig']]):
                    raise ValueError('Invalid ORF selection region')
                selected.append((r['contig'],a,b))
            calls = [r for r in calls if any(r['contig']==n and r['start']<b and r['end']>a for n,a,b in selected)]
        proteins,records = {},[]
        for i,r in enumerate(calls):
            protein = r.pop('protein')
            identity = 'p_'+hashlib.sha256(protein.encode()).hexdigest()
            proteins[identity] = protein
            records.append({'orf_id':f'orf_{i:09d}','protein_id':identity,**r,'genetic_code':code})
        table(out/'orfs.tsv',records,ORF_FIELDS)
        write_fasta(out/'proteins.faa',dict(sorted(proteins.items())))
        dump(out/'summary.json',{'orfs':len(records),'unique_proteins':len(proteins),
                                'status':'candidate ORFs, not validated gene annotation' if not imported else 'imported calls'})
        manifest['target_sha256'] = sha(contigs)
