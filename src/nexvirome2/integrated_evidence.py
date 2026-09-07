"""Coordinate-preserving evidence adapters; associations are not causal claims."""
from collections import Counter, defaultdict
from .common import rows, stage, table, dump, fasta, union_length
from .masking import MaskSet


def link(settings, reference, mask_file, annotations, output):
    """Annotation TSV: accession,start,end,evidence_id,kind,label (original nt coordinates)."""
    with stage(output,settings,[reference,mask_file,annotations]) as (out,_):
        mask = MaskSet.load(mask_file,reference)
        links, by_kind = [], defaultdict(list)
        for row in rows(annotations):
            acc,a,b = row['accession'],int(row['start']),int(row['end'])
            if acc not in mask.lengths or not 0<=a<b<=mask.lengths[acc]:
                raise ValueError('Annotation coordinates do not match original reference')
            if not row['evidence_id'] or not row['kind']:
                raise ValueError('Annotation requires evidence_id and kind')
            overlaps = [(max(a,x),min(b,y)) for x,y in mask.merged(acc) if max(a,x)<min(b,y)]
            amount = union_length(overlaps)
            links.append(dict(**row,masked_bases=amount,annotation_fraction_masked=amount/(b-a)))
            by_kind[(acc,row['kind'])].extend(overlaps)
        table(out/'links.tsv',links,['accession','start','end','evidence_id','kind','label',
                                   'masked_bases','annotation_fraction_masked'])
        table(out/'summary.tsv',[dict(accession=acc,kind=kind,masked_annotated_bases=union_length(intervals))
                                 for (acc,kind),intervals in sorted(by_kind.items())],
              ['accession','kind','masked_annotated_bases'])
        dump(out/'interpretation.json',dict(annotation_overlap_is_functional_validation=False,
             note='Evidence on masked intervals is retained for detection/annotation, not used to force fine taxonomy'))


def domains(settings, reference, orfs, hits, output):
    """Normalize curated protein evidence: orf_id,aa_start,aa_end,evidence_id,kind,label.

Both amino-acid and nucleotide coordinates are 0-based half-open. Strand is
    required; imported ORF boundaries must describe complete codons. Evidence
    providers must exclude stop codons from protein domain annotations.
"""
    with stage(output,settings,[reference,orfs,hits]) as (out,_):
        sequences, index = fasta(reference), {}
        for row in rows(orfs):
            if row['orf_id'] in index:
                raise ValueError('Duplicate ORF ID')
            acc,a,b = row['contig'],int(row['start']),int(row['end'])
            if acc not in sequences or not 0<=a<b<=len(sequences[acc]) or (b-a)%3 or row['strand'] not in ('+','-'):
                raise ValueError('Invalid original-coordinate ORF')
            index[row['orf_id']] = row
        result=[]
        for row in rows(hits):
            orf=index[row['orf_id']]
            a,b = int(row['aa_start']),int(row['aa_end'])
            start,end = int(orf['start']),int(orf['end'])
            if not 0<=a<b<=(end-start)//3:
                raise ValueError('Domain outside ORF')
            left,right = (start+3*a,start+3*b) if orf['strand']=='+' else (end-3*b,end-3*a)
            result.append(dict(accession=orf['contig'],start=left,end=right,evidence_id=row['evidence_id'],
                               kind=row['kind'],label=row['label'],orf_id=row['orf_id'],strand=orf['strand']))
        table(out/'annotations.tsv',result,['accession','start','end','evidence_id','kind','label','orf_id','strand'])


def joint(settings, reference, observations, output):
    """Disjoint assay units: unit,accession,start,end,classification,assembly.

classification: ambiguous/informative/unknown; assembly: chimeric/consistent/unknown.
Unresolved calls remain outside the 2x2 denominator. Distinct overlapping
observations are rejected rather than counted as independent base evidence.
"""
    with stage(output,settings,[reference,observations]) as (out,_):
        sequences=fasta(reference)
        seen=defaultdict(list)
        counts=Counter()
        for row in rows(observations):
            acc,a,b=row['accession'],int(row['start']),int(row['end'])
            if not row['unit'] or acc not in sequences or not 0<=a<b<=len(sequences[acc]):
                raise ValueError('Invalid joint observation coordinates')
            key=(row['unit'],acc)
            if any(max(a,x)<min(b,y) for x,y in seen[key]):
                raise ValueError('Overlapping joint observations within an assay unit')
            seen[key].append((a,b))
            c,assembly=row['classification'],row['assembly']
            if c not in ('ambiguous','informative','unknown') or assembly not in ('chimeric','consistent','unknown'):
                raise ValueError('Invalid joint observation status')
            counts[(c,assembly)]+=b-a
        table(out/'joint.tsv',[dict(classification=c,assembly=a,bases=n) for (c,a),n in sorted(counts.items())],
              ['classification','assembly','bases'])
        dump(out/'summary.json',dict(known_bases=sum(n for (c,a),n in counts.items() if 'unknown' not in (c,a)),
             unresolved_bases=sum(n for (c,a),n in counts.items() if 'unknown' in (c,a)),
             causal_interpretation=False,independent_units=len({unit for unit,_ in seen})))


def discovery(settings, profiles, output, motifs=None, context=None):
    """Four separate novelty axes. No-hit and not-searched are distinct states."""
    with stage(output,settings,[p for p in (profiles,motifs,context) if p]) as (out,_):
        records=rows(profiles)
        ids={r['orf_id'] for r in records}
        if len(ids)!=len(records):
            raise ValueError('Duplicate ORF profiles')
        def evidence(path):
            result=defaultdict(list)
            if path:
                for r in rows(path):
                    if r['orf_id'] not in ids or r['status'] not in ('supported','unsupported','uncertain'):
                        raise ValueError('Invalid discovery evidence')
                    result[r['orf_id']].append(r)
            return result
        motif_map,context_map=evidence(motifs),evidence(context)
        def axis(entries):
            statuses={r['status'] for r in entries}
            if not statuses:
                return 'not_assessed'
            return next(iter(statuses)) if len(statuses)==1 else 'uncertain'
        def flag(row,key):
            value=str(row[key]).lower()
            if value not in ('true','false'):
                raise ValueError('Profile flags must be true/false')
            return value=='true'
        cards=[]
        for row in records:
            seq,fold=flag(row,'sequence_related'),flag(row,'structural_related')
            ma,ca=axis(motif_map[row['orf_id']]),axis(context_map[row['orf_id']])
            sequence_assessed=settings.get('sequence_assessed',False)
            fold_assessed=settings.get('fold_assessed',False)
            if seq:
                novelty='sequence_supported'
            elif not sequence_assessed:
                novelty='insufficient_assessment'
            elif fold:
                novelty='remote_fold_candidate'
            elif not fold_assessed:
                novelty='insufficient_assessment'
            elif ma=='supported':
                novelty='motif_only_candidate'
            elif ma=='unsupported' and ca=='supported':
                novelty='dark_context_candidate'
            else:
                novelty='unresolved'
            cards.append(dict(orf_id=row['orf_id'],contig=row['contig'],
                sequence='supported' if seq else 'unsupported' if sequence_assessed else 'not_assessed',
                fold='supported' if fold else 'unsupported' if fold_assessed else 'not_assessed',
                motif=ma,context=ca,novelty=novelty,origin=row.get('origin','uncertain')))
        table(out/'cards.tsv',cards,['orf_id','contig','sequence','fold','motif','context','novelty','origin'])
        dump(out/'summary.json',dict(counts=dict(Counter(r['novelty'] for r in cards)),
             new_family_confirmed=False,interpretation='Candidate evidence categories; no-hit is not proof of novelty'))
