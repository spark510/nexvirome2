"""Expected-correct behavior for the five independently reproduced audit bugs."""
import pytest
from nexvirome2.common import table, dump, rows, config
from nexvirome2.io.sam import measure_sam
from nexvirome2.correction.support import pair_evidence
from nexvirome2.correction.policy import decision
from nexvirome2.multiview import fit, propose


def pair(name, target='c', position=51, mapq=42, secondary=False, cigar='50M', score=0):
    a,b=99+(256 if secondary else 0),147+(256 if secondary else 0)
    tag=f'\tAS:i:{score}' if score is not None else ''
    return (f'{name}\t{a}\t{target}\t{position}\t{mapq}\t{cigar}\tc\t{position+200}\t250\t*\t*{tag}\n'
            f'{name}\t{b}\t{target}\t{position+200}\t{mapq}\t50M\tc\t{position}\t-250\t*\t*{tag}\n')


@pytest.mark.parametrize('mapq,accepted',[(255,0),(19,0),(20,1),(42,1)])
def test_unavailable_mapq_is_not_unique_coverage(tmp_path,mapq,accepted):
    sam=tmp_path/'reads.sam'
    sam.write_text('@SQ\tSN:c\tLN:500\n'+pair('q',mapq=mapq))
    depth,stats=measure_sam(sam,{'c':'A'*500},[dict(contig='c',start=0,end=500)])
    assert stats['accepted_fragments']==accepted
    assert depth==[.2*accepted]


@pytest.mark.parametrize('score,accepted',[(0,0),(1,0),(None,0),(-1,1)])
def test_explicit_secondary_competition_without_xs(tmp_path,score,accepted):
    sam=tmp_path/'reads.sam'
    sam.write_text('@SQ\tSN:c\tLN:500\n@SQ\tSN:d\tLN:500\n'+pair('q')+pair('q',target='d',secondary=True,score=score))
    _,stats=measure_sam(sam,{'c':'A'*500,'d':'A'*500},[dict(contig='c',start=0,end=500)])
    assert stats['accepted_fragments']==accepted


@pytest.mark.parametrize('cigar,expected',[('10M20D140M',0),('10M20N140M',0),('20M20D130M',3)])
def test_only_aligned_anchor_bases_can_trigger_split(tmp_path,cigar,expected):
    sam=tmp_path/'reads.sam'
    sam.write_text(''.join(pair(f'q{i}',target='alt',position=76+i,cigar=cigar).replace('\tc\t','\talt\t') for i in range(3)))
    support,_=pair_evidence(sam,{'alt':(100,200),'current':(100,200)})
    assert support['alt']==expected
    assert decision(0,support['alt'])==('split' if expected else 'unresolved')


@pytest.mark.parametrize('command',[
    ['mask','export','--reference','r','--mask','m'],
    ['assemble','--assembler','megahit','--r1','a','--r2','b'],
    ['multiview','fit','--features','f'],
])
def test_no_subparser_accepts_config_abbreviation(command):
    from nexvirome2.cli import parser
    with pytest.raises(SystemExit):
        parser().parse_args([*command,'--config','frozen.yaml','--output','out','--conf','override.yaml'])
    assert parser().parse_args([*command,'--config','frozen.yaml','--output','out']).config=='frozen.yaml'


def test_zero_entity_survives_fit_and_proposal_without_imputation(tmp_path):
    data=[]
    for entity,values in [('zero',[0,0]),('a',[10,1]),('b',[1,10]),('c',[9,1]),('d',[1,9])]:
        for feature,value in zip(('f1','f2'),values):
            data.append(dict(entity=entity,block='coverage',feature=feature,value=value))
    table(tmp_path/'features.tsv',data,['entity','block','feature','value'])
    fit({'rank':2,'seeds':[17,29]},tmp_path/'features.tsv',tmp_path/'fit')
    zero=next(r for r in rows(tmp_path/'fit/loadings.tsv') if r['entity']=='zero')
    assert zero['status']=='inactive'
    assert float(zero['component_0'])==float(zero['component_1'])==0
    assert config(tmp_path/'fit/diagnostics.json')['inactive_entities']==['zero']
    table(tmp_path/'coordinates.tsv',[dict(entity=e,contig='contig',start=i*100,end=(i+1)*100)
          for i,e in enumerate(('a','zero','b','c','d'))],['entity','contig','start','end'])
    propose({},tmp_path/'fit',tmp_path/'coordinates.tsv',tmp_path/'propose')
    assert config(tmp_path/'propose/manifest.json')['status']=='complete'
    assert all(int(r['start']) not in (100,200) for r in rows(tmp_path/'propose/proposals.tsv'))


def test_inactive_window_cannot_be_bridged(tmp_path):
    model=tmp_path/'model'
    model.mkdir()
    table(model/'loadings.tsv',[
        dict(entity='a',status='active',component_0=1,component_1=0),
        dict(entity='zero',status='inactive',component_0=0,component_1=0),
        dict(entity='b',status='active',component_0=0,component_1=1)],
        ['entity','status','component_0','component_1'])
    dump(model/'diagnostics.json',dict(inactive_entities=['zero'],proposal_eligible=True,stability=1))
    table(tmp_path/'coords.tsv',[dict(entity=e,contig='c',start=i*100,end=(i+1)*100)
         for i,e in enumerate(('a','zero','b'))],['entity','contig','start','end'])
    propose({},model,tmp_path/'coords.tsv',tmp_path/'out')
    assert rows(tmp_path/'out/proposals.tsv')==[]
