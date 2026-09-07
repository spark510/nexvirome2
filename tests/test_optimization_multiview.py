import pytest
from nexvirome2.common import write_fasta,table,dump,config,sha,rows
from nexvirome2.masking import MaskSet
from nexvirome2.mask_optimization import optimize
from nexvirome2.multiview import FeatureMatrix,normalize,factorize,fit,propose


def test_optimizer_rejects_recall_damage_and_requires_paired_queries(tmp_path):
    ref=tmp_path/'ref.fa'
    write_fasta(ref,{'A':'ACGT'*50})
    for name in ('good','harm'):
        MaskSet.create(ref,[dict(accession='A',start=1,end=10)],
                       dict(version='v1',rank='species',conditions={})).save(tmp_path/f'{name}.json')
    dump(tmp_path/'candidates.json',dict(conditions={'pe150':dict(rank='species',read_length=150,
         library='genome',database_sha256=sha(ref))},policies={'unmasked':None,'good':'good.json','harm':'harm.json'}))
    data=[]
    for i in range(4):
        for name,wrong,correct,breadth in [('unmasked',.1,.8,.9),('good',.04,.8,.9),('harm',0,.4,.5)]:
            data.append(dict(condition='pe150',policy=name,unit=f'u{i}',queries_sha256='a'*64,
                             wrong=wrong,correct=correct,unclassified=.1,retained_breadth=breadth))
    table(tmp_path/'m.tsv',data,list(data[0]))
    settings={'partition':'development'}
    optimize(settings,ref,tmp_path/'candidates.json',tmp_path/'m.tsv',tmp_path/'out')
    selected=config(tmp_path/'out/selection.json')
    assert selected['selection']['pe150']['policy']=='good'
    assert selected['production_approved'] is False
    with pytest.raises(ValueError,match='development'):
        optimize({'partition':'heldout'},ref,tmp_path/'candidates.json',tmp_path/'m.tsv',tmp_path/'holdout')
    data[1]['queries_sha256']='b'*64
    table(tmp_path/'m.tsv',data,list(data[0]))
    with pytest.raises(ValueError,match='identical paired'):
        optimize(settings,ref,tmp_path/'candidates.json',tmp_path/'m.tsv',tmp_path/'bad')


def test_masked_nmf_does_not_fit_missing_values():
    import numpy as np
    x=np.array([[1.,0.,.2],[.1,1.,.3],[.8,.1,.4]])
    weights=np.ones_like(x)
    weights[0,1]=0
    alternate=x.copy()
    alternate[0,1]=999
    a,b,diag=factorize(x,weights,2,17,max_iter=200)
    c,d,other=factorize(alternate,weights,2,17,max_iter=200)
    assert np.allclose(a,c) and np.allclose(b,d)
    assert diag['loss']==pytest.approx(other['loss'])


def test_multiview_fit_and_missing_feature_weight(tmp_path):
    data=[]
    for entity,values in [('a',[10,1,1,0]),('b',[9,1,1,0]),('c',[1,10,0,1]),('d',[1,9,0,1])]:
        for (block,feature),value in zip([('coverage','s1'),('coverage','s2'),('domain','d1'),('domain','d2')],values):
            data.append(dict(entity=entity,block=block,feature=feature,value=value))
    data[0]['value']=''
    table(tmp_path/'features.tsv',data,list(data[0]))
    matrix=FeatureMatrix.load(tmp_path/'features.tsv')
    x,weights,scales=normalize(matrix,{})
    assert weights[0,0]==0
    assert weights[0,1]==1
    fit({'rank':2,'max_iter':2000},tmp_path/'features.tsv',tmp_path/'fit')
    assert config(tmp_path/'fit/diagnostics.json')['missing_entries']==1
    assert len(rows(tmp_path/'fit/loadings.tsv'))==4


def test_multiview_proposals_need_stability_and_contiguous_coordinates(tmp_path):
    model=tmp_path/'model'
    model.mkdir()
    table(model/'loadings.tsv',[dict(entity='a',component_0=.95,component_1=.05),
                               dict(entity='b',component_0=.05,component_1=.95)],
          ['entity','component_0','component_1'])
    table(tmp_path/'coords.tsv',[dict(entity='a',contig='c',start=0,end=100),
                                dict(entity='b',contig='c',start=100,end=200)],
          ['entity','contig','start','end'])
    for stable in (False,True):
        dump(model/'diagnostics.json',{'proposal_eligible':stable,'stability':.95})
        propose({},model,tmp_path/'coords.tsv',tmp_path/str(stable))
        result=rows(tmp_path/str(stable)/'proposals.tsv')
        assert len(result)==int(stable)
        if stable:
            assert result[0]['start']==result[0]['end']=='100'
