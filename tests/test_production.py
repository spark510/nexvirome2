import pytest
from nexvirome2.common import dump,table,config,write_fasta
from nexvirome2.masking import MaskSet
from nexvirome2.decisions import freeze,select
from nexvirome2.production import execute


def test_approved_frozen_mask_executes_and_tampering_is_rejected(tmp_path):
    write_fasta(tmp_path/'ref.fa',{'A':'ACGTACGT'})
    MaskSet.create(tmp_path/'ref.fa',[dict(accession='A',start=2,end=4)],
                   dict(version='v1',rank='species',conditions={})).save(tmp_path/'mask.json')
    parameters={'modules':{'mask':{'implementation':{'command':['mask','export'],'settings':{'mode':'hard'}}}}}
    table(tmp_path/'dev.tsv',[dict(unit='dev')],['unit'])
    freeze({'parameters':parameters},tmp_path/'dev.tsv',tmp_path/'frozen')
    evidence=[dict(module='mask',unit=f'u{i}',baseline_error=.2,error=.01,baseline_recovery=.9,
                   recovery=.9,resource_ratio=1) for i in range(4)]
    table(tmp_path/'e.tsv',evidence,list(evidence[0]))
    select({'parameters':parameters,'evidence_kind':'natural_heldout'},tmp_path/'e.tsv',
           tmp_path/'frozen/frozen.json',tmp_path/'selected')
    dump(tmp_path/'recipes.json',{'stages':[dict(module='mask',arguments={
         'reference':str(tmp_path/'ref.fa'),'mask':str(tmp_path/'mask.json')})]})
    execute({},tmp_path/'selected/production.json',tmp_path/'selected/decisions.json',tmp_path/'recipes.json',tmp_path/'run')
    assert 'ACNNACGT' in (tmp_path/'run/mask/reference.fasta').read_text()
    assert config(tmp_path/'run/execution.json')['truth_used'] is False
    for key,value in [('conf',str(tmp_path/'ref.fa')),('reference',['--config',str(tmp_path/'ref.fa')])]:
        recipe={'stages':[dict(module='mask',arguments={
            'reference':str(tmp_path/'ref.fa'),'mask':str(tmp_path/'mask.json'),key:value})]}
        dump(tmp_path/'override_recipes.json',recipe)
        with pytest.raises(ValueError,match='production argument|option tokens'):
            execute({},tmp_path/'selected/production.json',tmp_path/'selected/decisions.json',
                    tmp_path/'override_recipes.json',tmp_path/('override_'+key))
        assert not (tmp_path/('override_'+key)/'mask/manifest.json').exists()
    policy=config(tmp_path/'selected/production.json')
    policy['implementations']['mask']['settings']['mode']='unmasked'
    dump(tmp_path/'changed.json',policy)
    with pytest.raises(ValueError,match='Implementation differs'):
        execute({},tmp_path/'changed.json',tmp_path/'selected/decisions.json',tmp_path/'recipes.json',tmp_path/'bad')
