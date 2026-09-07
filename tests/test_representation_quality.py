from pathlib import Path
import pytest
from nexvirome2.common import write_fasta,table,rows
from nexvirome2.search.representations import prepare
from nexvirome2.quality_results import genomad


def test_representation_cache_reuse_and_corruption(tmp_path):
    model=tmp_path/'model'
    model.write_text('weights')
    calls=[]
    def runner(args,out,manifest,name):
        calls.append(args)
        for prefix in (str(args[3]),str(args[3])+'_ss'):
            for suffix in ('','.dbtype','.index'):
                Path(prefix+suffix).write_text('mock_db')
    first,identity=prepare(tmp_path/'cache',{'p':'AAAA'},model,[model],'v1',{},runner)
    second,reused=prepare(tmp_path/'cache',{'p':'AAAA'},model,[model],'v1',{},runner)
    assert first==second and len(calls)==1 and reused['cache_hit']
    Path(str(first)+'.dbtype').write_text('corrupt')
    with pytest.raises(ValueError,match='checksum'):
        prepare(tmp_path/'cache',{'p':'AAAA'},model,[model],'v1',{},runner)


def test_foldseek_reuses_prediction_for_a_second_target_database(tmp_path):
    from nexvirome2.search.service import SearchService
    from nexvirome2.search.backends import FoldseekBackend
    from nexvirome2.common import config
    write_fasta(tmp_path/'protein.fa',{'p':'AAAA'})
    model=tmp_path/'model'
    model.write_text('weights')
    for name in ('db1','db2'):
        (tmp_path/name).write_text(name)
    class Runner:
        def __init__(self):
            self.predictions=0
        def version(self,tool,flag,out,manifest):
            (out/'foldseek_version.stdout.log').write_text('v1')
        def __call__(self,args,out,manifest,name):
            if args[1]=='createdb':
                self.predictions+=1
                for prefix in (str(args[3]),str(args[3])+'_ss'):
                    for suffix in ('','.dbtype','.index'):
                        Path(prefix+suffix).write_text('mock_db')
            elif args[1]=='convertalis':
                Path(args[5]).write_text('p\tt\t0.5\t1\t1\t0.00001\t50\n')
    runner=Runner()
    settings={'representation_cache':str(tmp_path/'representations')}
    for name in ('db1','db2'):
        SearchService(FoldseekBackend(),runner).run(settings,tmp_path/'protein.fa',tmp_path/name,
                  tmp_path/(name+'_search'),model=model,cache=tmp_path/'hits.sqlite')
    assert runner.predictions==1
    assert config(tmp_path/'db2_search/manifest.json')['representation_cache']['cache_hit']


def test_genomad_provirus_coordinates_and_missing_not_negative(tmp_path):
    write_fasta(tmp_path/'ref.fa',{'host':'ACGT'*25,'other':'ACGT'})
    record=dict(seq_name='host|provirus_11_30',length=20,coordinates='11-30',virus_score=.95,
                topology='Provirus',taxonomy='Viruses;')
    table(tmp_path/'summary.tsv',[record],list(record))
    genomad({},tmp_path/'ref.fa',tmp_path/'summary.tsv',tmp_path/'qc')
    result=rows(tmp_path/'qc/quality.tsv')
    assert (result[0]['start'],result[0]['end'])==('10','30')
    assert result[0]['provirus_is_chimera']=='False'
    assert result[1]['status']=='not_reported'
    assert result[1]['virus_score']==''
