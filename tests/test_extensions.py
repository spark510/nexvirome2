"""Synthetic numerical tests and mocked search wiring; no biological validation claims."""
import json
from pathlib import Path
import numpy as np
import pytest
pytest.importorskip('sklearn',reason='Install the extensions extra or extended Conda environment')
pytest.importorskip('Bio',reason='Install the extensions extra or extended Conda environment')
from nexvirome2.common import table, rows, dump, fasta, write_fasta, revcomp


def matrix_fixture(tmp_path):
    regions = [{'window_id':f'w{i}','contig':c,'start':j*100,'end':(j+1)*100}
               for i,(c,j) in enumerate([('a',0),('a',1),('b',0),('b',1),('mixed',0),('mixed',1)])]
    values = [[20,10,1,2],[40,20,2,4],[20,10,1,2],[40,20,2,4],[20,10,1,2],[1,2,20,10]]
    names = ['s1','s2','s3','s4']
    table(tmp_path/'windows.tsv',regions,['window_id','contig','start','end'])
    table(tmp_path/'coverage.tsv',[dict(window_id=r['window_id'],**dict(zip(names,v))) for r,v in zip(regions,values)],['window_id',*names])
    table(tmp_path/'libraries.tsv',[{'sample':n,'biological_sample':n} for n in names],['sample','biological_sample'])
    dump(tmp_path/'composition.json',[{'window_id':r['window_id'],'kmers':{'AAAA':1.0}} for r in regions])
    return [tmp_path/n for n in ['coverage.tsv','windows.tsv','libraries.tsv','composition.json']]


def test_coverage_pairs_and_ambiguity(tmp_path):
    from nexvirome2.features import measure_sam, windows, blocks
    target = {'a':'A'*300}
    lines = ['@SQ\tSN:a\tLN:300\n']
    for name,extra in [('good',''),('tied','\tXS:i:0')]:
        for flag,pos,mate in [(99,1,51),(147,51,1)]:
            lines.append(f'{name}\t{flag}\ta\t{pos}\t42\t100M\t=\t{mate}\t150\t*\t*\tAS:i:0{extra}\n')
    sam = tmp_path/'reads.sam'; sam.write_text(''.join(lines))
    depth,stats = measure_sam(sam,target,windows(target,100,100))
    assert depth == [1,0.5,0]
    assert stats == {'accepted_fragments':1,'excluded_fragments':1,'aligned_bases':150}
    assert list(blocks(0,'10M5D10M')) == [(0,10),(15,25)]
    with pytest.raises(ValueError,match='reference'):
        measure_sam(sam,{'a':'A'*301},[])


def test_nmf_changes_and_simple_baseline(tmp_path):
    from nexvirome2.latent import fit
    from nexvirome2.features import propose
    cov,win,lib,comp = matrix_fixture(tmp_path)
    fit({'min_stability':0.8},cov,win,lib,tmp_path/'nmf')
    candidates = rows(tmp_path/'nmf/proposals.tsv')
    assert len(candidates)==1 and candidates[0]['contig']=='mixed'
    assert (candidates[0]['start'],candidates[0]['end'])==('0','200')
    assert json.loads((tmp_path/'nmf/diagnostics.json').read_text())['selected_rank']==2
    propose({},cov,win,lib,comp,tmp_path/'simple')
    assert [(r['contig'],r['method']) for r in rows(tmp_path/'simple/proposals.tsv')]==[('mixed','coverage')]


def test_nmf_identifiability_and_replicates(tmp_path):
    from nexvirome2.latent import factorize, aggregate
    with pytest.raises(ValueError,match='unidentifiable'):
        factorize(np.ones((4,4)),{})
    with pytest.raises(ValueError,match='nonzero coverage'):
        factorize(np.array([[1,2,0],[2,1,0]]),{})
    lib = tmp_path/'libs.tsv'
    table(lib,[{'sample':f's{i}','biological_sample':'same'} for i in range(4)],['sample','biological_sample'])
    with pytest.raises(ValueError,match='three biological'):
        aggregate(np.ones((2,4)),[f's{i}' for i in range(4)],lib)


def test_orf_coordinates_dedup_and_stops(tmp_path):
    from nexvirome2.orfs import predict, six_frame
    from Bio.Seq import Seq
    dna = 'ATGAAACCCGGG'
    path = tmp_path/'c.fa'; write_fasta(path,{'a':dna,'b':revcomp(dna)})
    imported = tmp_path/'orfs.tsv'
    table(imported,[{'contig':c,'start':0,'end':12,'strand':s} for c,s in [('a','+'),('b','-')]],['contig','start','end','strand'])
    predict({},path,tmp_path/'orfs',imported)
    assert len(fasta(tmp_path/'orfs/proteins.faa'))==1
    assert len(rows(tmp_path/'orfs/orfs.tsv'))==2
    for r in six_frame('TAA'+revcomp(dna)+'TAA',minimum=2):
        fragment = ('TAA'+revcomp(dna)+'TAA')[r['start']:r['end']]
        assert str(Seq(fragment if r['strand']=='+' else revcomp(fragment)).translate())==r['protein']


def test_search_cache_keys_and_commands(tmp_path,monkeypatch):
    import nexvirome2.searches as search
    proteins = tmp_path/'p.faa'; write_fasta(proteins,{'p':'MKKL'})
    database = tmp_path/'db'; database.write_text('database-v1')
    model = tmp_path/'weights'; model.write_text('weights')
    calls = []
    def version(tool,flag,out,manifest):
        (out/f'{tool}_version.stdout.log').write_text('test-version')
    def command(args,out,manifest,name):
        calls.append(args)
        Path(args[4]).write_text('p\ttarget\t0.4\t0.9\t0.8\t1e-10\t80\n')
    monkeypatch.setattr(search,'version',version)
    monkeypatch.setattr(search,'command',command)
    cache = tmp_path/'cache.sqlite'
    for i in range(2):
        search.search({},proteins,database,tmp_path/f'run{i}',model=model,cache=cache)
    assert len(calls)==1 and '--prostt5-model' in calls[0]
    assert rows(tmp_path/'run0/hits.tsv')==rows(tmp_path/'run1/hits.tsv')
    database.write_text('database-v2')
    search.search({},proteins,database,tmp_path/'run2',model=model,cache=cache)
    assert len(calls)==2
    with pytest.raises(ValueError,match='exactly one'):
        search.search({},proteins,database,tmp_path/'bad')


def test_structural_lineages_and_functional_families(tmp_path):
    from nexvirome2.structural import annotation
    from nexvirome2.searches import HIT_FIELDS
    orfs = tmp_path/'orfs.tsv'; labels = tmp_path/'labels.tsv'; hits = tmp_path/'hits.tsv'
    table(orfs,[dict(orf_id=f'o{i}',protein_id=f'p{i}',contig='a',start=i*100,end=i*100+90,strand='+') for i in range(4)],
          ['orf_id','protein_id','contig','start','end','strand'])
    meta = [dict(target=f't{i}',origin='viral',lineage='A' if i<2 else 'B',family=f'function{i}',hallmark='true') for i in range(4)]
    table(labels,meta,['target','origin','lineage','family','hallmark'])
    table(hits,[dict(query=f'p{i}',target=f't{i}',fident=0.3,qcov=0.9,tcov=0.9,evalue=1e-12,bits=100) for i in range(4)],HIT_FIELDS)
    annotation({},orfs,labels,tmp_path/'annotation',structure_hits=hits)
    proposals = rows(tmp_path/'annotation/proposals.tsv')
    assert [(r['start'],r['end'],r['method']) for r in proposals]==[('190','200','structural_lineage')]
    assert rows(tmp_path/'annotation/contigs.tsv')[0]['novelty_class']=='remote_structure_candidate'
    for r in meta: r['lineage']='A'
    table(labels,meta,['target','origin','lineage','family','hallmark'])
    annotation({},orfs,labels,tmp_path/'same_lineage',structure_hits=hits)
    assert rows(tmp_path/'same_lineage/proposals.tsv')==[]


def test_competitive_unknown_and_cellular_hits():
    from nexvirome2.structural import summarize
    hit = dict(target='viral',fident=.3,qcov=.9,tcov=.9,evalue=1e-10,bits=100)
    labels = {'viral':{'origin':'viral','lineage':'A','family':'capsid','hallmark':'true'},'cell':{'origin':'cellular'}}
    for competitor in ['unknown','cell']:
        result = summarize([hit,{**hit,'target':competitor}],labels,{})
        assert result['origin']=='uncertain' and not result['hallmark']


def test_proposal_gate_never_creates_split():
    from nexvirome2.structural import gate
    candidate = {'contig':'a','cut':100}
    proposal = dict(contig='a',start=90,end=110,method='nmf',score=.99)
    for policy in ['annotate','require_any','require_all']:
        settings = {'proposal_policy':policy,'required_methods':['nmf','structural_lineage']}
        assert gate('unresolved',candidate,[proposal],settings)[0]=='unresolved'
        assert gate('retain',candidate,[proposal],settings)[0]=='retain'
    assert gate('split',candidate,[],{'proposal_policy':'require_any'})[0]=='unresolved'
    assert gate('split',candidate,[proposal],{'proposal_policy':'require_any'})[0]=='split'


def test_bins_preserve_mosaic_as_unresolved(tmp_path):
    from nexvirome2.binning import build
    cov,win,lib,comp = matrix_fixture(tmp_path)
    build({},cov,win,lib,comp,tmp_path/'bins')
    result = {r['contig']:r for r in rows(tmp_path/'bins/bins.tsv')}
    assert result['a']['bin']==result['b']['bin']
    assert result['mixed']['status']=='unresolved'
    assert not list((tmp_path/'bins').glob('*.fasta'))


@pytest.mark.parametrize('tool,subcommand',[('genomad','end-to-end'),('checkv','end_to_end')])
def test_quality_wrapper_preserves_inputs(tmp_path,monkeypatch,tool,subcommand):
    import nexvirome2.quality as quality
    contigs=tmp_path/'contigs.fa'; write_fasta(contigs,{'a':'ACGT'})
    database=tmp_path/'db'; database.mkdir(); (database/'version').write_text('1')
    calls=[]
    monkeypatch.setattr(quality,'version',lambda *a:None)
    monkeypatch.setattr(quality,'command',lambda args,*a:calls.append(args))
    quality.run({'threads':2},tool,contigs,database,tmp_path/'qc')
    assert calls[0][:2]==[tool,subcommand]
    assert contigs in calls[0] and database in calls[0]
    assert fasta(contigs)=={'a':'ACGT'}


def test_local_neighborhoods_merge_and_missing_paths(tmp_path):
    from nexvirome2.features import neighborhoods
    gfa=tmp_path/'g.gfa'
    gfa.write_text('S\t1\tAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTT\nS\t5\tACAC\n'
                   'L\t1\t+\t2\t+\t0M\nL\t5\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\n')
    contigs=tmp_path/'c.fa'; write_fasta(contigs,{'c':'AAAACCCCGGGG','missing':'ATGC'})
    paths=tmp_path/'paths'; paths.write_text('c\n1+,2+,3+\n')
    neighborhoods({'neighborhood_bp':2},contigs,gfa,paths,tmp_path/'local')
    assert rows(tmp_path/'local/regions.tsv')==[{'contig':'c','start':'2','end':'10'}]
    assert any(not r['usable'] for r in json.loads((tmp_path/'local/candidates.json').read_text()))


def test_imported_coverage_pipeline_and_target_digest(tmp_path,monkeypatch):
    import shutil
    import nexvirome2.features as features
    from nexvirome2.common import sha
    contigs=tmp_path/'target.fa'; write_fasta(contigs,{'a':'A'*300})
    sam=tmp_path/'reads.sam'
    sam.write_text('@SQ\tSN:a\tLN:300\n'
        'q\t99\ta\t1\t42\t100M\t=\t51\t150\t*\t*\tAS:i:0\n'
        'q\t147\ta\t51\t42\t100M\t=\t1\t-150\t*\t*\tAS:i:0\n')
    samples=tmp_path/'samples.tsv'
    records=[dict(sample=n,biological_sample=n,sam='reads.sam',target_sha256=sha(contigs)) for n in ['a','b','c']]
    fields=['sample','biological_sample','sam','target_sha256']
    table(samples,records,fields)
    def external(args,*a):
        assert args[:3]==['samtools','sort','-n']
        shutil.copyfile(args[-1],args[args.index('-o')+1])
    monkeypatch.setattr(features,'version',lambda *a:None)
    monkeypatch.setattr(features,'command',external)
    features.coverage({'window_bp':100,'step_bp':100},contigs,samples,tmp_path/'features')
    coverage=rows(tmp_path/'features/normalized_coverage.tsv')
    assert float(coverage[0]['a'])==pytest.approx(1e6/150)
    assert float(coverage[1]['a'])==pytest.approx(0.5e6/150)
    assert rows(tmp_path/'features/libraries.tsv')[0]['accepted_fragments']=='1'
    records[0]['target_sha256']='wrong'; table(samples,records,fields)
    with pytest.raises(ValueError,match='target_sha256'):
        features.coverage({},contigs,samples,tmp_path/'invalid')


def test_orf_selection_keeps_original_coordinates(tmp_path):
    from nexvirome2.orfs import predict
    contigs=tmp_path/'contigs.fa'; write_fasta(contigs,{'a':'ATGAAA'*20,'b':'ATGCCC'*20})
    imported=tmp_path/'orfs.tsv'
    table(imported,[dict(contig=c,start=30,end=90,strand='+') for c in ['a','b']],['contig','start','end','strand'])
    region=tmp_path/'regions.tsv'; table(region,[dict(contig='a',start=50,end=60)],['contig','start','end'])
    predict({},contigs,tmp_path/'orfs',imported,region)
    records=rows(tmp_path/'orfs/orfs.tsv')
    assert len(records)==1 and (records[0]['contig'],records[0]['start'],records[0]['end'])==('a','30','90')


def test_coordinate_tools_and_quality_download(tmp_path,monkeypatch):
    import nexvirome2.searches as search
    import nexvirome2.quality as quality
    structures=tmp_path/'structures'; structures.mkdir()
    for n in ['a','b']: (structures/(n+'.pdb')).write_text('fixture: not real atomic coordinates')
    calls=[]
    def external(args,out,manifest,name,stdout=None):
        calls.append(args)
        if stdout: Path(stdout).write_text('mock motif output\n')
    monkeypatch.setattr(search,'version',lambda *a:None)
    monkeypatch.setattr(search,'command',external)
    with pytest.raises(ValueError,match='motif required'):
        search.motif({},structures,tmp_path/'bad',query=structures/'a.pdb')
    search.motif({},structures,tmp_path/'motif',query=structures/'a.pdb',residues='A1,A2')
    assert calls[0][:2]==['folddisco','index'] and calls[1][:2]==['folddisco','query']
    assert 'A1,A2' in calls[1]
    search.align_structures({},structures,tmp_path/'align')
    assert calls[-1][:2]==['foldmason','easy-msa']
    monkeypatch.setattr(quality,'version',lambda *a:None)
    def download(args,*a):
        destination=Path(args[-1]); assert destination.is_dir()
        (destination/'db.version').write_text('mock-version')
    monkeypatch.setattr(quality,'command',download)
    quality.download({},'checkv',tmp_path/'download')
    assert json.loads((tmp_path/'download/database.json').read_text())['files']
