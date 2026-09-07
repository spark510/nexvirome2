"""Behavior contracts for typed values, composed services and execution boundaries."""
from dataclasses import FrozenInstanceError, fields
import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from nexvirome2.domain import SequenceInterval, CutInterval, BreakpointProposal, PairSupport, CorrectionDecision
from nexvirome2.domain.settings import CorrectionSettings
from nexvirome2.common import write_fasta, fasta


def test_coordinate_semantics_and_immutable_proposals():
    with pytest.raises(ValueError):
        SequenceInterval(10,10)
    assert CutInterval(10,10).contains(10)
    assert not SequenceInterval(0,10).overlaps(SequenceInterval(10,20))
    proposal=BreakpointProposal('c',10,10,'nmf',.9)
    assert proposal.matches('c',10,.5)
    assert not proposal.matches('other',10,.5)
    with pytest.raises(FrozenInstanceError):
        proposal.score=.1
    with pytest.raises(ValueError):
        BreakpointProposal('c',0,1,'nmf',float('nan'))


def test_policy_results_enforce_read_evidence_boundary():
    from nexvirome2.correction.policy import EvidencePolicy
    policy=EvidencePolicy(CorrectionSettings(proposal_policy='require_any'))
    candidate={'contig':'c','cut':10}
    proposal=BreakpointProposal('c',5,15,'nmf',.9)
    assert policy.decide(PairSupport(0,3),candidate,[proposal]).decision=='split'
    assert policy.decide(PairSupport(0,3),candidate,[]).decision=='unresolved'
    assert policy.decide(PairSupport(0,0),candidate,[proposal]).decision=='unresolved'
    with pytest.raises(ValueError,match='read-unsupported'):
        CorrectionDecision('unresolved','split')


def test_run_context_and_injected_executor_preserve_failure(tmp_path):
    from nexvirome2.runtime import RunContext, CommandRunner
    def execute(argv,stdout,stderr,check):
        stdout.write('captured output\n')
        if argv[0]=='/usr/bin/time':
            Path(argv[argv.index('-o')+1]).write_text('Maximum resident set size (kbytes): 12\n')
        return SimpleNamespace(returncode=7)
    output=tmp_path/'run'
    with pytest.raises(RuntimeError,match='exited 7'):
        with RunContext(output,{'seed':1}) as context:
            CommandRunner(execute).run(['fixture'],context.output,context.manifest,'external')
    manifest=json.loads((output/'manifest.json').read_text())
    assert manifest['status']=='failed'
    assert manifest['commands'][0]['returncode']==7
    assert (output/'external.stdout.log').read_text()=='captured output\n'


def test_cache_no_hit_rollback_and_connection_lifetime(tmp_path):
    from nexvirome2.search.cache import SearchCache
    path=tmp_path/'cache.sqlite'
    with SearchCache(path) as cache:
        assert cache.get('unknown') is None
        cache.put('empty',[])
    assert cache.connection is None
    with pytest.raises(RuntimeError):
        with SearchCache(path) as cache:
            cache.put('failed',[{'query':'x'}])
            raise RuntimeError('abort transaction')
    with SearchCache(path) as cache:
        assert cache.get('empty')==[]
        assert cache.get('failed') is None


def test_matrix_owns_sample_alignment_and_technical_aggregation():
    import numpy as np
    from nexvirome2.coverage import CoverageMatrix, SampleSet
    regions=[dict(window_id='w',contig='c',start=0,end=10)]
    matrix=CoverageMatrix(regions,['a1','a2','b','c'],np.array([[2,4,6,8]]))
    samples=SampleSet([dict(sample=s,biological_sample=g) for s,g in [('a1','a'),('a2','a'),('b','b'),('c','c')]])
    grouped=matrix.aggregate(samples)
    assert grouped.samples==['a','b','c']
    assert grouped.values.tolist()==[[3,6,8]]
    with pytest.raises(ValueError,match='shape'):
        CoverageMatrix(regions,['a','b'],np.ones((1,3)))


def test_application_search_uses_injected_backend_and_reuses_no_hits(tmp_path):
    from nexvirome2.application import SearchRequest,run_search
    from nexvirome2.search.service import SearchService
    from nexvirome2.search.backends import MMseqsBackend
    proteins=tmp_path/'p.faa'; write_fasta(proteins,{'p':'MKKL'})
    database=tmp_path/'db'; database.write_text('fixture')
    class Runner:
        def __init__(self): self.calls=[]
        def version(self,tool,flag,out,manifest):
            (out/f'{tool}_version.stdout.log').write_text('fixture')
        def __call__(self,args,out,manifest,name):
            self.calls.append(args)
            Path(args[4]).write_text('')
    runner=Runner()
    service=SearchService(MMseqsBackend(),runner)
    request=SearchRequest(str(proteins),str(database),'mmseqs',cache=str(tmp_path/'cache.sqlite'))
    for name in ['first','second']:
        run_search({},request,tmp_path/name,service)
    assert len(runner.calls)==1
    manifest=json.loads((tmp_path/'second/manifest.json').read_text())
    assert manifest['cache_hits']==1 and manifest['searched_proteins']==0


def test_application_correction_composes_a_policy_without_truth(tmp_path):
    from nexvirome2.application import CorrectionRequest,run_correction
    from nexvirome2.correction.service import CorrectionService
    assert not {'truth','sources','breakpoints'} & {f.name for f in fields(CorrectionRequest)}
    graph=tmp_path/'g.gfa'
    segments={'1':'A'*100,'2':'C'*50,'3':'G'*100,'4':'T'*100,'5':'AC'*50}
    graph.write_text(''.join(f'S\t{k}\t{v}\n' for k,v in segments.items())+
        ''.join(f'L\t{a}\t+\t{b}\t+\t0M\n' for a,b in [('1','2'),('5','2'),('2','3'),('2','4')]))
    contigs=tmp_path/'c.fa'; write_fasta(contigs,{'c':segments['1']+segments['2']+segments['3']})
    paths=tmp_path/'paths'; paths.write_text('c\n1+,2+,3+\n')
    reads=tmp_path/'reads.fq'; reads.write_text('')
    class Runner:
        def version(self,*args): pass
        def __call__(self,args,out,manifest,name,stdout=None):
            if args[0]=='bowtie2':
                lines=[]
                for i in range(3):
                    for flag,pos,mate in [(99,11+i,161+i),(147,161+i,11+i)]:
                        lines.append(f'q{i}\t{flag}\tpath_1\t{pos}\t42\t50M\t=\t{mate}\t200\t*\t*\tAS:i:0\n')
                Path(args[args.index('-S')+1]).write_text(''.join(lines))
            if stdout: Path(stdout).write_text('')
    class DeferPolicy:
        def decide(self,support,candidate,proposals):
            assert support.alternative==3
            return CorrectionDecision('split','unresolved')
    request=CorrectionRequest(str(contigs),str(graph),str(paths),str(reads),str(reads))
    for name,policy,count in [('baseline',None,2),('deferred',DeferPolicy(),1)]:
        run_correction({},request,tmp_path/name,CorrectionService(Runner(),policy=policy))
        result=fasta(tmp_path/name/'corrected.fasta')
        assert len(result)==count
        assert ''.join(result.values())==fasta(contigs)['c']


def test_ablation_settings_do_not_mutate_base():
    from nexvirome2.application import ablation_settings
    base={'proposal_policy':'annotate','threads':1}
    combined=ablation_settings(base,'combined',8)
    assert combined['required_methods']==['nmf','structural_lineage']
    assert combined['proposal_policy']=='require_all'
    assert base=={'proposal_policy':'annotate','threads':1}
