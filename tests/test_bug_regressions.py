"""Reproductions from the post-implementation bug audit."""
import json
import pytest
from nexvirome2.common import stage, rows, table, dump, write_fasta


def test_scheduler_created_empty_output_is_accepted_but_completed_run_is_not(tmp_path):
    out=tmp_path/'precreated'
    out.mkdir()
    with stage(out,{}) as (folder,manifest):
        assert folder==out
    assert json.loads((out/'manifest.json').read_text())['status']=='complete'
    with pytest.raises(FileExistsError):
        with stage(out,{}):
            pytest.fail('Completed stage must not be reused')


def test_existing_untracked_output_is_not_overwritten(tmp_path):
    out=tmp_path/'existing'
    out.mkdir()
    (out/'user.txt').write_text('keep')
    with pytest.raises(FileExistsError):
        with stage(out,{}):
            pytest.fail('Existing output must not be overwritten')
    assert (out/'user.txt').read_text()=='keep'


def test_missing_input_still_records_failed_stage(tmp_path):
    out = tmp_path/'failed'
    with pytest.raises(FileNotFoundError):
        with stage(out,{},[tmp_path/'missing.fa']):
            pytest.fail('Missing input must prevent stage execution')
    manifest = json.loads((out/'manifest.json').read_text())
    assert manifest['status']=='failed'
    assert 'FileNotFoundError' in manifest['error']


def test_duplicate_tsv_columns_are_not_silently_lost(tmp_path):
    path=tmp_path/'coverage.tsv'
    path.write_text('window_id\tsample_a\tsample_a\nw0\t5\t90\n')
    with pytest.raises(ValueError,match='Duplicate'):
        rows(path)


@pytest.mark.parametrize('line',['w0\t5','w0\t5\t6\t7'])
def test_tsv_wrong_column_counts_are_rejected(tmp_path,line):
    path=tmp_path/'coverage.tsv'; path.write_text('window_id\ta\tb\n'+line+'\n')
    with pytest.raises(ValueError,match='column count mismatch'):
        rows(path)


def test_imported_orf_rejects_two_terminal_stops(tmp_path):
    pytest.importorskip('Bio')
    from nexvirome2.orfs import predict
    path=tmp_path/'contigs.fa'; write_fasta(path,{'a':'ATGAAATAATAA'})
    calls=tmp_path/'orfs.tsv'
    table(calls,[dict(contig='a',start=0,end=12,strand='+')],['contig','start','end','strand'])
    with pytest.raises(ValueError,match='internal stop'):
        predict({},path,tmp_path/'output',calls)


def test_single_terminal_stop_is_still_valid(tmp_path):
    pytest.importorskip('Bio')
    from nexvirome2.orfs import predict
    from nexvirome2.common import fasta
    path=tmp_path/'contigs.fa'; write_fasta(path,{'a':'ATGAAATAA'})
    calls=tmp_path/'orfs.tsv'
    table(calls,[dict(contig='a',start=0,end=9,strand='+')],['contig','start','end','strand'])
    predict({},path,tmp_path/'output',calls)
    assert list(fasta(tmp_path/'output/proteins.faa').values())==['MK']
    assert rows(tmp_path/'output/orfs.tsv')[0]['end']=='9'


def test_missing_structural_annotation_is_not_negative_evidence(tmp_path):
    from nexvirome2.binning import build
    regions=[dict(window_id=c,contig=c,start=0,end=100) for c in ['a','b']]
    win=tmp_path/'windows.tsv'; table(win,regions,['window_id','contig','start','end'])
    cov=tmp_path/'coverage.tsv'
    table(cov,[dict(window_id=c,s1=1,s2=2,s3=3) for c in ['a','b']],['window_id','s1','s2','s3'])
    lib=tmp_path/'libraries.tsv'
    table(lib,[dict(sample=n,biological_sample=n) for n in ['s1','s2','s3']],['sample','biological_sample'])
    comp=tmp_path/'composition.json'; dump(comp,[dict(window_id=c,kmers={'AAAA':1}) for c in ['a','b']])
    profiles=tmp_path/'profiles.tsv'; table(profiles,[dict(contig='a',structure_family='capsid')],['contig','structure_family'])
    build({},cov,win,lib,comp,tmp_path/'bins',profiles)
    bins=rows(tmp_path/'bins/bins.tsv')
    assert bins[0]['bin']==bins[1]['bin']
    assert rows(tmp_path/'bins/edges.tsv')[0]['structural_jaccard']==''


@pytest.mark.parametrize('reverse',[False,True])
@pytest.mark.parametrize('query_gap',[False,True])
def test_alignment_gap_is_not_recovered_sequence(tmp_path,reverse,query_gap):
    from nexvirome2.evaluation import evaluate
    qseq='A'*150+('-' if query_gap else 'A')+'A'*150
    sseq='A'*150+('A' if query_gap else '-')+'A'*150
    contigs=tmp_path/'c.fa'; write_fasta(contigs,{'c':qseq.replace('-','')})
    source_length=len(sseq.replace('-',''))
    sources=tmp_path/'s.fa'; write_fasta(sources,{'s':('T' if reverse else 'A')*source_length})
    alignment=tmp_path/'hits.tsv'
    alignment.write_text('\t'.join(map(str,['c','s',99.667,301,1,len(qseq.replace('-','')),
        source_length if reverse else 1,1 if reverse else source_length,600,qseq,sseq]))+'\n')
    evaluate({},contigs,sources,tmp_path/'evaluation',alignment)
    metrics=json.loads((tmp_path/'evaluation/metrics.json').read_text())['0']
    assert metrics['per_source']['s']['recovered_bp']==300
    assert metrics['per_source']['s']['longest_consistent_block_bp']==150
