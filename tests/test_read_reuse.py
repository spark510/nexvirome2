import gzip
import pytest
from nexvirome2.common import table,dump,config,sha,rows,write_fasta
from nexvirome2.read_preparation import prepare
from nexvirome2.classification import classify
from nexvirome2.masking import MaskSet
from nexvirome2.conserved_reuse import link


@pytest.mark.parametrize('header,plus',[('@','+'),('@q/1','+other')])
def test_invalid_fastq_identifier_is_explicit_error(tmp_path,header,plus):
    from nexvirome2.read_preparation import records
    (tmp_path/'bad.fq').write_text(f'{header}\nACGT\n{plus}\nIIII\n')
    with pytest.raises(ValueError,match='identifier'):
        list(records(tmp_path/'bad.fq'))


def test_read_qc_and_optional_host_subtraction_preserve_every_pair(tmp_path):
    for mate in (1,2):
        (tmp_path/f'r{mate}.fq').write_text(''.join(f'@{name}/{mate}\nACGTACGT\n+\n{quality}\n'
                  for name,quality in [('host','IIIIIIII'),('ambiguous','IIIIIIII'),('poor','!!!!!!!!')]))
    table(tmp_path/'calls.tsv',[dict(fragment='host',status='host_supported'),dict(fragment='ambiguous',status='ambiguous')],
          ['fragment','status'])
    dump(tmp_path/'metadata.json',dict(r1_sha256=sha(tmp_path/'r1.fq'),r2_sha256=sha(tmp_path/'r2.fq'),
         method='competitive_host_viral_alignment',reference_sha256='a'*64))
    prepare(dict(min_length=4,remove_host=True),tmp_path/'r1.fq',tmp_path/'r2.fq',tmp_path/'out',tmp_path/'calls.tsv',tmp_path/'metadata.json')
    assert config(tmp_path/'out/summary.json')['counts']==dict(retained=1,host_removed=1,quality_filtered=1)
    for category,name in [('retained','ambiguous'),('host_removed','host'),('quality_filtered','poor')]:
        with gzip.open(tmp_path/f'out/{category}_R1.fastq.gz','rt') as handle:
            assert handle.read().startswith(f'@{name}/1')


def test_read_pair_mismatch_fails_without_completed_manifest(tmp_path):
    (tmp_path/'r1.fq').write_text('@a/1\nACGT\n+\nIIII\n')
    (tmp_path/'r2.fq').write_text('@b/2\nACGT\n+\nIIII\n')
    with pytest.raises(ValueError,match='out of order'):
        prepare({'min_length':1},tmp_path/'r1.fq',tmp_path/'r2.fq',tmp_path/'out')
    assert config(tmp_path/'out/manifest.json')['status']=='failed'


def test_conserved_support_is_reused_without_fine_assignment(tmp_path):
    ref,query,tax=tmp_path/'r.fa',tmp_path/'q.fa',tmp_path/'tax.tsv'
    write_fasta(ref,{'A':'ACGTACGT','B':'ACGTACGT'})
    write_fasta(query,{'q':'ACGTACGT'})
    table(tax,[dict(accession='A',species='1'),dict(accession='B',species='2')],['accession','species'])
    MaskSet.create(ref,[dict(accession=a,start=0,end=8) for a in ('A','B')],
                   dict(version='v1',rank='species',conditions={})).save(tmp_path/'mask.json')
    settings={'rank':'species','k':3}
    classify(settings,query,ref,tax,tmp_path/'masked',tmp_path/'mask.json')
    classify(settings,query,ref,tax,tmp_path/'full')
    link({},tmp_path/'masked',tmp_path/'full',tmp_path/'reuse')
    assert rows(tmp_path/'reuse/predictions.tsv')[0]['status']=='unclassified'
    assert {r['target'] for r in rows(tmp_path/'reuse/reference_support.tsv')}=={'A','B'}
    assert config(tmp_path/'reuse/summary.json')['conserved_only_queries']==1
