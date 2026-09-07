import gzip
import pytest
from nexvirome2.common import write_fasta, table, rows, config
from nexvirome2.catalog_clustering import cluster
from nexvirome2.coassembly import prepare


def test_catalog_does_not_chain_distinct_representatives(tmp_path):
    sequences={'A':'AAAAAAAAAA','B':'CAAAAAAAAA','C':'CCAAAAAAAA'}
    write_fasta(tmp_path/'ref.fa',sequences)
    hits=[]
    for a,b in [('A','B'),('B','C')]:
        hits.append(f'{a}\t{b}\t90\t10\t1\t10\t1\t10\t50\t{sequences[a]}\t{sequences[b]}')
    (tmp_path/'hits.tsv').write_text('\n'.join(hits)+'\n')
    cluster({'min_identity':90},tmp_path/'ref.fa',tmp_path/'out',tmp_path/'hits.tsv')
    mapping={r['member']:r['representative'] for r in rows(tmp_path/'out/clusters.tsv')}
    assert mapping=={'A':'A','B':'A','C':'C'}


def test_catalog_rejects_short_partial_match_and_corrupt_sequence(tmp_path):
    write_fasta(tmp_path/'ref.fa',{'A':'ACGT'*5,'B':'ACGT'})
    line='A\tB\t100\t4\t1\t4\t1\t4\t10\tACGT\tACGT\n'
    (tmp_path/'hits.tsv').write_text(line)
    cluster({},tmp_path/'ref.fa',tmp_path/'out',tmp_path/'hits.tsv')
    assert config(tmp_path/'out/summary.json')['representatives']==2
    (tmp_path/'hits.tsv').write_text(line.replace('ACGT\tACGT','ACGT\tAAAA'))
    with pytest.raises(ValueError,match='differs from original'):
        cluster({},tmp_path/'ref.fa',tmp_path/'bad',tmp_path/'hits.tsv')


def sample_inputs(tmp_path):
    for sample in ('s1','s2'):
        for mate in (1,2):
            (tmp_path/f'{sample}_{mate}.fq').write_text(f'@q/{mate}\nACGT\n+\nIIII\n@unused/{mate}\nTGCA\n+\nIIII\n')
    table(tmp_path/'samples.tsv',[dict(sample=s,biological_sample='same_bio',library='PE150',r1=f'{s}_1.fq',r2=f'{s}_2.fq')
                                 for s in ('s1','s2')],['sample','biological_sample','library','r1','r2'])
    table(tmp_path/'members.tsv',[dict(group='g',sample=s,fragment='q') for s in ('s1','s2')],['group','sample','fragment'])


def test_targeted_pool_preserves_lineage_and_technical_grouping(tmp_path):
    sample_inputs(tmp_path)
    prepare({},tmp_path/'samples.tsv',tmp_path/'members.tsv',tmp_path/'out')
    with gzip.open(tmp_path/'out/g/reads_R1.fastq.gz','rt') as handle:
        text=handle.read()
    assert '@s1:q/1' in text and '@s2:q/1' in text and 'unused' not in text
    assert rows(tmp_path/'out/pools.tsv')[0]['biological_samples']=='1'
    assert len(rows(tmp_path/'out/lineage.tsv'))==2


def test_unknown_selected_fragment_does_not_complete_pool(tmp_path):
    sample_inputs(tmp_path)
    table(tmp_path/'members.tsv',[dict(group='g',sample='s1',fragment='missing')],['group','sample','fragment'])
    with pytest.raises(ValueError,match='absent from reads'):
        prepare({},tmp_path/'samples.tsv',tmp_path/'members.tsv',tmp_path/'out')
    assert config(tmp_path/'out/manifest.json')['status']=='failed'
