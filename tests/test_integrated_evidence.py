import pytest
from nexvirome2.common import write_fasta,table,rows,config,revcomp
from nexvirome2.masking import MaskSet,propose
from nexvirome2.catalog import build,split
from nexvirome2.integrated_evidence import domains,link,joint,discovery


def test_reverse_domain_projection_and_mask_overlap(tmp_path):
    ref=tmp_path/'ref.fa'
    write_fasta(ref,{'A':'ACGT'*30})
    table(tmp_path/'orfs.tsv',[dict(orf_id='o',contig='A',start=12,end=72,strand='-')],
          ['orf_id','contig','start','end','strand'])
    table(tmp_path/'hits.tsv',[dict(orf_id='o',aa_start=2,aa_end=5,evidence_id='d',kind='domain',label='polymerase')],
          ['orf_id','aa_start','aa_end','evidence_id','kind','label'])
    domains({},ref,tmp_path/'orfs.tsv',tmp_path/'hits.tsv',tmp_path/'domains')
    assert (rows(tmp_path/'domains/annotations.tsv')[0]['start'],rows(tmp_path/'domains/annotations.tsv')[0]['end'])==('57','66')
    MaskSet.create(ref,[dict(accession='A',start=60,end=70)],dict(version='v1',rank='species',conditions={})).save(tmp_path/'m.json')
    link({},ref,tmp_path/'m.json',tmp_path/'domains/annotations.tsv',tmp_path/'linked')
    assert rows(tmp_path/'linked/links.tsv')[0]['masked_bases']=='6'


def test_joint_unknowns_and_overlaps(tmp_path):
    write_fasta(tmp_path/'ref.fa',{'A':'ACGT'*30})
    records=[dict(unit='u',accession='A',start=0,end=10,classification='ambiguous',assembly='chimeric'),
             dict(unit='u',accession='A',start=10,end=20,classification='unknown',assembly='consistent')]
    table(tmp_path/'calls.tsv',records,list(records[0]))
    joint({},tmp_path/'ref.fa',tmp_path/'calls.tsv',tmp_path/'joint')
    assert config(tmp_path/'joint/summary.json')['unresolved_bases']==10
    records[1]['start']=9
    table(tmp_path/'calls.tsv',records,list(records[0]))
    with pytest.raises(ValueError,match='Overlapping'):
        joint({},tmp_path/'ref.fa',tmp_path/'calls.tsv',tmp_path/'bad')


def test_catalog_orientation_and_transitive_holdout_groups(tmp_path):
    sequence='AAAACCGT'
    write_fasta(tmp_path/'a.fa',{'a':sequence})
    write_fasta(tmp_path/'b.fa',{'b':revcomp(sequence)})
    table(tmp_path/'in.tsv',[dict(source='meta',fasta='a.fa'),dict(source='mega',fasta='b.fa')],['source','fasta'])
    build({},tmp_path/'in.tsv',tmp_path/'catalog')
    assert config(tmp_path/'catalog/summary.json')['representatives']==1
    assert {r['orientation'] for r in rows(tmp_path/'catalog/lineage.tsv')}=={'+','-'}
    members=[dict(accession='A',cluster='c1',genome_group='g1'),dict(accession='B',cluster='c1',genome_group='g2'),
             dict(accession='C',cluster='c2',genome_group='g2'),dict(accession='D',cluster='c3',genome_group='g3')]
    table(tmp_path/'members.tsv',members,list(members[0]))
    split({},tmp_path/'members.tsv',tmp_path/'split')
    result={r['accession']:r for r in rows(tmp_path/'split/split.tsv')}
    assert result['A']['unit']==result['B']['unit']==result['C']['unit']
    assert result['A']['split']!=result['D']['split']


def test_unsearched_is_not_dark_discovery(tmp_path):
    profile=dict(orf_id='o',contig='c',sequence_related='False',structural_related='False',origin='uncertain')
    table(tmp_path/'profiles.tsv',[profile],list(profile))
    discovery({},tmp_path/'profiles.tsv',tmp_path/'unknown')
    assert rows(tmp_path/'unknown/cards.tsv')[0]['novelty']=='insufficient_assessment'
    table(tmp_path/'motif.tsv',[dict(orf_id='o',status='unsupported')],['orf_id','status'])
    table(tmp_path/'context.tsv',[dict(orf_id='o',status='supported')],['orf_id','status'])
    discovery(dict(sequence_assessed=True,fold_assessed=True),tmp_path/'profiles.tsv',tmp_path/'dark',
              tmp_path/'motif.tsv',tmp_path/'context.tsv')
    assert rows(tmp_path/'dark/cards.tsv')[0]['novelty']=='dark_context_candidate'
