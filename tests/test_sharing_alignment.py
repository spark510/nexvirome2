import pytest
from nexvirome2.common import write_fasta,table,rows
from nexvirome2.sharing import build
from nexvirome2.sharing_alignment import build as align


def test_circular_seed_events_split_at_origin(tmp_path):
    write_fasta(tmp_path/'ref.fa',{'A':'ACGTTGCA','B':'TGCAACGT'})
    table(tmp_path/'tax.tsv',[dict(accession='A',species='1'),dict(accession='B',species='2')],['accession','species'])
    table(tmp_path/'meta.tsv',[dict(accession='A',topology='circular'),dict(accession='B',topology='circular')],['accession','topology'])
    build({'k':5},tmp_path/'ref.fa',tmp_path/'tax.tsv',tmp_path/'atlas',tmp_path/'meta.tsv')
    evidence=rows(tmp_path/'atlas/atlas.tsv')
    wrap=[r for r in evidence if r['evidence_id']=='A:6:5']
    assert {(r['start'],r['end']) for r in wrap}=={('6','8'),('0','3')}
    assert all(r['taxa_count']=='2' for r in evidence)


def test_approximate_sharing_without_exact_full_length_seed(tmp_path):
    write_fasta(tmp_path/'ref.fa',{'A':'ACGTTACG','B':'ACGTGACG'})
    table(tmp_path/'tax.tsv',[dict(accession='A',species='1'),dict(accession='B',species='2')],['accession','species'])
    hit='A\tB\t87.5\t8\t1\t8\t1\t8\t20\tACGTTACG\tACGTGACG\n'
    (tmp_path/'hits.tsv').write_text(hit)
    align({'min_identity':80,'min_alignment_bp':8},tmp_path/'ref.fa',tmp_path/'tax.tsv',tmp_path/'atlas',tmp_path/'hits.tsv')
    result=rows(tmp_path/'atlas/atlas.tsv')
    assert {r['accession'] for r in result}=={'A','B'}
    assert all(r['taxa_count']=='2' and r['method']=='competitive_blast' for r in result)
    (tmp_path/'hits.tsv').write_text(hit.replace('ACGTGACG','ACGTAACG'))
    with pytest.raises(ValueError,match='differs from original'):
        align({'min_identity':80,'min_alignment_bp':8},tmp_path/'ref.fa',tmp_path/'tax.tsv',tmp_path/'bad',tmp_path/'hits.tsv')
