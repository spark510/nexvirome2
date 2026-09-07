import pytest
from nexvirome2.common import write_fasta,dump,table,rows,config
from nexvirome2.experimental_design import multisample,rank_paths
from nexvirome2.correction.support import pair_evidence


def test_biological_samples_vary_but_technical_seeds_share_abundance(tmp_path):
    write_fasta(tmp_path/'a.fa',{'A':'ACGT'*40})
    write_fasta(tmp_path/'b.fa',{'B':'TGCA'*40})
    dump(tmp_path/'sources.json',{'sources':[dict(accession='A',fasta='a.fa'),dict(accession='B',fasta='b.fa')]})
    settings=dict(biological_samples=4,technical_seeds=[17,29],seed=17)
    multisample(settings,tmp_path/'sources.json',tmp_path/'independent')
    first=config(tmp_path/'independent/b0_s17.sources.json')['sources']
    repeat=config(tmp_path/'independent/b0_s29.sources.json')['sources']
    other=config(tmp_path/'independent/b1_s17.sources.json')['sources']
    assert first==repeat
    assert first[0]['coverage']/first[1]['coverage']!=other[0]['coverage']/other[1]['coverage']
    multisample({**settings,'mode':'correlated'},tmp_path/'sources.json',tmp_path/'correlated')
    for r in range(4):
        source=config(tmp_path/f'correlated/b{r}_s17.sources.json')['sources']
        assert source[0]['coverage']==source[1]['coverage']


def test_missing_path_feature_is_not_negative_evidence(tmp_path):
    table(tmp_path/'features.tsv',[dict(candidate='c',path='p',feature='coverage',value=.9),
                                  dict(candidate='c',path='p',feature='structure',value='')],
          ['candidate','path','feature','value'])
    rank_paths({'weights':{'coverage':1,'structure':1}},tmp_path/'features.tsv',tmp_path/'rank')
    row=rows(tmp_path/'rank/ranking.tsv')[0]
    assert float(row['score'])==.9
    assert row['missing']=='structure'
    assert row['decision']=='ranking_only'


def test_ambiguous_fragment_candidates_are_preserved(tmp_path):
    lines=[]
    for target in ('p0','p1'):
        lines.extend([f'q\t99\t{target}\t51\t0\t50M\t=\t201\t200\t'+('A'*50)+'\t*\tAS:i:0',
                      f'q\t147\t{target}\t201\t0\t50M\t=\t51\t-200\t'+('T'*50)+'\t*\tAS:i:0'])
    (tmp_path/'a.sam').write_text('\n'.join(lines)+'\n')
    counts,ambiguous=pair_evidence(tmp_path/'a.sam',{'p0':(100,200),'p1':(100,200)},audit_file=tmp_path/'audit.json')
    assert ambiguous==1 and sum(counts.values())==0
    evidence=config(tmp_path/'audit.json')['fragments'][0]
    assert evidence['status']=='ambiguous'
    assert {r['path'] for r in evidence['candidates']}=={'p0','p1'}
