import json
import pytest
from nexvirome2.common import write_fasta, table, config, dump, sha, rows
from nexvirome2.masking import MaskSet, import_mask, export_mask, propose
from nexvirome2.sharing import build
from nexvirome2.classification import classify, evaluate
from nexvirome2.mask_benchmark import run
from nexvirome2.decisions import freeze, select


@pytest.fixture
def inputs(tmp_path):
    ref, tax = tmp_path/'ref.fa', tmp_path/'tax.tsv'
    write_fasta(ref, {'A.1':'ACGTACGATCGGATCCGATAGCTAGCTA','B.1':'ACGTACGATCGGATCCGATTTTTTTTT'})
    table(tax,[dict(accession='A.1',species='1',genus='10'),dict(accession='B.1',species='2',genus='10')],
          ['accession','species','genus'])
    return ref,tax


POLICY = dict(version='test-v1',rank='species',conditions={'k':5})


@pytest.mark.parametrize('mode',['hard','soft'])
def test_mask_round_trip(inputs,tmp_path,mode):
    ref,_ = inputs
    before = sha(ref)
    mask = MaskSet.create(ref,[dict(accession='A.1',start=2,end=7),dict(accession='A.1',start=5,end=10)],POLICY)
    mask.save(tmp_path/'mask.json')
    export_mask({'mode':mode},ref,tmp_path/'mask.json',tmp_path/'export')
    import_mask({'format':mode,'policy':POLICY},ref,tmp_path/'export/reference.fasta',tmp_path/'import')
    back = MaskSet.load(tmp_path/'import/mask.json',ref)
    assert back.merged('A.1')==[(2,10)]
    assert back.render(ref,'unmasked')==mask.render(ref,'unmasked')
    assert sha(ref)==before


def test_mask_wrong_reference_and_coordinates(inputs,tmp_path):
    ref,_ = inputs
    with pytest.raises(ValueError,match='coordinate'):
        MaskSet.create(ref,[dict(accession='A.1',start=9,end=99)],POLICY)
    mask = MaskSet.create(ref,[],POLICY)
    mask.save(tmp_path/'m.json')
    ref.write_text('>A.1\nAAA\n')
    with pytest.raises(ValueError,match='checksum'):
        MaskSet.load(tmp_path/'m.json',ref)


def test_distinct_taxon_sharing_and_random_amount(inputs,tmp_path):
    ref,tax=inputs
    build({'k':5,'rank':'species'},ref,tax,tmp_path/'species')
    build({'k':5,'rank':'genus'},ref,tax,tmp_path/'genus')
    assert any(int(r['taxa_count'])==2 for r in rows(tmp_path/'species/atlas.tsv'))
    assert all(int(r['taxa_count'])==1 for r in rows(tmp_path/'genus/atlas.tsv'))
    propose({'policy':POLICY},ref,tmp_path/'species/atlas.tsv',tmp_path/'proposal')
    actual=MaskSet.load(tmp_path/'proposal/mask.json',ref)
    random=MaskSet.load(tmp_path/'proposal/random_mask.json',ref)
    for acc in actual.lengths:
        assert sorted(b-a for a,b in actual.merged(acc))==sorted(b-a for a,b in random.merged(acc))


def test_mask_benchmark_keeps_all_queries_and_truth_denominator(inputs,tmp_path):
    ref,tax=inputs
    query=tmp_path/'queries.fa'
    write_fasta(query,{'q1':'ACGTACGATCGGATCCGAT','q2':'TTTTTTTT','q3':'NNNNNNNN'})
    truth=tmp_path/'truth.tsv'
    table(truth,[dict(query='q1',taxon='1',unit='u1'),dict(query='q2',taxon='2',unit='u2'),
                 dict(query='q3',taxon='1',unit='u3')],['query','taxon','unit'])
    mask=MaskSet.create(ref,[dict(accession='B.1',start=0,end=27)],POLICY)
    mask.save(tmp_path/'mask.json')
    dump(tmp_path/'arms.json',{'unmasked':None,'masked':'mask.json'})
    run({'k':5,'rank':'species'},query,ref,tax,truth,tmp_path/'arms.json',tmp_path/'bench')
    baseline=config(tmp_path/'bench/unmasked/evaluation/metrics.json')
    assert baseline['total']==3
    assert baseline['ambiguous']==pytest.approx(1/3)
    assert baseline['unclassified']==pytest.approx(1/3)
    assert baseline['correct']==pytest.approx(1/3)
    info=config(tmp_path/'bench/comparison.json')
    assert info['biological_performance_validated'] is False
    assert info['queries_sha256']==sha(query)


def test_duplicate_reference_taxon_and_unknown_do_not_force_assignment(inputs,tmp_path):
    ref,tax=inputs
    write_fasta(tmp_path/'q.fa',{'q':'ACGTACGATCGGATCCGAT'})
    classify({'rank':'genus','k':5},tmp_path/'q.fa',ref,tax,tmp_path/'g')
    assert rows(tmp_path/'g/predictions.tsv')[0]['status']=='assigned'
    table(tax,[dict(accession='A.1',species='1'),dict(accession='B.1',species='')],['accession','species'])
    classify({'rank':'species','k':5},tmp_path/'q.fa',ref,tax,tmp_path/'s')
    assert rows(tmp_path/'s/predictions.tsv')[0]['status']=='ambiguous'


def test_adoption_gate_fixture_leakage_and_harm(tmp_path):
    params={'modules':{'mask':{'max_recovery_loss':.02,'min_relative_error_reduction':.5,
                             'implementation':{'command':['mask','export'],'settings':{'mode':'hard'}}}}}
    table(tmp_path/'dev.tsv',[{'unit':'dev'}],['unit'])
    freeze({'parameters':params},tmp_path/'dev.tsv',tmp_path/'frozen')
    records=[dict(module='mask',unit=f'u{i}',baseline_error=.1,error=.01,baseline_recovery=.9,
                  recovery=.9,resource_ratio=1.2) for i in range(4)]
    fields=list(records[0])
    table(tmp_path/'e.tsv',records,fields)
    for kind,expected in [('fixture','research'),('natural_heldout','default')]:
        select({'parameters':params,'evidence_kind':kind},tmp_path/'e.tsv',tmp_path/'frozen/frozen.json',tmp_path/kind)
        assert config(tmp_path/kind/'decisions.json')['decisions'][0]['status']==expected
    records[0]['unit']='dev'
    table(tmp_path/'e.tsv',records,fields)
    with pytest.raises(ValueError,match='leakage'):
        select({'parameters':params},tmp_path/'e.tsv',tmp_path/'frozen/frozen.json',tmp_path/'bad')
    for i,r in enumerate(records):
        r.update(unit=f'u{i}',recovery=.5)
    table(tmp_path/'e.tsv',records,fields)
    select({'parameters':params,'evidence_kind':'natural_heldout'},tmp_path/'e.tsv',tmp_path/'frozen/frozen.json',tmp_path/'harm')
    assert config(tmp_path/'harm/decisions.json')['decisions'][0]['status']=='excluded'


def test_taxdump_merged_deleted_and_missing_rank(tmp_path):
    from nexvirome2.taxonomy import build,TaxonomyResolver
    (tmp_path/'nodes.dmp').write_text('1 | 1 | no rank |\n10 | 1 | genus |\n11 | 10 | species |\n')
    (tmp_path/'merged.dmp').write_text('12 | 11 |\n')
    (tmp_path/'del.dmp').write_text('99 |\n')
    table(tmp_path/'map.tsv',[dict(accession='A',taxid='12'),dict(accession='B',taxid='99')],['accession','taxid'])
    build({},tmp_path/'map.tsv',tmp_path/'nodes.dmp',tmp_path/'tax',tmp_path/'merged.dmp',tmp_path/'del.dmp')
    resolver=TaxonomyResolver.load(tmp_path/'tax/taxonomy.tsv')
    assert resolver.taxon('A','species')=='11'
    assert resolver.taxon('A','genus')=='10'
    assert resolver.taxon('A','family') is None
    assert resolver.taxon('B','species') is None
