from nexvirome2.graph import Graph, candidates, contig_paths
from nexvirome2.correction import decision, pair_evidence, split_sequences
import pytest


def test_gfa_orientation_and_overlap(tmp_path):
    gfa = tmp_path/'g.gfa'
    gfa.write_text('S\t1\tAACCGG\nS\t2\tGGTTAA\nL\t1\t+\t2\t+\t2M\n')
    graph = Graph.load(gfa)
    assert graph.spell(['1+', '2+']) == ('AACCGGTTAA', [0, 4])
    assert graph.spell(['2-', '1-'])[0] == 'TTAACCGGTT'


def test_invalid_overlap_rejected(tmp_path):
    gfa = tmp_path/'g.gfa'
    gfa.write_text('S\t1\tAAAA\nS\t2\tCCCC\nL\t1\t+\t2\t+\t2M\n')
    with pytest.raises(ValueError, match='overlap'):
        Graph.load(gfa)


def test_repeat_candidate_and_paths(tmp_path):
    graph = Graph()
    graph.sequences = {'1':'A'*100, '2':'C'*50, '3':'G'*100, '4':'T'*100, '5':'AC'*50}
    from nexvirome2.graph import flip
    for a,b in [('1+','2+'), ('5+','2+'), ('2+','3+'), ('2+','4+')]:
        graph.links[a][b] = 0
        graph.links[flip(b)][flip(a)] = 0
    paths = tmp_path/'contigs.paths'
    paths.write_text('contig\n1+,2+,3+\ncontig\'\n3-,2-,1-\n')
    sequence = graph.spell(['1+','2+','3+'])[0]
    found = candidates(graph, {'contig': sequence}, contig_paths(paths))
    assert found[0]['cut'] == 150
    assert found[0]['alternatives'] == [['1+','2+','4+']]


def test_no_support_is_unresolved():
    assert decision(0, 0) == 'unresolved'
    assert decision(0, 2) == 'unresolved'
    assert decision(0, 3) == 'split'
    assert decision(10, 3) == 'retain'


def sam_pair(name, template, start=10, score=0, secondary=False):
    a = 99 + (256 if secondary else 0)
    b = 147 + (256 if secondary else 0)
    return (f'{name}\t{a}\t{template}\t{start+1}\t42\t50M\t=\t{start+201}\t250\t*\t*\tAS:i:{score}\n'
            f'{name}\t{b}\t{template}\t{start+201}\t42\t50M\t=\t{start+1}\t-250\t*\t*\tAS:i:{score}\n')


def test_multimapping_and_duplicate_fragments(tmp_path):
    sam = tmp_path/'x.sam'
    sam.write_text(sam_pair('unique','alt') + sam_pair('duplicate','alt') +
                   sam_pair('other','alt', 20) + sam_pair('tie','alt', 30) + sam_pair('tie','current',30, secondary=True))
    support, ambiguous = pair_evidence(sam, {'current':(100,200),'alt':(100,200)})
    assert support == {'current':0,'alt':2}
    assert ambiguous == 1


def test_repeat_only_mapping_not_counted(tmp_path):
    sam = tmp_path/'x.sam'
    sam.write_text(sam_pair('repeat','p'))
    assert pair_evidence(sam, {'p':(0,500)})[0]['p'] == 0


def test_split_preserves_all_bases_and_coordinates():
    sequences, mapping = split_sequences({'c':'ACGTACGT'}, {'c':[2,6]})
    assert ''.join(sequences.values()) == 'ACGTACGT'
    assert [(m['start'],m['end']) for m in mapping] == [(0,2),(2,6),(6,8)]
