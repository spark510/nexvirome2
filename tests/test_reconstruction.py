import json
import pytest
from nexvirome2.common import fasta
from nexvirome2.graph import Graph, flip
from nexvirome2.reconstruction import reconstruct, bounded_paths, coordinate_rows


class MappingFixture:
    def __init__(self, tied=False, count=3, verification=True):
        self.tied, self.count, self.verification = tied, count, verification

    def version(self, *args):
        pass

    def __call__(self, args, *rest):
        if args[0] != 'bowtie2':
            return
        sam = args[args.index('-S')+1]
        verifying = sam.parent.name == 'verification'
        if verifying and not self.verification:
            sam.write_text('')
            return
        lines = []
        for i in range(self.count):
            for target in (['c'] if verifying else ['path_1','path_0'] if self.tied else ['path_1']):
                for flag, pos, mate, length in [(99,51+i,191+i,190),(147,191+i,51+i,-190)]:
                    lines.append(f'q{i}\t{flag}\t{target}\t{pos}\t42\t50M\t=\t{mate}\t{length}\t*\t*\tAS:i:0\n')
        sam.write_text(''.join(lines))


@pytest.fixture
def inputs(tmp_path):
    gfa = tmp_path/'g.gfa'
    gfa.write_text('S\t1\t'+'A'*100+'\nS\t2\t'+'C'*60+'\nS\t3\t'+'T'*80+'\nS\t4\t'+'G'*100+
                   '\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n')
    contigs = tmp_path/'c.fa'
    contigs.write_text('>c\n'+'A'*100+'C'*60+'G'*100+'\n')
    paths = tmp_path/'c.paths'
    paths.write_text('c\n1+,2+,4+\n')
    reads = tmp_path/'r.fastq'
    reads.write_text('@q\nA\n+\nI\n')
    return contigs,gfa,paths,reads,reads


@pytest.mark.parametrize('apply,tied,count,changed',[(True,False,3,True),(False,False,3,False),(True,True,3,False),(True,False,2,False)])
def test_real_graph_spelling_with_mocked_mapping(tmp_path,inputs,apply,tied,count,changed):
    out = tmp_path/'run'
    reconstruct({'apply':apply},*inputs,out,runner=MappingFixture(tied,count))
    result = fasta(out/'reconstructed.fasta')['c']
    assert result == 'A'*100+('T'*80 if changed else 'C'*60)+'G'*100
    summary = json.loads((out/'summary.json').read_text())
    assert summary['applied'] == int(changed)
    assert json.loads((out/'manifest.json').read_text())['status']=='complete'
    if changed:
        from nexvirome2.common import rows
        mapping = rows(out/'graph_coordinate_map.tsv')
        assert [(int(r['output_start']),int(r['output_end'])) for r in mapping]==[(0,100),(100,180),(180,280)]
        assert mapping[1]['original_start']==''
        assert mapping[2]['original_start']=='160'
        assert json.loads((out/'output_paths.json').read_text())['c']==['1+','3+','4+']


def test_cycle_or_truncated_enumeration_is_not_actionable(inputs):
    graph = Graph.load(inputs[1])
    assert bounded_paths(graph,'1+','4+',2,32,100) is None
    assert bounded_paths(graph,'1+','4+',8,1,100) is None
    graph.links['3+']['1+']=0
    assert bounded_paths(graph,'1+','4+',8,32,100) is None


def test_reverse_overlap_coordinate_map(tmp_path):
    gfa=tmp_path/'g.gfa'
    gfa.write_text('S\t1\tAACCGG\nS\t2\tGGTTAA\nL\t1\t+\t2\t+\t2M\n')
    graph=Graph.load(gfa)
    path=['2-','1-']
    rows=coordinate_rows(graph,'c',path,path)
    assert [(r['output_start'],r['output_end']) for r in rows]==[(0,6),(6,10)]
    assert rows[1]['node_start']==0 and rows[1]['node_end']==4
    rebuilt=''
    from nexvirome2.common import revcomp
    for row in rows:
        sequence=graph.sequences[row['graph_node']][row['node_start']:row['node_end']]
        rebuilt+=revcomp(sequence) if row['strand']=='-' else sequence
    assert rebuilt==graph.spell(path)[0]


@pytest.mark.parametrize('path_text', ['c\n1+,3+,4+\n','c\n1+,999+,4+\n','c\n1+;4+\n','c\n'])
def test_mismatched_path_preserves_original(tmp_path,inputs,path_text):
    inputs[2].write_text(path_text)
    out=tmp_path/'out'
    reconstruct({'apply':True},*inputs,out,runner=MappingFixture())
    assert fasta(out/'reconstructed.fasta')==fasta(inputs[0])
    assert json.loads((out/'candidates.json').read_text())['skipped']
    assert json.loads((out/'output_paths.json').read_text())['c'] is None
    assert json.loads((out/'output_path_validation.json').read_text())['c']['valid'] is False


def test_failed_remapping_rolls_back_output(tmp_path,inputs):
    out=tmp_path/'out'
    reconstruct({'apply':True},*inputs,out,runner=MappingFixture(verification=False))
    assert fasta(out/'reconstructed.fasta')==fasta(inputs[0])
    assert json.loads((out/'decisions.json').read_text())[0]['reason']=='batch_remapping_failed'


def test_multiple_rewrites_are_withheld(tmp_path,inputs,monkeypatch):
    import nexvirome2.reconstruction as module
    original=module.discover
    def duplicates(*args):
        candidates,skipped=original(*args)
        return candidates+candidates,skipped
    monkeypatch.setattr(module,'discover',duplicates)
    out=tmp_path/'out'
    reconstruct({'apply':True},*inputs,out,runner=MappingFixture())
    assert fasta(out/'reconstructed.fasta')==fasta(inputs[0])
    assert all(r['reason']=='multiple_rewrites_require_joint_validation'
               for r in json.loads((out/'decisions.json').read_text()))


def test_cli_accepts_reconstruction_without_truth_input():
    from nexvirome2.cli import parser
    args=parser().parse_args(['reconstruct','run','--config','x','--output','o',
                             '--contigs','c','--graph','g','--paths','p','--r1','a','--r2','b'])
    assert args.command=='reconstruct' and not hasattr(args,'truth')


def test_retained_repeated_node_maps_each_occurrence():
    graph=Graph()
    graph.sequences={'1':'AAAA'}
    graph.links['1+']['1+']=0
    path=['1+','1+']
    rows=coordinate_rows(graph,'c',path,path)
    assert [r['original_start'] for r in rows]==[0,4]


def test_direct_edge_can_be_replaced_by_supported_detour(tmp_path,inputs):
    inputs[0].write_text('>c\n'+'A'*100+'G'*100+'\n')
    inputs[2].write_text('c\n1+,4+\n')
    inputs[1].write_text('S\t1\t'+'A'*100+'\nS\t3\t'+'T'*80+'\nS\t4\t'+'G'*100+
                        '\nL\t1\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n')
    out=tmp_path/'out'
    reconstruct({'apply':True},*inputs,out,runner=MappingFixture())
    assert fasta(out/'reconstructed.fasta')['c']=='A'*100+'T'*80+'G'*100
    path=json.loads((out/'output_paths.json').read_text())['c']
    assert path==['1+','3+','4+']
    assert Graph.load(inputs[1]).spell(path)[0]==fasta(out/'reconstructed.fasta')['c']


def test_direct_edge_candidates_in_reverse_orientation(inputs):
    from nexvirome2.reconstruction import discover
    graph=Graph.load(inputs[1])
    graph.links['1+']['4+']=0
    graph.links['4-']['1-']=0
    path=['4-','1-']
    candidates,_=discover(graph,{'c':graph.spell(path)[0]},{'c':path})
    assert len(candidates)==1
    assert ['4-','3-','1-'] in candidates[0]['alternatives']
