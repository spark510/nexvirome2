"""External command fixtures test wiring, not biological tool performance."""
import gzip
import json
from pathlib import Path
import random
import pytest

from nexvirome2.common import write_fasta, fasta, dump


def test_art_coordinates_circular_and_neutral_names(tmp_path, monkeypatch):
    import nexvirome2.simulation as sim
    source = tmp_path/'source.fa'
    write_fasta(source, {'NC_123456.1':'ACGT'*250})
    spec = tmp_path/'sources.json'
    dump(spec, {'sources':[{'accession':'NC_123456.1','fasta':str(source),'coverage':50,'topology':'circular'}]})
    def art(args, out, manifest, name):
        prefix = str(args[args.index('-o')+1])
        seq, quality = 'A'*150, 'I'*150
        Path(prefix+'1.fq').write_text(f'@template-1/1\n{seq}\n+\n{quality}\n')
        Path(prefix+'2.fq').write_text(f'@template-1/2\n{seq}\n+\n{quality}\n')
        Path(prefix+'.sam').write_text('template-1\t99\ttemplate\t901\t60\t150M\t=\t1051\t300\t*\t*\n'
                                      'template-1\t147\ttemplate\t1051\t60\t150M\t=\t901\t-300\t*\t*\n')
    monkeypatch.setattr(sim, 'command', art)
    monkeypatch.setattr(sim, 'version', lambda *a: None)
    sim.simulate({}, spec, tmp_path/'sim')
    with gzip.open(tmp_path/'sim/reads_R1.fastq.gz','rt') as handle:
        assert handle.readline().startswith('@fragment_000000000000/1')
    from nexvirome2.common import rows
    truth = rows(tmp_path/'sim/truth/read_origin.tsv')
    assert truth[0]['start'] == '900' and truth[0]['wraps_origin'] == 'True'
    assert truth[1]['start'] == '50'


@pytest.mark.parametrize('proposal_mode,expected_parts',[('baseline',2),('missing',1),('matching',2)])
def test_correction_end_to_end_without_truth(tmp_path, monkeypatch,proposal_mode,expected_parts):
    import nexvirome2.correction as correction
    rng = random.Random(4)
    seq = lambda n: ''.join(rng.choices('ACGT', k=n))
    segments = {str(i):seq(n) for i,n in [(1,100),(2,50),(3,100),(4,100),(5,100)]}
    gfa = tmp_path/'graph.gfa'
    gfa.write_text(''.join(f'S\t{k}\t{v}\n' for k,v in segments.items())+
                   ''.join(f'L\t{a}\t+\t{b}\t+\t0M\n' for a,b in [(1,2),(5,2),(2,3),(2,4)]))
    contigs = tmp_path/'contigs.fa'
    write_fasta(contigs, {'contig':segments['1']+segments['2']+segments['3']})
    paths = tmp_path/'paths'
    paths.write_text('contig\n1+,2+,3+\n')
    r1, r2 = tmp_path/'r1.fq', tmp_path/'r2.fq'
    r1.write_text(''); r2.write_text('')
    def external(args, out, manifest, name, stdout=None):
        if args[0] == 'bowtie2':
            target = Path(args[args.index('-S')+1])
            lines = []
            for i in range(3):
                start, end = 10+i, 160+i
                lines.append(f'r{i}\t99\tpath_1\t{start+1}\t42\t50M\t=\t{end+1}\t200\t*\t*\tAS:i:0\n')
                lines.append(f'r{i}\t147\tpath_1\t{end+1}\t42\t50M\t=\t{start+1}\t-200\t*\t*\tAS:i:0\n')
            target.write_text(''.join(lines))
        if stdout:
            Path(stdout).write_text('')
    monkeypatch.setattr(correction, 'command', external)
    monkeypatch.setattr(correction, 'version', lambda *a: None)
    from nexvirome2.common import table
    settings,proposals = {},[]
    if proposal_mode!='baseline':
        settings['proposal_policy']='require_any'
        proposals=[tmp_path/'proposals.tsv']
        table(proposals[0],[dict(contig='contig',start=140,end=160,method='nmf',score=.9)] if proposal_mode=='matching' else [],
              ['contig','start','end','method','score'])
    correction.correct(settings, contigs, gfa, paths, r1, r2, tmp_path/'correction',proposals)
    out = tmp_path/'correction'
    result = fasta(out/'corrected.fasta')
    assert list(map(len, result.values())) == ([150,100] if expected_parts==2 else [250])
    assert ''.join(result.values()) == fasta(contigs)['contig']
    assert len(fasta(out/'random_control.fasta')) == expected_parts
    assert json.loads((out/'summary.json').read_text())['truth_used'] is False


def test_design_counts_and_independence(tmp_path):
    from nexvirome2.benchmark import prepare
    panel = tmp_path/'panel.json'
    dump(panel, {'pairs':[{'a':'A','b':'B'},{'a':'C','b':'D'},{'a':'E','b':'F'}]})
    prepare({}, panel, tmp_path/'design')
    experiments = json.loads((tmp_path/'design/experiments.json').read_text())['experiments']
    assert len(experiments) == 27
    mixed = json.loads((tmp_path/'design/p0_c50_r1_s20260905_mix.sources.json').read_text())
    single_b = json.loads((tmp_path/'design/p0_c50_r1_s20260905_single_b.sources.json').read_text())
    assert mixed['sources'][1] == single_b['sources'][0]
    assert single_b['sources'][0]['seed_offset'] == 1
    with pytest.raises(ValueError, match='overlaps'):
        prepare({'partition':'evaluation','exclude_accessions':['A']}, panel, tmp_path/'heldout')


def test_report_zero_baseline_is_not_success(tmp_path):
    from nexvirome2.report import report
    for method in ['metaspades','megahit','corrected','random_control']:
        directory = tmp_path/'runs'/'sample'/'evaluation'/method
        directory.mkdir(parents=True)
        dump(directory/'metrics.json', {'0':{'contigs':1,'wrong_connections':0,'genome_recovery':1.0}})
    report({}, tmp_path/'runs', tmp_path/'report')
    data = json.loads((tmp_path/'report/targets.json').read_text())
    assert data['comparisons'][0]['connection_reduction'] is None
    assert not data['comparisons'][0]['target_met']
    assert (tmp_path/'report/benchmark.png').stat().st_size > 0
