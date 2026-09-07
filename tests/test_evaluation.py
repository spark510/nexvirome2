import random
from nexvirome2.evaluation import classify, aligned_bases, evaluate
from nexvirome2.common import write_fasta
import json


def hit(source, sequence, start=0, subject_start=0, identity=100):
    return dict(qseqid='contig', sseqid=source, pident=identity, length=len(sequence),
                qstart=start+1, qend=start+len(sequence), sstart=subject_start+1,
                send=subject_start+len(sequence), bitscore=2*len(sequence), qseq=sequence, sseq=sequence)


def test_source_and_shared():
    seq = 'ACGT'*100
    assert classify(seq, [hit('A', seq)])['status'] == 'consistent'
    assert classify(seq, [hit('A', seq), hit('B', seq)])['status'] == 'ambiguous'
    assert classify(seq, [])['status'] == 'unaligned'


def test_true_chimera_and_breakpoint_interval():
    seq = 'ACGT'*150
    a, b = hit('A', seq[:350]), hit('B', seq[250:], 250)
    result = classify(seq, [a, b])
    assert result['status'] == 'chimeric'
    assert result['breakpoints'] == [{'left_source': 'A', 'right_source': 'B', 'start': 250, 'end': 350}]


def test_real_recombinant_in_source_wins():
    seq = 'ACGT'*150
    assert classify(seq, [hit('A', seq[:300]), hit('B', seq[300:], 300), hit('R', seq)])['status'] != 'chimeric'


def test_circular_origin_and_single_source_rearrangement():
    seq = 'ACGT'*100
    assert classify(seq, [hit('A', seq[:200], subject_start=800), hit('A', seq[200:], 200)])['status'] == 'consistent'


def test_reverse_complement_coordinates():
    h = hit('A', 'ACGT'*50)
    h.update(sstart=200, send=1)
    bases = list(aligned_bases(h))
    assert bases[0][:2] == (0, 199)
    assert bases[-1][:2] == (199, 0)
    assert classify('ACGT'*50, [h])['status'] == 'consistent'


def test_short_or_low_identity_is_not_chimera():
    seq = 'ACGT'*100
    assert classify(seq, [hit('A', seq[:100]), hit('B', seq[100:], 100, identity=98)])['status'] == 'unaligned'


def test_evaluate_files_and_filtered_recovery(tmp_path):
    seq = 'ACGT'*150
    q, s, tsv = tmp_path/'q.fa', tmp_path/'s.fa', tmp_path/'hits.tsv'
    write_fasta(q, {'contig': seq})
    write_fasta(s, {'A': seq[:300], 'B': seq[300:]})
    from nexvirome2.alignment import FIELDS
    hs = [hit('A', seq[:300]), hit('B', seq[300:], 300)]
    tsv.write_text('\n'.join('\t'.join(str(h[f]) for f in FIELDS.split()) for h in hs)+'\n')
    evaluate({}, q, s, tmp_path/'out', tsv)
    metrics = json.loads((tmp_path/'out/metrics.json').read_text())
    assert metrics['0']['wrong_connections'] == 1
    assert metrics['500']['genome_recovery'] == 1
